#!/usr/bin/env python3
"""Kirsch plate-with-hole benchmark: HPC scalability sweep across four solve
pipelines (legacy Fortran, future Gauss+direct, future matrix-free, pymfiga
matrix-free), plus an Abaqus deck generator for external comparison.

Same problem as examples/future/12_plate_with_hole.ipynb / 16 / 17: quarter
annulus (Rin=1, Rex=2), exact Kirsch (1898) traction on the outer arc, symmetry
BCs on the two straight edges. This script packages that already-validated
geometry/physics into a standalone, parametrized, fault-tolerant tool meant for
a compute cluster (large `nbel`, long runtimes, partial failures).

Usage
-----
    # one (degree, nbel, pipeline) case -- normally only invoked by `sweep`
    python kirsch_sweep.py worker --pipeline future-mf --degree 3 --nbel 64 \\
        --n-cpus 1 --result-file /tmp/result.json

    # the full grid, one subprocess per case, resumable, incremental CSV
    python kirsch_sweep.py sweep --degrees 2 3 4 --nbels 1 2 4 8 16 32 \\
        --n-cpus 1 --results-csv results/sweep_results.csv

    # Abaqus decks only, degree-1, one per nbel (independent of the degree sweep)
    python kirsch_sweep.py abaqus --nbels 1 2 4 8 16 32 --out-dir results/abaqus

    # ... or write the decks AND launch Abaqus (ABAQUS_LAUNCH_CMD below) AND
    # the odb extraction (needs `abaqus python`, i.e. Abaqus itself, on PATH)
    # for each one, appending each nbel's H1 error to --results-csv:
    python kirsch_sweep.py abaqus --nbels 1 2 4 8 16 32 --out-dir results/abaqus \\
        --run --n-cpus 4 --results-csv results/sweep_results.csv

    # or, if you ran Abaqus yourself and already have a CSV from
    # `abaqus python extract_abaqus_results.py <job>.odb <job>_results.csv`
    # (next to this file -- needs Abaqus's own Python, not this one):
    python kirsch_sweep.py abaqus-results --csv-in job_results.csv --nbel 8 \\
        --results-csv results/sweep_results.csv

    # the 3 plots (H1 error, wall time, peak memory), one curve per (pipeline, degree)
    # -- includes "abaqus" rows from --run/abaqus-results alongside the 4 pipelines
    python kirsch_sweep.py plot --results-csv results/sweep_results.csv --out-dir results/plots

Fill in ABAQUS_LAUNCH_CMD below before using `abaqus --run` -- see that
constant's own comment for the one constraint on it (the job name).
"""
import sys

# Thread pinning MUST happen before numpy/scipy/yeti_iga are imported (these
# libraries read the thread count once, at first use). Each `worker` runs as
# its own fresh interpreter (see the "why subprocess-per-case" note in
# `cmd_sweep`), so setting this at module import time is enough -- no need to
# thread --n-cpus through as an env var to a already-running process.
if "--n-cpus" in sys.argv:
    _n_cpus = sys.argv[sys.argv.index("--n-cpus") + 1]
else:
    _n_cpus = "1"
import os
for _env in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_env] = _n_cpus

import argparse
import json
import logging
import re
import resource
import subprocess
import time
import traceback
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as spla

from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager, Patch as FPatch,
    GlobalDOFManager, PatchDOFManager, Material, PlaneStress,
    IGABasis1D, PatchIntegrator, Traction, ScalarLocalOperator,
    WQMatrixFreeStiffness,
)
from yeti_iga.future.fast_diagonalization import FastDiagonalizationPreconditioner
from yeti_iga.future.matrix_free_solver import solve as mf_solve

from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel

from yeti_iga import IgaModel, Patch as LPatch, ElasticMaterial

# ---------------------------------------------------------------------------
# Fill this in before running `abaqus --run`: it IS invoked (via shell,
# `cmd_abaqus`'s `run_abaqus_job`), formatted with `.format(nbel=nbel,
# n_cpus=n_cpus)` (n_cpus from `abaqus --n-cpus`, default 1 -- same flag as
# `worker`/`sweep`).
# The `job=` name MUST stay `kirsch_nbel{nbel}` (matching the deck filename
# `write_abaqus_deck` writes, `kirsch_nbel<N>.inp`) -- Abaqus resolves the
# input deck from the job name when `input=` is omitted, so a different job
# name here means Abaqus won't find the deck. `interactive` is required too:
# without it `abaqus job=...` returns immediately instead of blocking until
# the analysis finishes, and this script runs the odb extraction right after.
# `cpus={n_cpus}` is what actually tells Abaqus how many cores to use --
# unrelated to this script's own OMP_NUM_THREADS/etc. pinning (top of file),
# which only affects future/pymfiga's own code, not the external Abaqus process.
ABAQUS_LAUNCH_CMD = ("singularity exec --nv /opt/img/abaqus2019_intel2019u6_centos7_v1.sif "
                    "abaqus job=kirsch_nbel{nbel} cpus={n_cpus} interactive ask_delete=OFF")

# ---------------------------------------------------------------------------
# Kirsch benchmark parameters (examples/future/12_plate_with_hole.ipynb)
Rin_DEFAULT, Rex_DEFAULT, Tx_DEFAULT = 1.0, 2.0, 1.0
E_DEFAULT, NU_DEFAULT = 1e3, 0.3

# pymfiga's NURBS arc builder asserts degree > 1 (a degree-1 rational curve is
# a polyline, not a circle) -- the main sweep only accepts degree >= this.
BASE_DEGREE = 2

PIPELINES = ["legacy", "future-gauss", "future-mf", "pymfiga-mf"]


# ---------------------------------------------------------------------------
# Kirsch physics: exact traction, exact displacement + gradient, error norms.
# Verbatim from examples/future/12_plate_with_hole.ipynb / 16, the validated
# accuracy benchmark this script reuses rather than reinvents -- except
# kirsch_displacement_gradient, switched to its closed form here (see that
# function's docstring): the notebooks only sweep degree<=5, but this script
# is meant to go higher, where the notebooks' finite-difference version would
# floor the H1 error instead of showing true convergence.
# ---------------------------------------------------------------------------

class KirschTraction(Traction):
    """Exact Kirsch traction at the outer arc of the quarter annulus."""

    def __init__(self, Rin, Tx=1.0):
        super().__init__()
        self.Rin = Rin
        self.Tx = Tx

    def evaluate(self, pt):
        x, y = pt[0], pt[1]
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        R2 = (self.Rin / r) ** 2
        R4 = (self.Rin / r) ** 4
        g1 = (self.Tx / 2) * (
            2 * np.cos(theta)
            - R2 * (2 * np.cos(theta) + 3 * np.cos(3 * theta))
            + 3 * R4 * np.cos(3 * theta)
        )
        g2 = (3 * self.Tx / 2) * np.sin(3 * theta) * (R4 - R2)
        return np.array([g1, g2])


def kirsch_surface_force(Rin, Tx):
    """pymfiga's `assemble_surface_force` convention: a callable receiving
    args['position'] of shape (2, n_pts) and returning a (2, n_pts) force
    array in the SAME global-Cartesian convention as future's Traction
    (confirmed in 17_future_vs_pymfiga_matrix_free.ipynb cell 8: a constant
    force [0,1] matches bit-for-bit between the two conventions)."""

    def fn(args):
        position = np.asarray(args["position"])
        x, y = position[0], position[1]
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        R2 = (Rin / r) ** 2
        R4 = (Rin / r) ** 4
        g1 = (Tx / 2) * (
            2 * np.cos(theta)
            - R2 * (2 * np.cos(theta) + 3 * np.cos(3 * theta))
            + 3 * R4 * np.cos(3 * theta)
        )
        g2 = (3 * Tx / 2) * np.sin(3 * theta) * (R4 - R2)
        return np.vstack([g1, g2])

    return fn


def kirsch_displacement(pt, Rin, Tx, E, nu):
    x, y = pt[0], pt[1]
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    a2 = Rin**2
    u_x = Tx / (2 * E) * (2 * (r + 2 * a2 / r) * np.cos(theta)
                          + (1 + nu) * a2 / r * (1 - a2 / r**2) * np.cos(3 * theta))
    u_y = Tx / (2 * E) * (-2 * (nu * r + (1 - nu) * a2 / r) * np.sin(theta)
                          + (1 + nu) * a2 / r * (1 - a2 / r**2) * np.sin(3 * theta))
    return np.array([u_x, u_y])


def kirsch_displacement_gradient(pt, Rin, Tx, E, nu):
    """Analytic gradient of kirsch_displacement, via the chain rule on (r, theta).

    Replaces an earlier central-difference version: at degree>=~8 the IGA
    solution converges (for this smooth exact field) to H1 errors below the
    finite-difference noise floor (~eps^(2/3) ~ 3e-11 at the scheme's
    optimal step, regardless of h), which plateaued the error curve instead
    of showing the true convergence. Verified against the old FD
    implementation (agreement ~1e-12) and against a 5-point stencil at
    h=1e-3 (agreement ~5e-14).

    u_x = C*(A(r) cos(theta) + B(r) cos(3 theta))
    u_y = C*(D(r) sin(theta) + B(r) sin(3 theta))
    with C = Tx/(2E), a2 = Rin**2,
    A(r) = 2r + 4 a2/r, D(r) = -2 nu r - 2(1-nu) a2/r,
    B(r) = (1+nu)(a2/r - a2**2/r**3).
    """
    x, y = pt[0], pt[1]
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    a2 = Rin**2
    C = Tx / (2 * E)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    cos_3t, sin_3t = np.cos(3 * theta), np.sin(3 * theta)

    A = 2 * r + 4 * a2 / r
    Ap = 2 - 4 * a2 / r**2
    B = (1 + nu) * (a2 / r - a2**2 / r**3)
    Bp = (1 + nu) * (3 * a2**2 / r**4 - a2 / r**2)
    D = -2 * nu * r - 2 * (1 - nu) * a2 / r
    Dp = -2 * nu + 2 * (1 - nu) * a2 / r**2

    dux_dr = C * (Ap * cos_t + Bp * cos_3t)
    dux_dth = C * (-A * sin_t - 3 * B * sin_3t)
    duy_dr = C * (Dp * sin_t + Bp * sin_3t)
    duy_dth = C * (D * cos_t + 3 * B * cos_3t)

    # chain rule: df/dx = f_r cos(theta) - f_theta sin(theta)/r
    #             df/dy = f_r sin(theta) + f_theta cos(theta)/r
    dux_dx = dux_dr * cos_t - dux_dth * sin_t / r
    dux_dy = dux_dr * sin_t + dux_dth * cos_t / r
    duy_dx = duy_dr * cos_t - duy_dth * sin_t / r
    duy_dy = duy_dr * sin_t + duy_dth * cos_t / r

    return np.array([[dux_dx, dux_dy], [duy_dx, duy_dy]])


class L2ErrorOperator(ScalarLocalOperator):
    """Integrand |u_h - u_ex|^2 for the L2 error norm."""

    def __init__(self, u_exact_fn):
        super().__init__()
        self.u_exact_fn = u_exact_fn

    def compute_scalar_integrand(self, R, dRdx, dRdy, physical_point, u_local):
        R = np.asarray(R)
        u_x, u_y = np.asarray(u_local)[0::2], np.asarray(u_local)[1::2]
        u_ex = self.u_exact_fn(physical_point)
        return (R @ u_x - u_ex[0]) ** 2 + (R @ u_y - u_ex[1]) ** 2


class H1SemiNormErrorOperator(ScalarLocalOperator):
    """Integrand for the H1 semi-norm error (gradient mismatch)."""

    def __init__(self, grad_exact_fn):
        super().__init__()
        self.grad_exact_fn = grad_exact_fn

    def compute_scalar_integrand(self, R, dRdx, dRdy, physical_point, u_local):
        dRdx, dRdy = np.asarray(dRdx), np.asarray(dRdy)
        u_x, u_y = np.asarray(u_local)[0::2], np.asarray(u_local)[1::2]
        grad_ex = self.grad_exact_fn(physical_point)
        return (
            (dRdx @ u_x - grad_ex[0, 0]) ** 2 + (dRdy @ u_x - grad_ex[0, 1]) ** 2
            + (dRdx @ u_y - grad_ex[1, 0]) ** 2 + (dRdy @ u_y - grad_ex[1, 1]) ** 2
        )


# ---------------------------------------------------------------------------
# Shared geometry: pymfiga's GeomdlGenerator is the single source of truth
# (examples/future/17_future_vs_pymfiga_matrix_free.ipynb cells 5-8, already
# validated bit-for-bit control-point-identical against a hand-built future
# patch). Convention: direction 0 (XI) = radial, direction 1 (ETA) = angular.
# ---------------------------------------------------------------------------

def build_pymfiga_patch(degree, nbel, Rin, Rex):
    geo_args = {"degree": degree, "nbel": nbel, "parameters": {"Rin": Rin, "Rex": Rex}}
    geometry = GeomdlGenerator(filename="nurbs_quarter_annulus", geo_args=geo_args).export_geometry()
    pp = SinglePatch(geometry, quadclass="wq")
    pp.generate()
    return pp


def build_future_from_pymfiga(pp):
    deg_u, deg_v = int(pp.degree[0]), int(pp.degree[1])
    kv_u, kv_v = pp.knotvector[0], pp.knotvector[1]
    nu_, nv_ = int(pp.nbctrlpts[0]), int(pp.nbctrlpts[1])
    ctrlpts, weights = pp.ctrlpts, pp.nurbs_weights  # (2, n), XI-fastest

    mgr = ControlPointManager(dim=2)
    for i in range(nu_ * nv_):
        mgr.add_point([float(ctrlpts[0, i]), float(ctrlpts[1, i])], float(weights[i]))

    su, sv = BSpline(deg_u, np.asarray(kv_u)), BSpline(deg_v, np.asarray(kv_v))
    patch0 = FPatch(BSplineSurface(su, sv), mgr, list(range(nu_ * nv_)), [nu_, nv_])
    n_cp = nu_ * nv_
    gdm = GlobalDOFManager([2] * n_cp)
    pdm = PatchDOFManager(2, list(range(n_cp)), gdm)
    pk = FPatch(patch0.tensor, patch0.cp_manager, list(patch0.global_indices), list(patch0.local_shape), pdm)
    return pk, pdm


def future_fixed_dofs(pk, pdm):
    """Symmetry BCs: ETA-min (theta=0, x-axis) blocks u_y, ETA-max (theta=90deg,
    y-axis) blocks u_x -- not a one-sided full clamp."""
    fixed = set()
    for cp in pk.boundary_control_points(1, 0):
        fixed.add(pdm.get_global_dof_indices(cp)[1])
    for cp in pk.boundary_control_points(1, 1):
        fixed.add(pdm.get_global_dof_indices(cp)[0])
    return np.array(sorted(fixed))


def find_cp_permutation(coords_a, coords_b, tol=1e-6):
    """perm[i] = index in coords_b matching coords_a[i]. Used to reconcile
    legacy's post-refine_patch() control-point order (an internal, untrusted
    Fortran algorithm) against the shared future/pymfiga order, by matching
    physical coordinates -- valid because NURBS degree elevation and knot
    insertion have a unique, algorithm-independent result for a given target
    (degree, knot vector), confirmed to ~1e-15 in practice."""
    n = coords_a.shape[0]
    if coords_b.shape[0] != n:
        raise ValueError(f"control point count mismatch: {n} vs {coords_b.shape[0]}")
    perm = np.full(n, -1, dtype=int)
    used = np.zeros(n, dtype=bool)
    max_dist = 0.0
    for i in range(n):
        d = np.linalg.norm(coords_b - coords_a[i], axis=1)
        d_masked = np.where(used, np.inf, d)
        j = int(np.argmin(d_masked))
        max_dist = max(max_dist, d_masked[j])
        perm[i] = j
        used[j] = True
    if max_dist > tol:
        raise ValueError(f"max matching distance {max_dist:.3e} exceeds tol {tol:.3e}")
    if len(set(perm.tolist())) != n:
        raise ValueError("coordinate matching is not a bijection")
    return perm


# ---------------------------------------------------------------------------
# Legacy pipeline.
#
# IMPORTANT: the legacy Fortran/f2py layer is not safe for more than one
# IGAparametrization / build_stiffness_matrix() per process -- confirmed
# empirically: a second legacy model built in the same process after a first
# one corrupts the heap (glibc abort), even via the officially-tested "3D
# solid"/"3D shell" paths, unrelated to anything specific to this script. The
# `worker` subcommand's one-case-per-subprocess design (see cmd_sweep) avoids
# this by construction -- never call build_legacy_model()/build_stiffness_matrix()
# more than once in the same process.
#
# `IgaModel("2D solid")` is itself an incomplete stub: __init__ never sets
# `_mcrd` for '2D solid', and add_patch() hard-codes `_TENSOR='THREED'` for a
# 'U1' element regardless of model_type. The fix below mirrors exactly what
# the working .inp-based reader path already does for `TYPE=U1,
# COORDINATES=2, TENSOR=PSTRESS` (see benchs/plateWithHole/plateWithHole.inp,
# which uses this same U1/2D/tensor combination via PSTRAIN) -- not a blind
# hack into unverified territory.
#
# Boundary conditions and the load vector are handled entirely in Python
# (see run_legacy below) rather than via add_boundary_condition()/
# add_distributed_load(): (a) legacy has no public API for a
# position-dependent VECTOR traction (Kirsch's exact traction has both a
# normal and a tangential component -- legacy's only native load mechanism is
# a scalar pressure along the outward normal) or for evaluating a solution at
# an arbitrary point (needed for the H1/L2 error), so the shared load vector
# is computed once via future's integrate_boundary_load() and injected
# directly, and legacy's solved vector is fed back into future's
# PatchIntegrator/error operators; (b) add_boundary_condition() must be
# called before refine_patch(), which conflicts with identifying boundary
# control points only AFTER refinement via the coordinate-matching
# permutation above -- so free/fixed DOFs are sliced in Python instead,
# exactly like the future/pymfiga pipelines already do.
# ---------------------------------------------------------------------------

def build_legacy_model(degree, nbel, Rin, Rex, E, nu):
    if nbel & (nbel - 1) != 0:
        raise ValueError("nbel must be a power of 2")
    if degree < BASE_DEGREE:
        raise ValueError(f"legacy pipeline requires degree >= {BASE_DEGREE}")
    n_levels = nbel.bit_length() - 1

    pp_base = build_pymfiga_patch(BASE_DEGREE, 1, Rin, Rex)
    nu_, nv_ = int(pp_base.nbctrlpts[0]), int(pp_base.nbctrlpts[1])
    n_cp0 = nu_ * nv_
    ctrlpts0, weights0 = pp_base.ctrlpts, pp_base.nurbs_weights

    control_points = np.zeros((n_cp0, 3))  # _COORDS is hard-coded 3D even for 2D patches
    control_points[:, 0] = ctrlpts0[0, :]
    control_points[:, 1] = ctrlpts0[1, :]
    weights = np.asarray(weights0, dtype=float).reshape(-1)
    connectivity = np.array([list(range(n_cp0))], dtype=int)  # single base element
    spans = np.array([[BASE_DEGREE, BASE_DEGREE]], dtype=int)
    degrees = np.array([BASE_DEGREE, BASE_DEGREE])
    knot_vectors = [np.asarray(pp_base.knotvector[0]), np.asarray(pp_base.knotvector[1])]

    material = ElasticMaterial(young_modulus=E, poisson_ratio=nu)
    patch = LPatch(element_type="U1", degrees=degrees, knot_vectors=knot_vectors,
                   control_points=control_points, weights=weights,
                   connectivity=connectivity, spans=spans, material=material)

    model = IgaModel("2D solid")
    model.iga_param._mcrd = 2
    model.add_patch(patch)
    model.iga_param._TENSOR[-1] = "PSTRESS"
    assert model.iga_param._TENSOR[-1] == "PSTRESS", "TENSOR field silently truncated"

    model.refine_patch(ipatch=0,
                        nb_degree_elevation=np.array([degree - BASE_DEGREE, degree - BASE_DEGREE]),
                        nb_subdivision=np.array([n_levels, n_levels]))
    return model


def run_legacy(pk, pdm, F_future, fixed_future, Rin, Rex, E, nu, degree, nbel):
    n_cp = pk.n_cp
    ndof = 2 * n_cp

    t0 = time.perf_counter()
    model = build_legacy_model(degree, nbel, Rin, Rex, E, nu)
    coords_legacy = model.cp_coordinates[:, :2]
    coords_future = pk.local_control_point_view()
    perm = find_cp_permutation(coords_legacy, coords_future)  # perm[legacy_i] = future_j
    inv_perm = np.empty(n_cp, dtype=int)
    inv_perm[perm] = np.arange(n_cp)
    dof_perm = np.empty(ndof, dtype=int)  # dof_perm[future_dof] = legacy_dof
    for j in range(n_cp):
        dof_perm[2 * j], dof_perm[2 * j + 1] = 2 * inv_perm[j], 2 * inv_perm[j] + 1

    stiff, _rhs_unused = model.build_stiffness_matrix()
    t_assembly = time.perf_counter() - t0

    F_legacy = F_future[dof_perm]
    fixed_legacy = dof_perm[fixed_future]
    free_legacy = np.setdiff1d(np.arange(ndof), fixed_legacy)

    t0 = time.perf_counter()
    K_ff = stiff.tocsr()[free_legacy][:, free_legacy].tocsc()
    u_legacy = np.zeros(ndof)
    u_legacy[free_legacy] = spla.spsolve(K_ff, F_legacy[free_legacy])
    t_solve = time.perf_counter() - t0

    u_future_order = np.empty(ndof)
    u_future_order[dof_perm] = u_legacy
    return dict(u=u_future_order, t_assembly=t_assembly, t_solve=t_solve, n_iter=1)


# ---------------------------------------------------------------------------
# `future` Gauss+direct and matrix-free pipelines.
# ---------------------------------------------------------------------------

def run_future_gauss_direct(pk, pdm, F, fixed, E, nu):
    law = PlaneStress(Material(E=E, nu=nu))
    ndof = 2 * pk.n_cp
    free = np.setdiff1d(np.arange(ndof), fixed)

    # basis_u/basis_v must stay alive in the same scope as the PatchIntegrator
    # using them: PatchIntegrator stores non-owning references (no pybind11
    # keep-alive), so once the Python objects that built them go out of scope
    # the references dangle (confirmed while building the source notebooks).
    basis_u = IGABasis1D.build(pk.tensor.components[0], pk.tensor.components[0].degree + 1)
    basis_v = IGABasis1D.build(pk.tensor.components[1], pk.tensor.components[1].degree + 1)
    integrator = PatchIntegrator(pk, basis_u, basis_v, law)

    t0 = time.perf_counter()
    K = integrator.integrate_stiffness()
    t_assembly = time.perf_counter() - t0

    K_ff = K.tocsr()[free][:, free].tocsc()
    t0 = time.perf_counter()
    u = np.zeros(ndof)
    u[free] = spla.spsolve(K_ff, F[free])
    t_solve = time.perf_counter() - t0
    return dict(u=u, t_assembly=t_assembly, t_solve=t_solve, n_iter=1)


def run_future_matrix_free(pk, pdm, F, fixed, E, nu, rtol, maxiter, restart):
    law = PlaneStress(Material(E=E, nu=nu))
    ndof = 2 * pk.n_cp

    t0 = time.perf_counter()
    wq = WQMatrixFreeStiffness(pk, law)
    precond = FastDiagonalizationPreconditioner(pk, pdm, n_dofs_per_cp=2, fixed_dofs=fixed, law=law)
    t_assembly = time.perf_counter() - t0

    it = [0]

    def cb(_xk):
        it[0] += 1

    t0 = time.perf_counter()
    u, info = mf_solve(wq.apply, ndof, F, fixed, method="gmres",
                        precondition_fn=precond.apply, rtol=rtol, maxiter=maxiter,
                        callback=cb, restart=restart)
    t_solve = time.perf_counter() - t0
    if info != 0:
        raise RuntimeError(f"future matrix-free GMRES did not converge (info={info})")
    return dict(u=u, t_assembly=t_assembly, t_solve=t_solve, n_iter=it[0])


def run_pymfiga_matrix_free(pp, Rin, Tx, E, nu, rtol, maxiter, restart):
    n_cp = int(pp.nbctrlpts[0]) * int(pp.nbctrlpts[1])
    ndof = 2 * n_cp

    # pymfiga's LinearElasticity implements the plane-STRAIN Lame formula
    # unconditionally (no plane-stress option) -- feed it the apparent
    # (E*, nu*) that makes its plane-strain formula reproduce the true
    # plane-stress (E, nu) constitutive matrix exactly (standard closed
    # form; verified against future's PlaneStress to ~1e-15 relative
    # Frobenius difference on the assembled stiffness operator):
    #   nu* = nu / (1 + nu),  E* = E * (1 + 2*nu) / (1 + nu)**2
    nu_star = nu / (1 + nu)
    E_star = E * (1 + 2 * nu) / (1 + nu) ** 2

    t0 = time.perf_counter()
    material = LinearElasticity({"elastic_modulus": E_star, "poisson_ratio": nu_star})
    boundary = BoundaryCondition(nbctrlpts=pp.nbctrlpts, dofs=(DOF.UX, DOF.UY))
    boundary.add_constraint(
        constraint_info=[
            {"direction": ParametricDirection.ETA, "face": BoundarySide.MIN, "dofs": (DOF.UY,)},
            {"direction": ParametricDirection.ETA, "face": BoundarySide.MAX, "dofs": (DOF.UX,)},
        ],
        constraint_type="dirichlet",
    )
    model = MechanicalModel(material, pp, boundary)

    F = model.assemble_surface_force({
        DOF.ALL: ({"direction": ParametricDirection.XI, "face": BoundarySide.MAX}, kirsch_surface_force(Rin, Tx))
    })  # radial-max (r=Rex): outer arc, matching future's direction=0, side=1

    _free, fixed = model.get_free_and_constraint_nodes()
    fixed = np.array(sorted(fixed))

    precond = model.preconditioner
    # Geometry+material-aware correction (matches future's FastDiagonalizationPreconditioner
    # default method="mean" when law= is given): pymfiga computes this mean
    # coefficient itself but never wires it in unless asked.
    precond.add_scalar_space_time_correctors(
        mass_corrector=np.array(model.scalar_mean_mass).tolist(),
        stiffness_corrector=np.array(model.scalar_mean_stiffness).tolist(),
    )
    precond.update_space_eigenvalues((0.0, 1.0))  # A = 0*M + 1*K: pure stiffness preconditioner
    t_assembly = time.perf_counter() - t0

    it = [0]

    def cb(_xk):
        it[0] += 1

    # future's own matrix-free solver drives pymfiga's operator directly:
    # MechanicalModel.solve_linearized_system() is not usable for a pure
    # static K.u=F solve as vendored in this repo (see git history for the
    # static/dynamic compatibility fix), so this mirrors
    # 17_future_vs_pymfiga_matrix_free.ipynb's own approach exactly.
    t0 = time.perf_counter()
    u_blocked, info = mf_solve(model.compute_mf_stiffness, ndof, F, fixed, method="gmres",
                               precondition_fn=precond.apply_spatial_preconditioner,
                               rtol=rtol, maxiter=maxiter, callback=cb, restart=restart)
    t_solve = time.perf_counter() - t0
    if info != 0:
        raise RuntimeError(f"pymfiga matrix-free GMRES did not converge (info={info})")

    u = np.empty(ndof)  # blocked [all ux, then all uy] -> interleaved [ux0,uy0,...]
    u[0::2] = u_blocked[:n_cp]
    u[1::2] = u_blocked[n_cp:]
    return dict(u=u, t_assembly=t_assembly, t_solve=t_solve, n_iter=it[0])


# ---------------------------------------------------------------------------
# Error norms (against the exact Kirsch solution), evaluated via future's
# PatchIntegrator regardless of which pipeline produced `u` -- every pipeline's
# solution is mapped into the shared future patch's (cp, dof) order first.
# ---------------------------------------------------------------------------

def compute_errors(pk, u, Rin, Tx, E, nu):
    basis_u = IGABasis1D.build(pk.tensor.components[0], pk.tensor.components[0].degree + 1)
    basis_v = IGABasis1D.build(pk.tensor.components[1], pk.tensor.components[1].degree + 1)
    law = PlaneStress(Material(E=E, nu=nu))
    integrator = PatchIntegrator(pk, basis_u, basis_v, law)
    u_exact_fn = lambda pt: kirsch_displacement(pt, Rin, Tx, E, nu)
    grad_exact_fn = lambda pt: kirsch_displacement_gradient(pt, Rin, Tx, E, nu)
    err_l2_sq = integrator.integrate_scalar_operator(L2ErrorOperator(u_exact_fn), u)
    err_h1_sq = integrator.integrate_scalar_operator(H1SemiNormErrorOperator(grad_exact_fn), u)
    return float(np.sqrt(err_l2_sq)), float(np.sqrt(err_l2_sq + err_h1_sq))


# ---------------------------------------------------------------------------
# One (degree, nbel, pipeline) case.
# ---------------------------------------------------------------------------

def run_one_case(pipeline, degree, nbel, Rin, Rex, Tx, E, nu, rtol, maxiter, restart):
    pp = build_pymfiga_patch(degree, nbel, Rin, Rex)
    pk, pdm = build_future_from_pymfiga(pp)
    ndof = 2 * pk.n_cp
    fixed = future_fixed_dofs(pk, pdm)

    basis_u = IGABasis1D.build(pk.tensor.components[0], pk.tensor.components[0].degree + 1)
    basis_v = IGABasis1D.build(pk.tensor.components[1], pk.tensor.components[1].degree + 1)
    law = PlaneStress(Material(E=E, nu=nu))
    F = PatchIntegrator(pk, basis_u, basis_v, law).integrate_boundary_load(0, 1, KirschTraction(Rin, Tx))

    if pipeline == "legacy":
        result = run_legacy(pk, pdm, F, fixed, Rin, Rex, E, nu, degree, nbel)
    elif pipeline == "future-gauss":
        result = run_future_gauss_direct(pk, pdm, F, fixed, E, nu)
    elif pipeline == "future-mf":
        result = run_future_matrix_free(pk, pdm, F, fixed, E, nu, rtol, maxiter, restart)
    elif pipeline == "pymfiga-mf":
        result = run_pymfiga_matrix_free(pp, Rin, Tx, E, nu, rtol, maxiter, restart)
    else:
        raise ValueError(f"unknown pipeline {pipeline!r}")

    l2_error, h1_error = compute_errors(pk, result["u"], Rin, Tx, E, nu)
    return dict(ndof=ndof, l2_error=l2_error, h1_error=h1_error,
                t_assembly=result["t_assembly"], t_solve=result["t_solve"],
                n_iter=result["n_iter"])


# ---------------------------------------------------------------------------
# `worker`: run exactly one case, write one JSON result.
# ---------------------------------------------------------------------------

def cmd_worker(args):
    logging.disable(logging.CRITICAL)  # both future and pymfiga are verbose by default
    row = dict(degree=args.degree, nbel=args.nbel, pipeline=args.pipeline,
               ndof=None, l2_error=None, h1_error=None, wall_time_s=None,
               assembly_time_s=None, solve_time_s=None, n_iter=None,
               peak_rss_mb=None, status="failed", error_message=None)
    t0 = time.perf_counter()
    try:
        r = run_one_case(args.pipeline, args.degree, args.nbel,
                         args.rin, args.rex, args.tx, args.E, args.nu,
                         args.rtol, args.maxiter, args.restart)
        row.update(ndof=r["ndof"], l2_error=r["l2_error"], h1_error=r["h1_error"],
                   assembly_time_s=r["t_assembly"], solve_time_s=r["t_solve"],
                   n_iter=r["n_iter"], status="ok")
    except Exception as e:
        row["error_message"] = f"{type(e).__name__}: {e}\n{traceback.format_exc()[-2000:]}"
    row["wall_time_s"] = time.perf_counter() - t0
    row["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

    with open(args.result_file, "w") as f:
        json.dump(row, f)
    print(json.dumps(row))
    sys.exit(0 if row["status"] == "ok" else 1)


# ---------------------------------------------------------------------------
# `sweep`: drive the (degree, nbel, pipeline) grid, one subprocess per case.
# ---------------------------------------------------------------------------

CSV_COLUMNS = ["degree", "nbel", "pipeline", "ndof", "l2_error", "h1_error",
               "wall_time_s", "assembly_time_s", "solve_time_s", "n_iter",
               "peak_rss_mb", "status", "error_message"]


def _read_done_cases(csv_path):
    done = set()
    if not csv_path.exists():
        return done
    import csv
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("status") == "ok":
                done.add((int(row["degree"]), int(row["nbel"]), row["pipeline"]))
    return done


def _append_csv_row(csv_path, row):
    import csv
    is_new = not csv_path.exists()
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        if is_new:
            writer.writeheader()
        writer.writerow({k: row.get(k) for k in CSV_COLUMNS})
        f.flush()
        os.fsync(f.fileno())


def cmd_sweep(args):
    if any(d < BASE_DEGREE for d in args.degrees):
        raise ValueError(f"degree must be >= {BASE_DEGREE} (a degree-1 NURBS circular "
                         f"arc cannot exist -- use the `abaqus` subcommand for degree 1)")
    for n in args.nbels:
        if n & (n - 1) != 0:
            raise ValueError(f"nbel must be a power of 2, got {n}")

    csv_path = Path(args.results_csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    done = _read_done_cases(csv_path) if args.resume else set()

    grid = [(d, n, p) for d in args.degrees for n in args.nbels for p in args.pipelines]
    for degree, nbel, pipeline in grid:
        if (degree, nbel, pipeline) in done:
            print(f"[skip, already ok] degree={degree} nbel={nbel} pipeline={pipeline}")
            continue
        print(f"[running] degree={degree} nbel={nbel} pipeline={pipeline} ...", flush=True)
        if args.dry_run:
            continue

        result_file = csv_path.parent / f".result_{degree}_{nbel}_{pipeline}.json"
        cmd = [sys.executable, str(Path(__file__).resolve()), "worker",
               "--pipeline", pipeline, "--degree", str(degree), "--nbel", str(nbel),
               "--n-cpus", str(args.n_cpus), "--rtol", str(args.rtol),
               "--maxiter", str(args.maxiter), "--restart", str(args.restart),
               "--rin", str(args.rin), "--rex", str(args.rex), "--tx", str(args.tx),
               "--E", str(args.E), "--nu", str(args.nu),
               "--result-file", str(result_file)]

        row = dict(degree=degree, nbel=nbel, pipeline=pipeline, status="failed")
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout_s)
            if result_file.exists():
                with open(result_file) as f:
                    row = json.load(f)
                result_file.unlink()
            else:
                row["error_message"] = (proc.stderr or "")[-2000:]
        except subprocess.TimeoutExpired:
            row["status"] = "timeout"
            row["error_message"] = f"exceeded --timeout-s={args.timeout_s}"
            if result_file.exists():
                result_file.unlink()

        _append_csv_row(csv_path, row)
        print(f"  -> status={row['status']}"
             + (f" h1_error={row.get('h1_error')}" if row.get("h1_error") is not None else "")
             + (f" wall_time_s={row.get('wall_time_s'):.2f}" if row.get("wall_time_s") else ""))

    print("SWEEP_DONE")


# ---------------------------------------------------------------------------
# `abaqus`: degree-1 deck, one per nbel (independent of the degree sweep).
#
# A genuine degree-1 NURBS circular arc cannot exist (BASE_DEGREE note above),
# so the degree-1 nodes are obtained by evaluating the shared (degree>=2)
# geometry on the uniform parametric grid (i/nbel, j/nbel) -- for a degree-1
# B-spline, control points coincide exactly with the curve at the knots
# (Greville abscissae), so this is the exact node placement a genuine
# degree-1 reparametrization at this knot spacing would have, not an
# approximation of one.
# ---------------------------------------------------------------------------

def _degree1_knots(nbel):
    if nbel == 1:
        return np.array([0.0, 0.0, 1.0, 1.0])
    interior = np.linspace(0.0, 1.0, nbel + 1)[1:-1]
    return np.concatenate([[0.0, 0.0], interior, [1.0, 1.0]])


def build_abaqus_patch(nbel, Rin, Rex):
    pp = build_pymfiga_patch(BASE_DEGREE, nbel, Rin, Rex)
    pk, pdm = build_future_from_pymfiga(pp)

    su, sv = pk.tensor.components
    params = np.array([[i / nbel, j / nbel] for j in range(nbel + 1) for i in range(nbel + 1)])
    spans_arr = np.array([[su.find_span(u), sv.find_span(v)] for u, v in params], dtype=np.int32)
    phys_pts = pk.evaluate_patch_nd_omp(spans_arr, params)  # (nbel+1)^2, XI-fastest

    n1 = nbel + 1
    mgr = ControlPointManager(dim=2)
    for i in range(n1 * n1):
        mgr.add_point([float(phys_pts[i, 0]), float(phys_pts[i, 1])], 1.0)
    su1 = BSpline(1, _degree1_knots(nbel))
    sv1 = BSpline(1, _degree1_knots(nbel))
    patch1 = FPatch(BSplineSurface(su1, sv1), mgr, list(range(n1 * n1)), [n1, n1])
    gdm1 = GlobalDOFManager([2] * (n1 * n1))
    pdm1 = PatchDOFManager(2, list(range(n1 * n1)), gdm1)
    pk1 = FPatch(patch1.tensor, patch1.cp_manager, list(patch1.global_indices), list(patch1.local_shape), pdm1)
    return pk1, pdm1


def write_abaqus_deck(path, pk1, pdm1, nbel, Rin, Tx, E, nu):
    n1 = nbel + 1
    coords = pk1.local_control_point_view()  # (n1*n1, 2), XI-fastest -> node id = idx+1
    fixed = future_fixed_dofs(pk1, pdm1)
    fixed_by_node = {}
    for dof in fixed:
        cp, comp = dof // 2, dof % 2
        fixed_by_node.setdefault(cp + 1, []).append(comp + 1)

    basis_u = IGABasis1D.build(pk1.tensor.components[0], pk1.tensor.components[0].degree + 1)
    basis_v = IGABasis1D.build(pk1.tensor.components[1], pk1.tensor.components[1].degree + 1)
    law = PlaneStress(Material(E=E, nu=nu))
    F = PatchIntegrator(pk1, basis_u, basis_v, law).integrate_boundary_load(0, 1, KirschTraction(Rin, Tx))
    loaded_nodes = [(cp + 1, F[2 * cp], F[2 * cp + 1]) for cp in range(n1 * n1) if abs(F[2 * cp]) + abs(F[2 * cp + 1]) > 0]

    # Standard structured Q4 connectivity, XI-fastest local numbering (same as
    # every geometry builder in this script). Abaqus's CPS4 wants
    # counterclockwise node order (a clockwise quad gives a negative
    # Jacobian); this quarter-annulus parametrization's (i, j) traversal
    # happens to be physically clockwise, so the orientation is checked via
    # the shoelace formula on the first element and every element is flipped
    # accordingly, rather than assuming either order.
    def node_id(i, j):
        return j * n1 + i + 1

    def quad_signed_area(n_a, n_b, n_c, n_d):
        pts = coords[[n_a - 1, n_b - 1, n_c - 1, n_d - 1]]
        x, y = pts[:, 0], pts[:, 1]
        return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)

    raw = [(node_id(i, j), node_id(i + 1, j), node_id(i + 1, j + 1), node_id(i, j + 1))
           for j in range(nbel) for i in range(nbel)]
    if quad_signed_area(*raw[0]) < 0:
        raw = [(a, d, c, b) for (a, b, c, d) in raw]
    elements = raw

    with open(path, "w") as f:
        f.write(f"** Kirsch benchmark, degree-1 (classical FE) mesh, nbel={nbel} per side\n")
        f.write(f"** ABAQUS_LAUNCH_CMD (fill in at the top of kirsch_sweep.py): {ABAQUS_LAUNCH_CMD!r}\n")
        f.write("*HEADING\nKirsch plate with hole, degree-1 mesh\n")
        f.write("*NODE\n")
        for idx in range(n1 * n1):
            f.write(f"{idx + 1}, {coords[idx, 0]:.10g}, {coords[idx, 1]:.10g}, 0.0\n")
        f.write("*ELEMENT, TYPE=CPS4, ELSET=ALLELS\n")
        for eid, (n_a, n_b, n_c, n_d) in enumerate(elements, start=1):
            f.write(f"{eid}, {n_a}, {n_b}, {n_c}, {n_d}\n")
        f.write("*MATERIAL, NAME=STEEL\n*ELASTIC\n")
        f.write(f"{E:.10g}, {nu:.10g}\n")
        f.write("*SOLID SECTION, ELSET=ALLELS, MATERIAL=STEEL\n1.0,\n")
        f.write("*BOUNDARY\n")
        for node, comps in sorted(fixed_by_node.items()):
            for c in comps:
                f.write(f"{node}, {c}, {c}, 0.0\n")
        # *CLOAD is history data (a load is always step-specific), not model
        # data -- it must be inside *STEP, unlike *BOUNDARY above which is
        # valid either as model data (the base/initial state, used here) or
        # inside a step. Abaqus itself confirmed this the only way it could
        # (a real run): "*CLOAD ... misplaced ... suboption for ... step".
        f.write("*STEP, PERTURBATION\n*STATIC\n")
        f.write("*CLOAD\n")
        for node, fx, fy in loaded_nodes:
            if abs(fx) > 0:
                f.write(f"{node}, 1, {fx:.10g}\n")
            if abs(fy) > 0:
                f.write(f"{node}, 2, {fy:.10g}\n")
        f.write("*OUTPUT, FIELD\n*NODE OUTPUT\nU\n*END STEP\n")


def parse_abaqus_dat_memory_mb(dat_text):
    """Abaqus's own "MEMORY TO MINIMIZE I/O" figure (MB) from the .dat
    file's MEMORY ESTIMATE table, summed across all PROCESS rows (Abaqus
    prints one row per process for a multi-process/MPI run; direct-solver
    thread-parallel runs, as here, print a single row regardless of
    `cpus=`). Returns None if the table isn't found (e.g. the job failed
    before writing it).

    This is used instead of an OS-level RUSAGE_CHILDREN measurement: that
    approach was tried first and found unreliable here -- ru_maxrss for
    RUSAGE_CHILDREN is a running MAXIMUM across every child the parent
    process has ever reaped, not a per-call gauge, so subtracting
    before/after snapshots around one Abaqus job silently reads back 0 (and
    was then coerced to None by an `or None`) whenever an earlier job or the
    odb-extraction step in the same `cmd_abaqus --nbels ...` loop had
    already pushed the watermark at or above this job's own peak -- exactly
    what happened for nbel=8 (see this file's own git history), not
    something a smaller/larger --nbels list would reliably dodge. Abaqus's
    own figure has no such cross-call contamination."""
    idx = dat_text.find("M E M O R Y   E S T I M A T E")
    if idx == -1:
        return None
    rows = re.findall(r"^\s*\d+\s+[\d.]+E[+-]\d+\s+\d+\s+(\d+)\s*$", dat_text[idx:], re.MULTILINE)
    if not rows:
        return None
    return float(sum(int(v) for v in rows))


def run_abaqus_job(out_dir, nbel, n_cpus):
    """Actually launches Abaqus (ABAQUS_LAUNCH_CMD, shell-formatted with
    nbel=nbel, n_cpus=n_cpus), blocking until it finishes. Returns
    (odb_path, job_name, wall_time_s, peak_rss_mb). Raises on a nonzero exit
    code OR a missing .odb afterwards -- Abaqus's own exit code isn't
    always a reliable success signal by itself."""
    if not ABAQUS_LAUNCH_CMD:
        raise RuntimeError("ABAQUS_LAUNCH_CMD is empty -- fill it in at the top of this file")
    job_name = f"kirsch_nbel{nbel}"
    cmd = ABAQUS_LAUNCH_CMD.format(nbel=nbel, n_cpus=n_cpus)
    print(f"  launching: {cmd}  (cwd={out_dir})")
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, shell=True, cwd=str(out_dir), capture_output=True, text=True, check=False)
    wall_time_s = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"Abaqus exited with code {proc.returncode}:\n"
                           f"{(proc.stdout or '')[-1000:]}\n{(proc.stderr or '')[-1000:]}")
    odb_path = out_dir / f"{job_name}.odb"
    if not odb_path.exists():
        raise RuntimeError(f"Abaqus reported success (exit 0) but {odb_path} was never "
                           f"created -- check {out_dir / (job_name + '.log')} / '.sta' / '.msg'")
    dat_path = out_dir / f"{job_name}.dat"
    peak_rss_mb = parse_abaqus_dat_memory_mb(dat_path.read_text()) if dat_path.exists() else None
    return odb_path, job_name, wall_time_s, peak_rss_mb


def _abaqus_launcher_prefix():
    """Everything ABAQUS_LAUNCH_CMD runs before the `job=...` arguments --
    e.g. `singularity exec --nv some/abaqus2019....sif abaqus` -- so the odb
    extraction (`abaqus python ...`) goes through the exact same
    singularity/module/wrapper as the solve itself, instead of assuming a
    bare `abaqus` is on PATH (it usually isn't, e.g. inside a container
    image). `\\babaqus\\b` matches only the standalone command word, not an
    image filename like `abaqus2019_....sif` that happens to start with it."""
    m = re.search(r"^(.*?)\babaqus\b", ABAQUS_LAUNCH_CMD)
    if not m:
        raise RuntimeError(f"could not find a standalone 'abaqus' word in "
                           f"ABAQUS_LAUNCH_CMD to reuse for odb extraction: {ABAQUS_LAUNCH_CMD!r}")
    return m.group(0)


def run_abaqus_extraction(out_dir, odb_path, job_name):
    """Launches extract_abaqus_results.py through Abaqus's own bundled
    Python (`abaqus python ...`, via the same launcher prefix as
    ABAQUS_LAUNCH_CMD), which is the only place `odbAccess` exists. Returns
    the CSV path it wrote."""
    csv_path = out_dir / f"{job_name}_results.csv"
    script_path = Path(__file__).resolve().parent / "extract_abaqus_results.py"
    # subprocess.run below sets cwd=out_dir (matching run_abaqus_job, in case
    # the launcher prefix itself assumes it -- e.g. singularity bind mounts),
    # so odb_path/csv_path (constructed relative to THIS process's cwd, not
    # out_dir) must be made absolute first or they'd resolve to out_dir/out_dir/...
    odb_abs = Path(odb_path).resolve()
    csv_abs = Path(csv_path).resolve()
    cmd = f"{_abaqus_launcher_prefix()} python {script_path} {odb_abs} {csv_abs}"
    print(f"  extracting: {cmd}")
    proc = subprocess.run(cmd, shell=True, cwd=str(out_dir), capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"odb extraction exited with code {proc.returncode}:\n"
                           f"{(proc.stdout or '')[-1000:]}\n{(proc.stderr or '')[-1000:]}")
    if not csv_path.exists():
        raise RuntimeError(f"extraction reported success (exit 0) but {csv_path} was never created")
    return csv_path


def cmd_abaqus(args):
    logging.disable(logging.CRITICAL)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for nbel in args.nbels:
        if nbel & (nbel - 1) != 0:
            raise ValueError(f"nbel must be a power of 2, got {nbel}")
        pk1, pdm1 = build_abaqus_patch(nbel, args.rin, args.rex)
        path = out_dir / f"kirsch_nbel{nbel}.inp"
        write_abaqus_deck(path, pk1, pdm1, nbel, args.rin, args.tx, args.E, args.nu)
        print(f"wrote {path}  ({(nbel + 1) ** 2} nodes, {nbel * nbel} elements)")

        if not args.run:
            continue
        # A failure here (Abaqus itself, or the odb extraction) is reported
        # and skipped rather than aborting the whole --nbels list, matching
        # `sweep`'s fault-tolerant per-case handling.
        try:
            odb_path, job_name, wall_time_s, peak_rss_mb = run_abaqus_job(out_dir, nbel, args.n_cpus)
            csv_path = run_abaqus_extraction(out_dir, odb_path, job_name)
            append_abaqus_result(csv_path, nbel, args.rin, args.rex, args.tx, args.E, args.nu,
                                 args.results_csv, wall_time_s=wall_time_s, peak_rss_mb=peak_rss_mb)
        except Exception as e:
            print(f"  FAILED (nbel={nbel}): {type(e).__name__}: {e}")

    if not args.run and not ABAQUS_LAUNCH_CMD:
        print("\nABAQUS_LAUNCH_CMD is empty -- fill it in at the top of this file "
             "before launching Abaqus yourself (this script never invokes it without --run).")


# ---------------------------------------------------------------------------
# `abaqus-results`: pull a finished Abaqus run's nodal displacements back into
# the same results CSV as the 4 internal pipelines, under pipeline="abaqus"
# (degree=1), so it shows up in `plot`'s 3 figures alongside them.
#
# Two-step process, since this script cannot run Abaqus or import odbAccess
# itself (odbAccess only exists inside Abaqus's own bundled Python):
#   1. On the machine/account with Abaqus, after the job finishes, run
#      `abaqus python extract_abaqus_results.py <job>.odb <job>_results.csv`
#      (extract_abaqus_results.py, next to this file) -- it writes a plain
#      CSV (node_id,x,y,ux,uy), no odbAccess needed on this side.
#   2. Here: `kirsch_sweep.py abaqus-results --csv-in <job>_results.csv
#      --nbel <N> --results-csv results/sweep_results.csv`.
#
# The error norm is computed by feeding Abaqus's own nodal displacements into
# the SAME build_abaqus_patch()/compute_errors() used to build the .inp deck
# in the first place -- not a separate, untested code path -- and node
# coordinates are cross-checked against that patch's own geometry (not just
# trusted) to catch a wrong --nbel or a re-meshed odb early.
# ---------------------------------------------------------------------------

def read_abaqus_results_csv(path):
    import csv
    data = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[int(row["node_id"])] = (float(row["x"]), float(row["y"]),
                                         float(row["ux"]), float(row["uy"]))
    return data


def append_abaqus_result(csv_in, nbel, rin, rex, tx, E, nu, results_csv,
                         wall_time_s=None, peak_rss_mb=None):
    """Shared by cmd_abaqus_results (manual, CLI-driven) and cmd_abaqus's
    --run (automatic): read an extract_abaqus_results.py CSV, compute its H1
    error against the exact Kirsch solution, append one row to results_csv."""
    if nbel & (nbel - 1) != 0:
        raise ValueError(f"nbel must be a power of 2, got {nbel}")

    pk1, _pdm1 = build_abaqus_patch(nbel, rin, rex)
    n1 = nbel + 1
    n_cp = n1 * n1
    ndof = 2 * n_cp
    coords = pk1.local_control_point_view()

    data = read_abaqus_results_csv(csv_in)
    if len(data) != n_cp:
        raise ValueError(f"expected {n_cp} nodes for nbel={nbel} (from the deck "
                         f"this script would generate), found {len(data)} in {csv_in}")

    u = np.zeros(ndof)
    for node_id, (x, y, ux, uy) in data.items():
        cp = node_id - 1  # write_abaqus_deck's node_id = pk1's flat CP index + 1
        if not (0 <= cp < n_cp):
            raise ValueError(f"node_id {node_id} out of range for nbel={nbel}")
        if abs(coords[cp, 0] - x) > 1e-6 or abs(coords[cp, 1] - y) > 1e-6:
            raise ValueError(
                f"node {node_id} coordinate mismatch: this script's own mesh has "
                f"{tuple(coords[cp])}, the Abaqus export reports ({x}, {y}) -- likely "
                f"wrong --nbel, or the .odb wasn't run from the exact .inp this script "
                f"generates for nbel={nbel} (--rin/--rex must match too)")
        u[2 * cp], u[2 * cp + 1] = ux, uy

    l2_error, h1_error = compute_errors(pk1, u, rin, tx, E, nu)
    print(f"  nbel={nbel}  ndof={ndof}  l2_error={l2_error:.6e}  h1_error={h1_error:.6e}")

    row = dict(degree=1, nbel=nbel, pipeline="abaqus", ndof=ndof,
               l2_error=l2_error, h1_error=h1_error,
               wall_time_s=wall_time_s, assembly_time_s=None, solve_time_s=None,
               n_iter=None, peak_rss_mb=peak_rss_mb, status="ok", error_message=None)

    csv_path = Path(results_csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    _append_csv_row(csv_path, row)
    print(f"  appended to {csv_path}")


def cmd_abaqus_results(args):
    logging.disable(logging.CRITICAL)
    append_abaqus_result(args.csv_in, args.nbel, args.rin, args.rex, args.tx, args.E, args.nu,
                         args.results_csv, wall_time_s=args.wall_time_s, peak_rss_mb=args.peak_rss_mb)


# ---------------------------------------------------------------------------
# `plot`: 3 figures from the results CSV.
# ---------------------------------------------------------------------------

def cmd_plot(args):
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.results_csv)
    df.loc[df["status"] != "ok", ["h1_error", "wall_time_s", "peak_rss_mb"]] = np.nan

    degrees = sorted(df["degree"].unique())
    # "abaqus" (degree=1, from `abaqus-results`) is not one of the 4 internal
    # sweep pipelines in PIPELINES, but it shares the same CSV/plot -- listed
    # here explicitly rather than added to PIPELINES itself, which also
    # drives `sweep`'s/`worker`'s own --pipeline choices.
    pipelines = [p for p in PIPELINES + ["abaqus"] if p in set(df["pipeline"])]
    # Matplotlib's default "CN" cycle (tab10, indexed discretely) -- sampling
    # tab10 continuously via np.linspace(0, 1, N) skips most of the palette
    # for small N (e.g. red, C3, never appears with only 2-3 degrees).
    colors = [f"C{i % 10}" for i in range(len(degrees))]
    markers = {"legacy": "^", "future-gauss": "o", "future-mf": "s", "pymfiga-mf": "D", "abaqus": "*"}
    linestyles = {"legacy": ":", "future-gauss": "-", "future-mf": "--", "pymfiga-mf": "-.", "abaqus": (0, (1, 1))}

    specs = [
        ("h1_error", "H1 error vs exact Kirsch solution", "H1 error", True, "kirsch_h1_error"),
        ("wall_time_s", "Wall-clock time (assembly + solve)", "time (s)", True, "kirsch_wall_time"),
        ("peak_rss_mb", "Peak memory (RSS)", "peak RSS (MB)", True, "kirsch_peak_memory"),
    ]
    for col, title, ylabel, logy, basename in specs:
        fig, ax = plt.subplots(figsize=(8, 6))
        for i, degree in enumerate(degrees):
            for pipeline in pipelines:
                rows = df[(df["degree"] == degree) & (df["pipeline"] == pipeline)].sort_values("nbel")
                if rows.empty:
                    continue
                ax.plot(rows["nbel"], rows[col], marker=markers.get(pipeline, "x"),
                       linestyle=linestyles.get(pipeline, "-"), color=colors[i],
                       label=f"{pipeline}, degree {degree}")
        ax.set_xscale("log")  # base 10
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel("elements per side")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        for fmt in args.formats:
            path = out_dir / f"{basename}.{fmt}"
            fig.savefig(path, dpi=150)  # dpi only affects raster formats (png)
            print(f"wrote {path}")
        plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _add_physics_args(p):
    p.add_argument("--rin", type=float, default=Rin_DEFAULT)
    p.add_argument("--rex", type=float, default=Rex_DEFAULT)
    p.add_argument("--tx", type=float, default=Tx_DEFAULT)
    p.add_argument("--E", type=float, default=E_DEFAULT)
    p.add_argument("--nu", type=float, default=NU_DEFAULT)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    p_worker = sub.add_parser("worker", help="run exactly one (degree, nbel, pipeline) case")
    p_worker.add_argument("--pipeline", choices=PIPELINES, required=True)
    p_worker.add_argument("--degree", type=int, required=True)
    p_worker.add_argument("--nbel", type=int, required=True)
    p_worker.add_argument("--n-cpus", type=int, default=1)
    p_worker.add_argument("--rtol", type=float, default=1e-10)
    p_worker.add_argument("--maxiter", type=int, default=5000)
    p_worker.add_argument("--restart", type=int, default=200)
    p_worker.add_argument("--result-file", required=True)
    _add_physics_args(p_worker)
    p_worker.set_defaults(func=cmd_worker)

    p_sweep = sub.add_parser("sweep", help="run the full (degree, nbel, pipeline) grid")
    p_sweep.add_argument("--degrees", type=int, nargs="+", required=True)
    p_sweep.add_argument("--nbels", type=int, nargs="+", required=True)
    p_sweep.add_argument("--pipelines", nargs="+", choices=PIPELINES, default=PIPELINES)
    p_sweep.add_argument("--n-cpus", type=int, default=1)
    p_sweep.add_argument("--rtol", type=float, default=1e-10)
    p_sweep.add_argument("--maxiter", type=int, default=5000)
    p_sweep.add_argument("--restart", type=int, default=200)
    p_sweep.add_argument("--results-csv", default="results/sweep_results.csv")
    p_sweep.add_argument("--timeout-s", type=float, default=3600.0)
    p_sweep.add_argument("--resume", dest="resume", action="store_true", default=True)
    p_sweep.add_argument("--no-resume", dest="resume", action="store_false")
    p_sweep.add_argument("--dry-run", action="store_true")
    _add_physics_args(p_sweep)
    p_sweep.set_defaults(func=cmd_sweep)

    p_abaqus = sub.add_parser("abaqus", help="write degree-1 Abaqus decks, one per nbel")
    p_abaqus.add_argument("--nbels", type=int, nargs="+", required=True)
    p_abaqus.add_argument("--out-dir", default="results/abaqus")
    p_abaqus.add_argument("--run", action="store_true",
                          help="also launch Abaqus (ABAQUS_LAUNCH_CMD) and the odb "
                               "extraction for each deck, appending to --results-csv")
    p_abaqus.add_argument("--results-csv", default="results/sweep_results.csv",
                          help="only used with --run")
    p_abaqus.add_argument("--n-cpus", type=int, default=1,
                          help="only used with --run: forwarded to ABAQUS_LAUNCH_CMD as "
                               "n_cpus (e.g. Abaqus's own cpus=), and pins this script's "
                               "own OMP_NUM_THREADS/etc. the same way worker/sweep do")
    _add_physics_args(p_abaqus)
    p_abaqus.set_defaults(func=cmd_abaqus)

    p_abaqus_results = sub.add_parser(
        "abaqus-results",
        help="parse a finished Abaqus run's exported displacements and append its "
             "H1 error to the results CSV (see extract_abaqus_results.py)")
    p_abaqus_results.add_argument("--csv-in", required=True,
                                  help="CSV written by extract_abaqus_results.py (node_id,x,y,ux,uy)")
    p_abaqus_results.add_argument("--nbel", type=int, required=True)
    p_abaqus_results.add_argument("--results-csv", default="results/sweep_results.csv")
    p_abaqus_results.add_argument("--wall-time-s", type=float, default=None,
                                  help="optional, e.g. from Abaqus's own job log -- omit for none")
    p_abaqus_results.add_argument("--peak-rss-mb", type=float, default=None,
                                  help="optional, e.g. from Abaqus's own job log -- omit for none")
    _add_physics_args(p_abaqus_results)
    p_abaqus_results.set_defaults(func=cmd_abaqus_results)

    p_plot = sub.add_parser("plot", help="produce the 3 summary plots from a results CSV")
    p_plot.add_argument("--results-csv", default="results/sweep_results.csv")
    p_plot.add_argument("--out-dir", default="results/plots")
    p_plot.add_argument("--formats", nargs="+", default=["pdf"],
                        help="file extensions to save each figure as, e.g. --formats pdf png "
                             "(default: pdf, vector -- any matplotlib-supported extension works)")
    p_plot.set_defaults(func=cmd_plot)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
