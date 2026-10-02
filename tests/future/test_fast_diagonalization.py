"""
FastDiagonalizationPreconditioner: ported from pymfiga's
iga/fastdiagonalization/single/cspace.py (SpaceFD) -- validated as (a) a
pure linear-algebra identity (matches a direct dense solve of the same
Kronecker-sum approximation it implements matrix-free) and (b) an actual
iteration-count reduction when used with matrix_free_solver.solve().
"""

import numpy as np
import pytest

from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager,
    Patch, PRefiner, SubdivisionRefiner,
    GlobalDOFManager, PatchDOFManager, Material, PlaneStress, PlaneStrain,
    WeightedQuadrature1D, WQMatrixFreeStiffness,
)
from yeti_iga.future.fast_diagonalization import (
    FastDiagonalizationPreconditioner, compute_stiffness_correctors,
    compute_tensor_decomposition_correctors, compute_operator_diagonal,
)
from yeti_iga.future.matrix_free_solver import solve as mf_solve


def test_matches_dense_kronecker_sum_solve():
    """
    Pure linear-algebra check: apply() on a scalar (1 dof/cp) problem must
    match np.linalg.solve() against the SAME Kronecker-sum matrix
    (kron(Kv_free, Mu_free) + kron(Mv_free, Ku_free)) explicitly assembled.
    """
    degree, nbel = 2, 3
    kv = np.array([0.] * (degree + 1) + list(np.linspace(0, 1, nbel + 1)[1:-1]) + [1.] * (degree + 1))
    su, sv = BSpline(degree, kv), BSpline(degree, kv)
    nu = nv = len(kv) - degree - 1

    mgr = ControlPointManager(dim=2)
    mapping = []
    for jv in range(nv):
        for iu in range(nu):
            mgr.add_point([float(iu), float(jv)])
            mapping.append(jv * nu + iu)
    gdm = GlobalDOFManager([1] * mgr.n_points)
    pdm = PatchDOFManager(1, mapping, gdm)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [nu, nv], pdm)

    fixed_cps = patch.boundary_control_points(0, 0)
    fixed_dofs = [pdm.get_global_dof_indices(cp)[0] for cp in fixed_cps]

    precond = FastDiagonalizationPreconditioner(patch, pdm, n_dofs_per_cp=1, fixed_dofs=fixed_dofs)

    wq_u = WeightedQuadrature1D.build(su, "2")
    wq_v = WeightedQuadrature1D.build(sv, "2")
    Mu, Ku = (wq_u.W00 @ wq_u.B0).toarray(), (wq_u.W10 @ wq_u.B1).toarray()
    Mv, Kv = (wq_v.W00 @ wq_v.B0).toarray(), (wq_v.W10 @ wq_v.B1).toarray()

    free_u, free_v = np.arange(nu)[1:], np.arange(nv)
    Ku_f, Mu_f = Ku[np.ix_(free_u, free_u)], Mu[np.ix_(free_u, free_u)]
    Kv_f, Mv_f = Kv[np.ix_(free_v, free_v)], Mv[np.ix_(free_v, free_v)]
    K_dense = np.kron(Kv_f, Mu_f) + np.kron(Mv_f, Ku_f)

    rng = np.random.default_rng(0)
    r_free = rng.standard_normal(len(free_u) * len(free_v))
    expected_free = np.linalg.solve(K_dense, r_free)

    grid = np.zeros((nv, nu))
    grid[np.ix_(free_v, free_u)] = r_free.reshape(len(free_v), len(free_u))
    out = precond.apply(grid.ravel())
    out_free = out.reshape(nv, nu)[np.ix_(free_v, free_u)].ravel()

    assert np.linalg.norm(out_free - expected_free) / np.linalg.norm(expected_free) < 1e-9


def _build_quarter_ring(degree, n_levels):
    Rin, Rex = 1.0, 2.0
    w_c = 1.0 / np.sqrt(2.0)
    mgr = ControlPointManager(dim=2)
    for (x, y, w) in [
        (Rin, 0., 1.0), (Rin, Rin, w_c), (0., Rin, 1.0),
        (Rex, 0., 1.0), (Rex, Rex, w_c), (0., Rex, 1.0),
    ]:
        mgr.add_point([x, y], w)
    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    patch = Patch(BSplineSurface(su, sv), mgr, list(range(6)), [3, 2])
    if degree > 2:
        PRefiner(direction=0, n_elevations=degree - 2).refine(patch)
    PRefiner(direction=1, n_elevations=degree - 1).refine(patch)
    SubdivisionRefiner(direction=0, n_levels=n_levels).refine(patch)
    SubdivisionRefiner(direction=1, n_levels=n_levels).refine(patch)

    n_cp = patch.n_cp
    gdm = GlobalDOFManager([2] * n_cp)
    pdm = PatchDOFManager(2, list(range(n_cp)), gdm)
    pk = Patch(patch.tensor, patch.cp_manager, list(patch.global_indices), list(patch.local_shape), pdm)
    law = PlaneStress(Material(E=210000., nu=0.3))
    return pk, pdm, law


def test_reduces_iteration_count_and_keeps_it_flat_under_refinement():
    """
    The actual point of a preconditioner: fewer iterations, and -- unlike
    the unpreconditioned count, which grows sharply with mesh refinement --
    roughly CONSTANT iteration counts across refinement levels.
    """
    degree = 3
    counts_unprecond, counts_precond = [], []
    for n_levels in [2, 4]:
        patch, pdm, law = _build_quarter_ring(degree, n_levels)
        wq_stiff = WQMatrixFreeStiffness(patch, law)
        ndof = 2 * patch.n_cp

        fixed = set()
        for cp in patch.boundary_control_points(0, 0):
            fixed.add(pdm.get_global_dof_indices(cp)[1])
        for cp in patch.boundary_control_points(0, 1):
            fixed.add(pdm.get_global_dof_indices(cp)[0])
        fixed = np.array(sorted(fixed))

        rng = np.random.default_rng(0)
        F = rng.standard_normal(ndof)

        for counts, precond_fn in [(counts_unprecond, None), (counts_precond, "fd")]:
            it = [0]

            def cb(_xk, it=it):
                it[0] += 1

            if precond_fn == "fd":
                precond = FastDiagonalizationPreconditioner(patch, pdm, n_dofs_per_cp=2, fixed_dofs=fixed)
                pf = precond.apply
            else:
                pf = None

            u, info = mf_solve(wq_stiff.apply, ndof, F, fixed, method="gmres",
                               precondition_fn=pf, rtol=1e-10, maxiter=10000, callback=cb)
            assert info == 0
            counts.append(it[0])

    # Preconditioned: always fewer iterations than unpreconditioned, at both levels.
    assert all(cp < cu for cp, cu in zip(counts_precond, counts_unprecond))
    # Unpreconditioned count grows sharply under refinement (at least 5x here);
    # preconditioned count stays close to flat (well under 2x).
    assert counts_unprecond[1] > 5 * counts_unprecond[0]
    assert counts_precond[1] < 2 * counts_precond[0]


def test_stiffness_correctors_on_unit_square_match_closed_form():
    """
    On an affine, UNIT-SQUARE, unit-weight patch, the Jacobian is the
    identity everywhere (detJ=1, invJ=I), so compute_stiffness_correctors()
    must reduce to the bare Lame constants: c[m, l] = lambda + 2*mu if l==m
    else mu -- exactly the diagonal elasticity coefficients Section 2.3.2 of
    Cornejo Fuentes' thesis calls "conductivity-like". A closed-form check,
    independent of the WQMatrixFreeStiffness machinery being preconditioned.
    """
    degree, nbel = 2, 3
    su = BSpline(degree, np.array([0.] * (degree + 1) + [1.] * (degree + 1)))
    sv = BSpline(degree, np.array([0.] * (degree + 1) + [1.] * (degree + 1)))
    nu = nv = degree + 1
    mgr = ControlPointManager(dim=2)
    mapping = []
    for jv in range(nv):
        for iu in range(nu):
            mgr.add_point([iu / degree, jv / degree])
            mapping.append(jv * nu + iu)
    gdm = GlobalDOFManager([2] * mgr.n_points)
    pdm = PatchDOFManager(2, mapping, gdm)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [nu, nv], pdm)

    E, nu_poisson = 210000.0, 0.3
    law = PlaneStrain(Material(E=E, nu=nu_poisson))
    mu = E / (2.0 * (1.0 + nu_poisson))
    lam = nu_poisson * E / ((1.0 + nu_poisson) * (1.0 - 2.0 * nu_poisson))

    c = compute_stiffness_correctors(patch, law)
    expected = np.array([[lam + 2 * mu, mu], [mu, lam + 2 * mu]])
    assert np.allclose(c, expected, rtol=1e-10)


def test_law_aware_correction_reduces_iteration_count_further():
    """
    On the curved, rational quarter-ring, passing `law=` to
    FastDiagonalizationPreconditioner (the thesis's geometry+material-aware
    correction, Eq. (2.15)) must need fewer GMRES iterations than the
    uncorrected ("classic FD") default -- the entire point of the
    correction. Regression guard for compute_stiffness_correctors() and its
    wiring into the Kronecker-sum eigenvalues.
    """
    degree, n_levels = 4, 3
    patch, pdm, law = _build_quarter_ring(degree, n_levels)
    wq_stiff = WQMatrixFreeStiffness(patch, law)
    ndof = 2 * patch.n_cp

    fixed = set()
    for cp in patch.boundary_control_points(0, 0):
        fixed.add(pdm.get_global_dof_indices(cp)[1])
    for cp in patch.boundary_control_points(0, 1):
        fixed.add(pdm.get_global_dof_indices(cp)[0])
    fixed = np.array(sorted(fixed))

    rng = np.random.default_rng(0)
    F = rng.standard_normal(ndof)

    counts = {}
    for label, law_arg in [("classic", None), ("corrected", law)]:
        it = [0]

        def cb(_xk, it=it):
            it[0] += 1

        precond = FastDiagonalizationPreconditioner(
            patch, pdm, n_dofs_per_cp=2, fixed_dofs=fixed, law=law_arg)
        u, info = mf_solve(wq_stiff.apply, ndof, F, fixed, method="gmres",
                           precondition_fn=precond.apply, rtol=1e-10, maxiter=10000,
                           callback=cb, restart=200)
        assert info == 0
        counts[label] = it[0]

    assert counts["corrected"] < counts["classic"]


def test_tensor_decomposition_method_reduces_iterations_further_than_mean():
    """
    Montardini's per-quadrature-point tensor decomposition ("tensor" method,
    ported from iga_wq_mf/fsrc/tensor_algebra.f90's tensor_decomposition_2d
    -- the actual Fortran code behind Cornejo Fuentes' thesis) should need
    no more GMRES iterations than the coarser single-mean-coefficient
    "mean" method (compute_stiffness_correctors) on the same curved,
    rational quarter-ring -- it is a strictly finer approximation of the
    same underlying coefficient field. A small case is used deliberately:
    compute_tensor_decomposition_correctors is O(nq_u * nq_v) in pure
    Python and gets slow at the mesh sizes the notebooks use.
    """
    degree, n_levels = 3, 1
    patch, pdm, law = _build_quarter_ring(degree, n_levels)
    wq_stiff = WQMatrixFreeStiffness(patch, law)
    ndof = 2 * patch.n_cp

    fixed = set()
    for cp in patch.boundary_control_points(0, 0):
        fixed.add(pdm.get_global_dof_indices(cp)[1])
    for cp in patch.boundary_control_points(0, 1):
        fixed.add(pdm.get_global_dof_indices(cp)[0])
    fixed = np.array(sorted(fixed))

    rng = np.random.default_rng(0)
    F = rng.standard_normal(ndof)

    counts = {}
    for label, method in [("mean", "mean"), ("tensor", "tensor")]:
        it = [0]

        def cb(_xk, it=it):
            it[0] += 1

        precond = FastDiagonalizationPreconditioner(
            patch, pdm, n_dofs_per_cp=2, fixed_dofs=fixed, law=law, method=method)
        u, info = mf_solve(wq_stiff.apply, ndof, F, fixed, method="gmres",
                           precondition_fn=precond.apply, rtol=1e-10, maxiter=10000,
                           callback=cb, restart=200)
        assert info == 0
        counts[label] = it[0]

    assert counts["tensor"] <= counts["mean"]


def test_operator_diagonal_matches_dense_probe():
    """compute_operator_diagonal() (Gauss-assembled K.diagonal(), used as
    the diag_physical= target for the "S" scaling correction) must match
    K's diagonal from an independent dense-probe of WQMatrixFreeStiffness's
    own apply() -- the two are different quadrature schemes (Gauss vs WQ)
    approximating the same operator, so they should agree closely, not
    exactly."""
    degree, n_levels = 2, 1
    patch, pdm, law = _build_quarter_ring(degree, n_levels)
    wq_stiff = WQMatrixFreeStiffness(patch, law)
    ndof = 2 * patch.n_cp

    diag_gauss = compute_operator_diagonal(patch, law)

    diag_probe = np.zeros(ndof)
    for i in range(ndof):
        e = np.zeros(ndof)
        e[i] = 1.0
        diag_probe[i] = wq_stiff.apply(e)[i]

    assert np.allclose(diag_gauss, diag_probe, rtol=2e-2)
