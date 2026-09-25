# -*- coding: utf-8 -*-
"""
Test de validation du solve statique (elasticite lineaire).

Regression pour MechanicalModel.solve_linearized_system(): avant correction,
cette methode resolvait toujours l'operateur de masse (compute_mf_mass), meme
pour un modele purement statique -- appeler solve_linearized_system() sans
scalar_coefs (comme le fait StaticElastoPlasticity/IncrementalElastoPlasticity)
resolvait donc la mauvaise equation. Ce test verifie qu'un tel appel resout
bien la raideur (K @ u = F), avec un resultat qui coincide avec une resolution
directe independante (sonde dense de compute_mf_stiffness + np.linalg.solve).
"""

import numpy as np

from yeti_iga.pymfiga.iga.boundary import (
    BoundaryCondition,
    DOF,
    ParametricDirection,
    BoundarySide,
)
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.single_model import MechanicalModel
from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.physics import StaticElastoPlasticity
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver


DEGREE = 2
NBEL = 4

E = 2.1e11
NU = 0.3

QUADCLASS = "wq"
QUADTYPE = "2"


def build_model():
    geometry = GeomdlGenerator(
        filename="quarter_annulus",
        geo_args={"degree": DEGREE, "nbel": NBEL},
    ).export_geometry()

    patch = SinglePatch(geometry, quadclass=QUADCLASS, quadtype=QUADTYPE)
    patch.generate()

    material = LinearElasticity({"elastic_modulus": E, "poisson_ratio": NU})

    # Symmetry BCs: theta=0 edge (ETA-min) blocks UY, theta=90deg edge
    # (ETA-max) blocks UX -- same convention as future/fast_diagonalization.py
    # and the yeti-iga notebooks' quarter-ring benchmark.
    boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX, DOF.UY))
    boundary.add_constraint(
        constraint_info=[
            {"direction": ParametricDirection.ETA, "face": BoundarySide.MIN, "dofs": (DOF.UY,)},
            {"direction": ParametricDirection.ETA, "face": BoundarySide.MAX, "dofs": (DOF.UX,)},
        ],
        constraint_type="dirichlet",
    )

    model = MechanicalModel(material, patch, boundary)
    return model


def surface_force(args):
    position = args["position"]
    force = np.zeros_like(position)
    force[1, :] = 1.0  # constant unit traction in y
    return force


def _assemble_force(model):
    return model.assemble_surface_force(
        {DOF.ALL: ({"direction": ParametricDirection.XI, "face": BoundarySide.MAX}, surface_force)}
    )


def _direct_reference(model, F, free):
    ndof = len(F)
    n_free = len(free)
    K_free = np.zeros((n_free, n_free))
    for j, dof_j in enumerate(free):
        e = np.zeros(ndof)
        e[dof_j] = 1.0
        K_free[:, j] = model.compute_mf_stiffness(e)[free]
    u_ref = np.zeros(ndof)
    u_ref[free] = np.linalg.solve(K_free, F[free])
    return u_ref


def test_solve_linearized_system_defaults_to_static_stiffness():
    model = build_model()
    F = _assemble_force(model)
    free, fixed = model.get_free_and_constraint_nodes()

    linear_solver = LinearSolver(tolerance=1e-12, maxiters=2000, linear_type="gmres")

    # No scalar_coefs passed -- must default to pure stiffness (0, 1), not
    # silently solve the mass operator.
    u = model.solve_linearized_system(F, linear_solver_backend=linear_solver)

    u_ref = _direct_reference(model, F, free)

    assert np.allclose(u[free], u_ref[free], rtol=1e-6, atol=1e-9)
    assert np.allclose(u[fixed], 0.0, atol=1e-12)


def test_solve_linearized_system_without_preconditioner_also_matches():
    model = build_model()
    F = _assemble_force(model)
    free, _ = model.get_free_and_constraint_nodes()

    linear_solver = LinearSolver(tolerance=1e-12, maxiters=5000, linear_type="gmres")
    u = model.solve_linearized_system(
        F, linear_solver_backend=linear_solver, use_preconditioner=False
    )

    u_ref = _direct_reference(model, F, free)
    assert np.allclose(u[free], u_ref[free], rtol=1e-6, atol=1e-9)


def test_static_elastoplasticity_physics_solves_correctly():
    """End-to-end check via the actual high-level API that was silently
    broken (StaticElastoPlasticity never passes scalar_coefs either)."""
    model = build_model()
    F = _assemble_force(model)
    free, _ = model.get_free_and_constraint_nodes()

    physics = StaticElastoPlasticity(
        tolerance_linear=1e-12,
        maxiters_linear=2000,
        tolerance_nonlinear=1e-8,
        maxiters_nonlinear=5,
    )
    u = np.zeros_like(F)
    physics.solve(model, u, F)

    u_ref = _direct_reference(model, F, free)
    assert np.allclose(u[free], u_ref[free], rtol=1e-6, atol=1e-9)
