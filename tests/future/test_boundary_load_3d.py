"""
Test the 3D boundary load kernel (computeLocalBoundaryLoadContribution3D /
integrate_boundary_load_3d()): a 3D boundary load integrates over a FACE (2
varying parametric directions), needing the surface Jacobian
dA = |T_p x T_q| (cross product of the two physical tangent vectors) instead
of the 2D case's curve-tangent norm ds = |T|.
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineVolume, Patch, ControlPointManager,
    GlobalDOFManager, PatchDOFManager, IGABasis1D, Material, Elasticity3D,
    PatchIntegrator, Traction3D, ConstantTraction3D,
)


def _build_box_patch(degree, nbel, Lx, Ly, Lz):
    """Non-rational BSplineVolume box, x in [0,Lx], y in [0,Ly], z in [0,Lz],
    Greville-abscissae control points (exact for this affine map)."""
    kv_interior = np.linspace(0.0, 1.0, nbel + 1)[1:-1]
    kv = np.concatenate([[0.0] * (degree + 1), kv_interior, [1.0] * (degree + 1)])
    splines = [BSpline(degree, kv) for _ in range(3)]
    tensor = BSplineVolume(*splines)

    n1 = nbel + degree
    n_cp = n1 ** 3
    n_basis = len(kv) - degree - 1
    gr = np.array([np.mean(kv[i+1:i+degree+1]) for i in range(n_basis)])

    mgr = ControlPointManager(dim=3)
    coords = np.zeros((n_cp, 3))
    idx = 0
    for kw in range(n1):
        for jv in range(n1):
            for iu in range(n1):
                x, y, z = gr[iu] * Lx, gr[jv] * Ly, gr[kw] * Lz
                coords[idx] = (x, y, z)
                mgr.add_point([x, y, z])
                idx += 1

    patch0 = Patch(tensor, mgr, list(range(n_cp)), [n1, n1, n1])
    gdm = GlobalDOFManager([3] * n_cp)
    pdm = PatchDOFManager(3, list(range(n_cp)), gdm)
    patch = Patch(patch0.tensor, patch0.cp_manager, list(patch0.global_indices),
                  list(patch0.local_shape), pdm)
    return patch, coords


def _build_integrator(patch, degree, law):
    # PatchIntegrator's pybind11 constructor uses py::keep_alive to keep
    # patch/basis_u/v/w alive for as long as the integrator is (see
    # test_patch_integrator_lifetime.py), so basis_u/v/w may safely go out
    # of scope here.
    basis_u = IGABasis1D.build(patch.tensor.components[0], degree + 1)
    basis_v = IGABasis1D.build(patch.tensor.components[1], degree + 1)
    basis_w = IGABasis1D.build(patch.tensor.components[2], degree + 1)
    return PatchIntegrator(patch, basis_u, basis_v, basis_w, law)


def test_integrate_boundary_load_3d_resultant_force_uniform_traction():
    """
    A constant traction (0, 0, -p) on a face of area A must produce a load
    vector whose z-components sum to -p*A (partition of unity), x/y
    components zero, and zero entries away from the loaded face.
    """
    degree, nbel = 2, 1
    Lx, Ly, Lz = 3.0, 2.0, 1.0
    patch, coords = _build_box_patch(degree, nbel, Lx, Ly, Lz)
    law = Elasticity3D(Material(E=210000.0, nu=0.3))
    integrator = _build_integrator(patch, degree, law)

    p = 1000.0
    traction = ConstantTraction3D(np.array([0.0, 0.0, -p]))

    # Top face (z=Lz, direction=2, side=1), area A = Lx*Ly.
    F = integrator.integrate_boundary_load_3d(2, 1, traction)

    A = Lx * Ly
    np.testing.assert_allclose(F[0::3].sum(), 0.0, atol=1.e-8)
    np.testing.assert_allclose(F[1::3].sum(), 0.0, atol=1.e-8)
    np.testing.assert_allclose(F[2::3].sum(), -p * A, rtol=1.e-10)

    loaded_cps = set(patch.boundary_control_points(2, 1))
    n_cp = coords.shape[0]
    for cp in range(n_cp):
        if cp not in loaded_cps:
            dofs = patch.dof_manager.get_global_dof_indices(cp)
            for d in dofs:
                assert F[d] == 0.0


def test_integrate_boundary_load_3d_linear_traction_matches_moment():
    """
    A traction linear in y, t = (M*y/Iz, 0, 0), on the face x=L of a beam
    cross-section [-h/2,h/2] x [-b/2,b/2] must produce zero net force (odd
    in y, symmetric section) and a virtual-work moment of exactly M under
    the rigid-rotation virtual field v = (-y, x, 0) -- validates the
    cross-product surface Jacobian (dA = |T_p x T_q|) against an
    independent, geometry-only expectation (int y^2 dA = Iz).
    """
    degree, nbel = 2, 1
    L, h, b = 4.0, 1.0, 1.0
    Iz = b * h**3 / 12.0
    M = 1.0

    kv_interior = np.linspace(0.0, 1.0, nbel + 1)[1:-1]
    kv = np.concatenate([[0.0] * (degree + 1), kv_interior, [1.0] * (degree + 1)])
    splines = [BSpline(degree, kv) for _ in range(3)]
    tensor = BSplineVolume(*splines)
    n1 = nbel + degree
    n_cp = n1 ** 3
    n_basis = len(kv) - degree - 1
    gr = np.array([np.mean(kv[i+1:i+degree+1]) for i in range(n_basis)])

    mgr = ControlPointManager(dim=3)
    coords = np.zeros((n_cp, 3))
    idx = 0
    for kw in range(n1):
        for jv in range(n1):
            for iu in range(n1):
                x = gr[iu] * L
                y = (gr[jv] - 0.5) * h
                z = (gr[kw] - 0.5) * b
                coords[idx] = (x, y, z)
                mgr.add_point([x, y, z])
                idx += 1

    patch0 = Patch(tensor, mgr, list(range(n_cp)), [n1, n1, n1])
    gdm = GlobalDOFManager([3] * n_cp)
    pdm = PatchDOFManager(3, list(range(n_cp)), gdm)
    patch = Patch(patch0.tensor, patch0.cp_manager, list(patch0.global_indices),
                  list(patch0.local_shape), pdm)

    law = Elasticity3D(Material(E=210000.0, nu=0.3))
    integrator = _build_integrator(patch, degree, law)

    class BendingTraction(Traction3D):
        def evaluate(self, pt):
            _, y, _ = pt
            return np.array([M * y / Iz, 0.0, 0.0])

    F = integrator.integrate_boundary_load_3d(0, 1, BendingTraction())

    np.testing.assert_allclose(F[0::3].sum(), 0.0, atol=1.e-10)
    np.testing.assert_allclose(F[1::3].sum(), 0.0, atol=1.e-10)
    np.testing.assert_allclose(F[2::3].sum(), 0.0, atol=1.e-10)

    # Rigid rotation v = omega x r about the z axis (omega=1): v_x=-y, v_y=x.
    # Virtual work = int t.v dA = int (M*y/Iz)*(-y) dA = -(M/Iz)*int y^2 dA
    # = -(M/Iz)*Iz = -M.
    v = np.zeros(3 * n_cp)
    for i in range(n_cp):
        x, y, z = coords[i]
        v[3*i + 0] = -y
        v[3*i + 1] = x
    np.testing.assert_allclose(F @ v, -M, rtol=1.e-10)


if __name__ == '__main__':
    test_integrate_boundary_load_3d_resultant_force_uniform_traction()
    test_integrate_boundary_load_3d_linear_traction_matches_moment()
