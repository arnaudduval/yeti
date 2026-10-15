"""
Test the 3D solid (BSplineVolume) consistent mass kernel
(computeLocalMassContribution3D / integrate_mass()).
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineVolume, Patch, ControlPointManager,
    GlobalDOFManager, PatchDOFManager, IGABasis1D, Material, Elasticity3D,
    PatchIntegrator,
)

RHO = 7800.0


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


def test_mass3d_unit_field_gives_rho_times_volume():
    """
    u^T M u == rho * Volume for a constant (rigid-translation) unit field --
    partition of unity: sum_a R_a(xi) == 1 everywhere, so the diagonal
    translational mass is exactly the total mass.
    """
    degree, nbel = 2, 2
    Lx, Ly, Lz = 3.0, 1.5, 2.0
    patch, coords = _build_box_patch(degree, nbel, Lx, Ly, Lz)

    law = Elasticity3D(Material(E=210000.0, nu=0.3, rho=RHO))
    basis_u = IGABasis1D.build(patch.tensor.components[0], degree + 1)
    basis_v = IGABasis1D.build(patch.tensor.components[1], degree + 1)
    basis_w = IGABasis1D.build(patch.tensor.components[2], degree + 1)
    integrator = PatchIntegrator(patch, basis_u, basis_v, basis_w, law)

    M = integrator.integrate_mass()

    n_cp = coords.shape[0]
    u = np.tile([1.0, 0.0, 0.0], n_cp)
    volume = Lx * Ly * Lz
    np.testing.assert_allclose(u @ (M @ u), RHO * volume, rtol=1.e-8)


def test_mass3d_symmetric():
    degree, nbel = 2, 1
    patch, _ = _build_box_patch(degree, nbel, 2.0, 1.0, 1.0)

    law = Elasticity3D(Material(E=210000.0, nu=0.3, rho=RHO))
    basis_u = IGABasis1D.build(patch.tensor.components[0], degree + 1)
    basis_v = IGABasis1D.build(patch.tensor.components[1], degree + 1)
    basis_w = IGABasis1D.build(patch.tensor.components[2], degree + 1)
    integrator = PatchIntegrator(patch, basis_u, basis_v, basis_w, law)

    M = integrator.integrate_mass().toarray()
    np.testing.assert_allclose(M, M.T, atol=1.e-10)


def test_mass3d_requires_rho():
    degree, nbel = 1, 1
    patch, _ = _build_box_patch(degree, nbel, 1.0, 1.0, 1.0)

    law = Elasticity3D(Material(E=210000.0, nu=0.3))  # rho not set
    basis_u = IGABasis1D.build(patch.tensor.components[0], degree + 1)
    basis_v = IGABasis1D.build(patch.tensor.components[1], degree + 1)
    basis_w = IGABasis1D.build(patch.tensor.components[2], degree + 1)
    integrator = PatchIntegrator(patch, basis_u, basis_v, basis_w, law)

    try:
        integrator.integrate_mass()
        assert False, "expected an exception when rho is not set"
    except Exception:
        pass


if __name__ == '__main__':
    test_mass3d_unit_field_gives_rho_times_volume()
    test_mass3d_symmetric()
    test_mass3d_requires_rho()
