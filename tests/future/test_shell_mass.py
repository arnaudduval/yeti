"""
Test Kirchhoff-Love shell mass matrix build.

No shell fixture in benchs/ defines a material density (all are static
analyses), so there is no legacy mass matrix to cross-check against.
Instead, this validates the consistent mass matrix against the analytical
total mass of a flat plate: for a partition-of-unity basis (B-spline or
NURBS), sum_{a,b} R_a * R_b integrated over the patch equals the patch
area, so summing one full translational block of M must equal rho * t * Area.
"""

import numpy as np

from yeti_iga.future.bspline import BSpline, BSplineSurface, ControlPointManager, \
    Patch, GlobalDOFManager, PatchDOFManager, IGABasis1D, PatchIntegrator, \
    Material, KirchhoffLoveShellLaw


def test_flat_plate_total_mass():
    """
    Flat 10x10 unit-square shell patch, degree 3x3, single Bezier element
    (same geometry as test_shell_stiffness.test_flat_plate_1_element_degree_3,
    with a synthetic density added). total_mass = rho * thickness * area.
    """

    coords = [
        (0.0, 0.0, 0.0), (3.33, 0.0, 0.0), (6.67, 0.0, 0.0), (10.0, 0.0, 0.0),
        (0.0, 3.33, 0.0), (3.33, 3.33, 0.0), (6.67, 3.33, 0.0), (10.0, 3.33, 0.0),
        (0.0, 6.67, 0.0), (3.33, 6.67, 0.0), (6.67, 6.67, 0.0), (10.0, 6.67, 0.0),
        (0.0, 10.0, 0.0), (3.33, 10.0, 0.0), (6.67, 10.0, 0.0), (10.0, 10.0, 0.0),
    ]
    mgr = ControlPointManager(dim=3)
    for c in coords:
        mgr.add_point(list(c))

    kv = np.array([0., 0., 0., 0., 1., 1., 1., 1.])
    su = BSpline(3, kv)
    sv = BSpline(3, kv)
    surf = BSplineSurface(su, sv)
    mapping = list(range(16))
    local_shape = [4, 4]

    dofs_per_cp = [3 for _ in range(mgr.n_points)]
    gdm = GlobalDOFManager(dofs_per_cp)
    pdm = PatchDOFManager(3, mapping, gdm)
    patch = Patch(surf, mgr, mapping, local_shape, pdm)

    basis_u = IGABasis1D.build(su, 4, 2)
    basis_v = IGABasis1D.build(sv, 4, 2)

    rho = 7850.0
    thickness = 0.1
    shell_law = KirchhoffLoveShellLaw(Material(E=210e9, nu=0.3, rho=rho, thickness=thickness))

    integrator = PatchIntegrator(patch, basis_u, basis_v)
    M = integrator.integrate_shell_mass(shell_law).toarray()

    area = 10.0 * 10.0
    expected_mass = rho * thickness * area

    x_dofs = slice(0, None, 3)
    total_mass = M[x_dofs, x_dofs].sum()

    assert abs(total_mass - expected_mass) / expected_mass < 1.e-10

    # The 3 translational directions must carry identical, uncoupled mass.
    y_dofs = slice(1, None, 3)
    z_dofs = slice(2, None, 3)
    assert np.allclose(M[x_dofs, x_dofs], M[y_dofs, y_dofs])
    assert np.allclose(M[x_dofs, x_dofs], M[z_dofs, z_dofs])
    assert np.allclose(M[x_dofs, y_dofs], 0.0)
    assert np.allclose(M[x_dofs, z_dofs], 0.0)
    assert np.allclose(M[y_dofs, z_dofs], 0.0)
