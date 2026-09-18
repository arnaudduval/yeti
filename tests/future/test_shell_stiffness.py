"""
Test Kirchhoff-Love shell stiffness matrix build.
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
from yeti_iga.future.bspline import BSpline, BSplineSurface, ControlPointManager, \
    Patch, GlobalDOFManager, PatchDOFManager, IGABasis1D, PatchIntegrator, \
    Material, KirchhoffLoveShellLaw


def stiffness_matrix_legacy(filename):
    """
    Return stiffness matrix computed by the legacy Fortran shell element,
    from an Abaqus-style .inp/.NB fixture (no extension, just base name).
    """

    iga_model = IGAparametrization(filename=filename)
    data, row, col, _ = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())

    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    return stiff_side + stiff_side.transpose()


def relative_frobenius_error(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def test_flat_plate_1_element_degree_3():
    """
    Flat 10x10 unit-square shell patch, degree 3x3, single Bezier element,
    pure B-spline (weights=1). Cross-checked against the legacy Fortran
    shell element (UELMAT3 / 'U3'). Curvature terms are zero here, so this
    exercises the membrane/bending B-matrices and matH but not the 2nd
    derivative NURBS quotient rule.
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    fixture = os.path.join(
        os.path.dirname(script_dir), '..', 'benchs',
        'shell_1_element_load', '1_element_degree_3')
    stiff_legacy = stiffness_matrix_legacy(fixture)

    # u-fastest CP order (4x4 grid), matching the fixture's *Node block.
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

    shell_law = KirchhoffLoveShellLaw(Material(E=210e9, nu=0.3, thickness=0.1))

    integrator = PatchIntegrator(patch, basis_u, basis_v)
    stiffness_matrix = integrator.integrate_shell_stiffness(shell_law)

    assert relative_frobenius_error(
        stiffness_matrix.toarray(), stiff_legacy.toarray()) < 1.e-8


def build_cylinder_patch(radius=2.0, length=3.0):
    """
    Exact NURBS quarter-cylinder shell patch: degree 2 (hoop) x degree 1
    (axial), weight != 1 in the hoop direction (exact circular arc).
    Curved and rational -- exercises the 2nd-derivative NURBS quotient rule
    and the curvature-dependent bending B-matrix, which the flat legacy
    fixtures above cannot (no curved single-patch shell fixture with
    non-planar geometry or weight != 1 exists in benchs/: catenary/shellArch,
    squareShellRoof(Disp) and Tbeam2cplg are all flat panels).
    """

    w1 = np.cos(np.pi / 4)
    hoop_pts = [(radius, 0.0), (radius, radius), (0.0, radius)]
    hoop_w = [1.0, w1, 1.0]

    coords, weights = [], []
    for z in (0.0, length):
        for (x, y), w in zip(hoop_pts, hoop_w):
            coords.append((x, y, z))
            weights.append(w)

    mgr = ControlPointManager(dim=3)
    for c, w in zip(coords, weights):
        mgr.add_point(list(c), w=w)

    ku = np.array([0., 0., 0., 1., 1., 1.])
    kv = np.array([0., 0., 1., 1.])
    su = BSpline(2, ku)
    sv = BSpline(1, kv)
    surf = BSplineSurface(su, sv)
    mapping = list(range(6))
    local_shape = [3, 2]

    dofs_per_cp = [3 for _ in range(mgr.n_points)]
    gdm = GlobalDOFManager(dofs_per_cp)
    pdm = PatchDOFManager(3, mapping, gdm)
    patch = Patch(surf, mgr, mapping, local_shape, pdm)

    basis_u = IGABasis1D.build(su, 4, 2)
    basis_v = IGABasis1D.build(sv, 4, 2)

    return patch, basis_u, basis_v, np.array(coords)


def test_curved_nurbs_cylinder_rigid_body_modes():
    """
    A consistent linear shell stiffness operator must exactly annihilate
    linearized rigid body motions (3 translations + 3 infinitesimal
    rotations): applying such a displacement field produces zero strain
    energy, K @ d == 0. Checked here on a curved, rational (NURBS,
    weight != 1) cylindrical shell patch, which is the only way available
    (in the absence of a curved legacy fixture) to validate the full
    curvature/bendingB/2nd-derivative-NURBS-quotient-rule machinery from
    Phase B/C end-to-end.
    """

    shell_law = KirchhoffLoveShellLaw(Material(E=210e9, nu=0.3, thickness=0.1))
    patch, basis_u, basis_v, coords = build_cylinder_patch()

    integrator = PatchIntegrator(patch, basis_u, basis_v)
    K = integrator.integrate_shell_stiffness(shell_law).toarray()
    k_norm = np.linalg.norm(K)

    n_cp = coords.shape[0]
    n_dof = 3 * n_cp
    center = coords.mean(axis=0)

    rigid_modes = []
    for i in range(3):  # translations
        d = np.zeros(n_dof)
        d[i::3] = 1.0
        rigid_modes.append(d)
    for axis in np.eye(3):  # infinitesimal rotations about the centroid
        d = np.zeros(n_dof)
        for a in range(n_cp):
            d[3 * a:3 * a + 3] = np.cross(axis, coords[a] - center)
        rigid_modes.append(d)

    for d in rigid_modes:
        assert np.linalg.norm(K @ d) / k_norm < 1.e-10
