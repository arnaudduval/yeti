"""
Test Kirchhoff-Love shell distributed surface load build.
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


def test_flat_plate_load_vs_legacy():
    """
    benchs/squareShellRoof/squareShellPlate defines a *Dload I1.EltToLoad,
    U66, -500. (legacy "snow load": a global-Z force per unit area scaled by
    the local normal's Z component). On this flat patch the normal is
    constant, so it reduces exactly to a constant vertical vector -- cross-
    checked here against integrate_shell_surface_load(direction=(0,0,1),
    magnitude=-500).

    Every one of this fixture's 4 control points is a corner, so
    IGAparametrization._update_dof_info() hits its nb_dof_free == 0 edge
    case (see test_shell_stiffness.py's analogous stiffness test) --
    ind_dof_free is overridden the same way to build the raw RHS vector.
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    fixture = os.path.join(
        os.path.dirname(script_dir), '..', 'benchs',
        'squareShellRoof', 'squareShellPlate')

    iga_model = IGAparametrization(filename=fixture)
    if iga_model._nb_dof_free == 0:
        iga_model._ind_dof_free = np.arange(1, iga_model._nb_dof_tot + 1)
        iga_model._nb_dof_free = iga_model._nb_dof_tot
    _, _, _, F_legacy = build_stiffmatrix(*iga_model.get_inputs4system_elemStorage())

    coords = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (10.0, 10.0, 0.0)]
    mgr = ControlPointManager(dim=3)
    for c in coords:
        mgr.add_point(list(c))

    kv = np.array([0., 0., 1., 1.])
    su = BSpline(1, kv)
    sv = BSpline(1, kv)
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    local_shape = [2, 2]

    dofs_per_cp = [3 for _ in range(mgr.n_points)]
    gdm = GlobalDOFManager(dofs_per_cp)
    pdm = PatchDOFManager(3, mapping, gdm)
    patch = Patch(surf, mgr, mapping, local_shape, pdm)

    basis_u = IGABasis1D.build(su, 2, 2)
    basis_v = IGABasis1D.build(sv, 2, 2)

    integrator = PatchIntegrator(patch, basis_u, basis_v)
    F_future = integrator.integrate_shell_surface_load(np.array([0., 0., 1.]), -500.0)

    assert np.linalg.norm(F_future - F_legacy) / np.linalg.norm(F_future) < 1.e-8


def test_flat_plate_load_resultant():
    """
    Analytical check on the degree-3 flat plate from test_shell_stiffness.py:
    for a partition-of-unity basis, integrating a constant force-per-unit-area
    vector against R_a over the whole patch must give a total resultant force
    of magnitude * direction * area, regardless of direction/degree/mesh.
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

    integrator = PatchIntegrator(patch, basis_u, basis_v)
    direction = np.array([0.0, 0.0, -1.0])
    magnitude = 250.0
    F = integrator.integrate_shell_surface_load(direction, magnitude)

    area = 10.0 * 10.0
    expected_resultant = magnitude * direction * area

    resultant = np.array([F[0::3].sum(), F[1::3].sum(), F[2::3].sum()])
    assert np.allclose(resultant, expected_resultant, rtol=1.e-10)
