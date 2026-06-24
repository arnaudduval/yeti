"""
Test stiffness matrix build
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
    MaterialProperties


def stiffness_matrix_legagy(filename):
    """
    Return stiffness matrix computed by legacy function
    Parameters
    ----------
        filename : str
            Path to input files .inp and .NB (no extension, just base name)
    """

    iga_model = IGAparametrization(filename=filename)
    data, row, col, _ = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())

    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    return stiff_side + stiff_side.transpose()


def test_integration_1elt_lin_square_1():
    """
    Test Gauss integration over a single patch, degree 1, dimension 1x1
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/1_elt_lin')

    # u-fastest CP order (nu=2, nv=2): iv outer, iu inner
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])  # flat 0 = (iu=0, iv=0)
    mgr.add_point([1.0, 0.0])  # flat 1 = (iu=1, iv=0)
    mgr.add_point([0.0, 1.0])  # flat 2 = (iu=0, iv=1)
    mgr.add_point([1.0, 1.0])  # flat 3 = (iu=1, iv=1)

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    local_shape = [2, 2]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 2)
    basis_v = IGABasis1D.build(sv, 2)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))
    stiffness_matrix = integrator.integrate_stiffness()

    assert np.allclose(stiffness_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)

def test_integration_1elt_lin_rect():
    """
    Test Gauss integration over a single patch, degree 1, dimension 1x3
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/1_elt_lin_rect')

    # u-fastest CP order (nu=2, nv=2): iv outer, iu inner
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])  # flat 0 = (iu=0, iv=0)
    mgr.add_point([3.0, 0.0])  # flat 1 = (iu=1, iv=0)
    mgr.add_point([0.0, 1.0])  # flat 2 = (iu=0, iv=1)
    mgr.add_point([3.0, 1.0])  # flat 3 = (iu=1, iv=1)

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    local_shape = [2, 2]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 2)
    basis_v = IGABasis1D.build(sv, 2)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))
    stiffness_matrix = integrator.integrate_stiffness()

    assert np.allclose(stiffness_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)

def test_integration_1elt_d2_square_1():
    """
    Test Gauss integration over a single patch, degree 2 in both directions, dimension 1x1
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/1_elt_d2')

    # u-fastest CP order (nu=3, nv=3): iv outer, iu inner
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])   # flat 0 = (iu=0, iv=0)
    mgr.add_point([0.5, 0.0])   # flat 1 = (iu=1, iv=0)
    mgr.add_point([1.0, 0.0])   # flat 2 = (iu=2, iv=0)
    mgr.add_point([0.0, 0.5])   # flat 3 = (iu=0, iv=1)
    mgr.add_point([0.5, 0.5])   # flat 4 = (iu=1, iv=1)
    mgr.add_point([1.0, 0.5])   # flat 5 = (iu=2, iv=1)
    mgr.add_point([0.0, 1.0])   # flat 6 = (iu=0, iv=2)
    mgr.add_point([0.5, 1.0])   # flat 7 = (iu=1, iv=2)
    mgr.add_point([1.0, 1.0])   # flat 8 = (iu=2, iv=2)

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    local_shape = [3, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))
    stiffness_matrix = integrator.integrate_stiffness()

    assert np.allclose(stiffness_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)

def test_integration_1elt_d2_rect():
    """
    Test Gauss integration over a single patch, degree 2 in both directions, dimension 3x1
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/1_elt_d2_rect')

    # u-fastest CP order (nu=3, nv=3): iv outer, iu inner
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])   # flat 0 = (iu=0, iv=0)
    mgr.add_point([1.5, 0.0])   # flat 1 = (iu=1, iv=0)
    mgr.add_point([3.0, 0.0])   # flat 2 = (iu=2, iv=0)
    mgr.add_point([0.0, 0.5])   # flat 3 = (iu=0, iv=1)
    mgr.add_point([1.5, 0.5])   # flat 4 = (iu=1, iv=1)
    mgr.add_point([3.0, 0.5])   # flat 5 = (iu=2, iv=1)
    mgr.add_point([0.0, 1.0])   # flat 6 = (iu=0, iv=2)
    mgr.add_point([1.5, 1.0])   # flat 7 = (iu=1, iv=2)
    mgr.add_point([3.0, 1.0])   # flat 8 = (iu=2, iv=2)

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    local_shape = [3, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))
    stiffness_matrix = integrator.integrate_stiffness()

    assert np.allclose(stiffness_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_integration_2_elements_C0():
    """
    Test Gauss integration over a Patch with 2 elements, C0 continuity.
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/2_elts_C0_d2_rect')

    # u-fastest CP order (nu=5, nv=3): iv outer, iu inner
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])   # flat 0  = (iu=0, iv=0)
    mgr.add_point([1.5, 0.0])   # flat 1  = (iu=1, iv=0)
    mgr.add_point([3.0, 0.0])   # flat 2  = (iu=2, iv=0)
    mgr.add_point([4.5, 0.0])   # flat 3  = (iu=3, iv=0)
    mgr.add_point([6.0, 0.0])   # flat 4  = (iu=4, iv=0)
    mgr.add_point([0.0, 0.5])   # flat 5  = (iu=0, iv=1)
    mgr.add_point([1.5, 0.5])   # flat 6  = (iu=1, iv=1)
    mgr.add_point([3.0, 0.5])   # flat 7  = (iu=2, iv=1)
    mgr.add_point([4.5, 0.5])   # flat 8  = (iu=3, iv=1)
    mgr.add_point([6.0, 0.5])   # flat 9  = (iu=4, iv=1)
    mgr.add_point([0.0, 1.0])   # flat 10 = (iu=0, iv=2)
    mgr.add_point([1.5, 1.0])   # flat 11 = (iu=1, iv=2)
    mgr.add_point([3.0, 1.0])   # flat 12 = (iu=2, iv=2)
    mgr.add_point([4.5, 1.0])   # flat 13 = (iu=3, iv=2)
    mgr.add_point([6.0, 1.0])   # flat 14 = (iu=4, iv=2)

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 0.5, 0.5, 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    local_shape = [5, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))
    stiffness_matrix = integrator.integrate_stiffness()

    assert np.allclose(stiffness_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_integration_2_elements_C1():
    """
    Test Gauss integration over a Patch with 2 elements, C1 continuity.
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/2_elts_C1_d2_rect')

    # u-fastest CP order (nu=4, nv=3): iv outer, iu inner
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])   # flat 0  = (iu=0, iv=0)
    mgr.add_point([1.5, 0.0])   # flat 1  = (iu=1, iv=0)
    mgr.add_point([4.5, 0.0])   # flat 2  = (iu=2, iv=0)
    mgr.add_point([6.0, 0.0])   # flat 3  = (iu=3, iv=0)
    mgr.add_point([0.0, 0.5])   # flat 4  = (iu=0, iv=1)
    mgr.add_point([1.5, 0.5])   # flat 5  = (iu=1, iv=1)
    mgr.add_point([4.5, 0.5])   # flat 6  = (iu=2, iv=1)
    mgr.add_point([6.0, 0.5])   # flat 7  = (iu=3, iv=1)
    mgr.add_point([0.0, 1.0])   # flat 8  = (iu=0, iv=2)
    mgr.add_point([1.5, 1.0])   # flat 9  = (iu=1, iv=2)
    mgr.add_point([4.5, 1.0])   # flat 10 = (iu=2, iv=2)
    mgr.add_point([6.0, 1.0])   # flat 11 = (iu=3, iv=2)

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 0.5, 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    local_shape = [4, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))
    stiffness_matrix = integrator.integrate_stiffness()

    assert np.allclose(stiffness_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


if __name__ == '__main__':
    test_integration_1elt_lin_square_1()
    test_integration_1elt_d2_rect()
    test_integration_1elt_d2_square_1()
    test_integration_1elt_d2_rect()
    test_integration_2_elements_C0()
    test_integration_2_elements_C1()
