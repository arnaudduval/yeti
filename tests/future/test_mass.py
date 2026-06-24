"""
Test mass matrix build
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.massmtrx import build_cmassmatrix
from yeti_iga.future.bspline import BSpline, BSplineSurface, ControlPointManager, \
    Patch, GlobalDOFManager, PatchDOFManager, IGABasis1D, PatchIntegrator, \
    MaterialProperties

RHO = 7800.0


def mass_matrix_legagy(filename, rho=RHO):
    """
    Return mass matrix computed by legacy function

    Parameters
    ----------
        filename : str
            Path to input files .inp and .NB (no extension, just base name)
        rho : float
            Mass density. The .inp fixtures used here carry no *Density
            section (density defaults to 0 on read), so it is set directly
            on the loaded model before calling get_inputs4massmat().
    """

    iga_model = IGAparametrization(filename=filename)
    iga_model._MATERIAL_PROPERTIES[2, :] = rho
    data, row, col = build_cmassmatrix(*iga_model.get_inputs4massmat())

    mass_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    return mass_side + mass_side.transpose()


def test_mass_1elt_lin_square_1():
    """
    Test Gauss integration over a single patch, degree 1, dimension 1x1
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    mass_legacy = mass_matrix_legagy(f'{script_dir}/1_elt_lin')

    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.0, 1.0])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    local_shape = [2, 2]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    basis_u = IGABasis1D.build(su, 2)
    basis_v = IGABasis1D.build(sv, 2)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3, rho=RHO))
    mass_matrix = integrator.integrate_mass()

    assert np.allclose(mass_matrix.toarray(), mass_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_mass_1elt_lin_rect():
    """
    Test Gauss integration over a single patch, degree 1, dimension 1x3
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    mass_legacy = mass_matrix_legagy(f'{script_dir}/1_elt_lin_rect')

    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([3.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([3.0, 1.0])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    local_shape = [2, 2]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    basis_u = IGABasis1D.build(su, 2)
    basis_v = IGABasis1D.build(sv, 2)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3, rho=RHO))
    mass_matrix = integrator.integrate_mass()

    assert np.allclose(mass_matrix.toarray(), mass_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_mass_1elt_d2_square_1():
    """
    Test Gauss integration over a single patch, degree 2 in both directions, dimension 1x1
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    mass_legacy = mass_matrix_legagy(f'{script_dir}/1_elt_d2')

    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([0.5, 0.0])
    mgr.add_point([1.0, 0.0])
    mgr.add_point([0.0, 0.5])
    mgr.add_point([0.5, 0.5])
    mgr.add_point([1.0, 0.5])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([0.5, 1.0])
    mgr.add_point([1.0, 1.0])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    local_shape = [3, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3, rho=RHO))
    mass_matrix = integrator.integrate_mass()

    assert np.allclose(mass_matrix.toarray(), mass_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_mass_1elt_d2_rect():
    """
    Test Gauss integration over a single patch, degree 2 in both directions, dimension 3x1
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    mass_legacy = mass_matrix_legagy(f'{script_dir}/1_elt_d2_rect')

    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.5, 0.0])
    mgr.add_point([3.0, 0.0])
    mgr.add_point([0.0, 0.5])
    mgr.add_point([1.5, 0.5])
    mgr.add_point([3.0, 0.5])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.5, 1.0])
    mgr.add_point([3.0, 1.0])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    local_shape = [3, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3, rho=RHO))
    mass_matrix = integrator.integrate_mass()

    assert np.allclose(mass_matrix.toarray(), mass_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_mass_2_elements_C0():
    """
    Test Gauss integration over a Patch with 2 elements, C0 continuity.
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    mass_legacy = mass_matrix_legagy(f'{script_dir}/2_elts_C0_d2_rect')

    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.5, 0.0])
    mgr.add_point([3.0, 0.0])
    mgr.add_point([4.5, 0.0])
    mgr.add_point([6.0, 0.0])
    mgr.add_point([0.0, 0.5])
    mgr.add_point([1.5, 0.5])
    mgr.add_point([3.0, 0.5])
    mgr.add_point([4.5, 0.5])
    mgr.add_point([6.0, 0.5])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.5, 1.0])
    mgr.add_point([3.0, 1.0])
    mgr.add_point([4.5, 1.0])
    mgr.add_point([6.0, 1.0])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 0.5, 0.5, 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    local_shape = [5, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3, rho=RHO))
    mass_matrix = integrator.integrate_mass()

    assert np.allclose(mass_matrix.toarray(), mass_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_mass_2_elements_C1():
    """
    Test Gauss integration over a Patch with 2 elements, C1 continuity.
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    mass_legacy = mass_matrix_legagy(f'{script_dir}/2_elts_C1_d2_rect')

    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.5, 0.0])
    mgr.add_point([4.5, 0.0])
    mgr.add_point([6.0, 0.0])
    mgr.add_point([0.0, 0.5])
    mgr.add_point([1.5, 0.5])
    mgr.add_point([4.5, 0.5])
    mgr.add_point([6.0, 0.5])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.5, 1.0])
    mgr.add_point([4.5, 1.0])
    mgr.add_point([6.0, 1.0])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 0.5, 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    local_shape = [4, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3, rho=RHO))
    mass_matrix = integrator.integrate_mass()

    assert np.allclose(mass_matrix.toarray(), mass_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_mass_requires_rho():
    """
    integrate_mass() must fail loudly (not silently return a zero/wrong
    matrix) if rho was not set on the material.
    """
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.0, 1.0])

    dof_manager = GlobalDOFManager([2] * mgr.n_points)
    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping = [0, 1, 2, 3]
    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [2, 2], dof_manager_patch)

    basis_u = IGABasis1D.build(su, 2)
    basis_v = IGABasis1D.build(sv, 2)

    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))
    try:
        integrator.integrate_mass()
        assert False, "expected an exception when rho is not set"
    except Exception:
        pass


if __name__ == '__main__':
    test_mass_1elt_lin_square_1()
    test_mass_1elt_lin_rect()
    test_mass_1elt_d2_square_1()
    test_mass_1elt_d2_rect()
    test_mass_2_elements_C0()
    test_mass_2_elements_C1()
    test_mass_requires_rho()
