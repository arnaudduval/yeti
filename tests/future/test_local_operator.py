"""
Test the generic LocalOperator mechanism (PatchIntegrator.integrate_operator /
assemble_operator): subclassing LocalOperator from Python to define just the
per-Gauss-point term to integrate (e.g. B^T*D*B for stiffness), while
PatchIntegrator does the Gauss-point loop and the Jacobian/gradient
computation, exactly like it does internally for the built-in stiffness/mass
kernels.
"""

import os

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (BSpline, BSplineSurface, ControlPointManager,
    Patch, GlobalDOFManager, PatchDOFManager, IGABasis1D, PatchAssembly,
    PatchIntegrator, Material, PlaneStress, LocalOperator)

from test_stiffness import stiffness_matrix_legagy


class PythonStiffnessOperator(LocalOperator):
    """
    Pure-Python re-implementation of PatchIntegrator's built-in B^T*D*B
    stiffness term (see PatchIntegrator.cpp::computeLocalStiffnessContribution).
    PatchIntegrator already did the Gauss-point loop and the Jacobian
    inversion -- only the algebraic term itself is written here.
    """

    def __init__(self, E, nu):
        super().__init__()
        self.E = E
        self.nu = nu

    def compute_integrand(self, R, dRdx, dRdy):
        nb_loc = len(R)
        factor = self.E / (1.0 - self.nu ** 2)
        D = np.array([
            [factor, factor * self.nu, 0.0],
            [factor * self.nu, factor, 0.0],
            [0.0, 0.0, factor * (1.0 - self.nu) / 2.0],
        ])

        B = np.zeros((3, 2 * nb_loc))
        B[0, 0::2] = dRdx
        B[1, 1::2] = dRdy
        B[2, 0::2] = dRdy
        B[2, 1::2] = dRdx

        return B.T @ D @ B


class LumpedMassOperator(LocalOperator):
    """
    Row-sum ("HRZ") lumped mass term -- deliberately NOT one of
    PatchIntegrator's built-in kernels (only the consistent mass matrix is
    implemented in C++), to demonstrate genuine extensibility rather than
    just reproducing an existing one.

    Row-summing a consistent mass term outer(R, R) at a single Gauss point
    gives R itself (partition of unity: the span's nb_loc active basis
    functions already sum to 1 there), so lumping reduces to a per-Gauss-
    point diagonal term -- no need to wait for the full span to be
    integrated before lumping.
    """

    def __init__(self, rho):
        super().__init__()
        self.rho = rho

    def compute_integrand(self, R, dRdx, dRdy):
        nb_loc = len(R)
        diag = self.rho * np.asarray(R)
        M_loc = np.zeros((2 * nb_loc, 2 * nb_loc))
        M_loc[0::2, 0::2] = np.diag(diag)
        M_loc[1::2, 1::2] = np.diag(diag)
        return M_loc


def _build_1elt_d2_rect_patch():
    """Same fixture/layout as test_stiffness.py::test_integration_1elt_d2_rect."""
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
    return patch, basis_u, basis_v


def test_integrate_operator_reproduces_stiffness():
    """
    integrate_operator() with a Python LocalOperator that re-implements the
    stiffness kernel must reproduce the legacy-validated stiffness matrix,
    exactly like integrate_stiffness() does.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/1_elt_d2_rect')

    patch, basis_u, basis_v = _build_1elt_d2_rect_patch()
    integrator = PatchIntegrator(patch, basis_u, basis_v, PlaneStress(Material(E=210000, nu=0.3)))

    op = PythonStiffnessOperator(210000, 0.3)
    operator_matrix = integrator.integrate_operator(op)
    stiffness_matrix = integrator.integrate_stiffness()

    np.testing.assert_allclose(
        operator_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)
    np.testing.assert_allclose(
        operator_matrix.toarray(), stiffness_matrix.toarray(), rtol=1.e-8, atol=1.e-8)


def test_integrate_operator_lumped_mass_conserves_total_mass():
    """
    A lumped mass kernel has no built-in counterpart (only the consistent
    mass matrix is implemented in C++) -- check it behaves sanely: diagonal,
    and total mass (sum of the ux-block diagonal) equal to rho * area
    (partition of unity: sum_a R_a(xi) == 1 everywhere).
    """
    patch, basis_u, basis_v = _build_1elt_d2_rect_patch()
    rho = 7800.0
    integrator = PatchIntegrator(patch, basis_u, basis_v, PlaneStress(Material(E=210000, nu=0.3)))

    mass_matrix = integrator.integrate_operator(LumpedMassOperator(rho)).toarray()

    off_diag = mass_matrix - np.diag(np.diag(mass_matrix))
    assert np.allclose(off_diag, 0.0, atol=1.e-10)

    area = 3.0 * 1.0  # patch spans x in [0, 3], y in [0, 1]
    ux_dofs = mass_matrix[0::2, 0::2]
    assert np.isclose(np.trace(ux_dofs), rho * area, rtol=1.e-6)


def test_assemble_operator_domain_split_matches_single_patch():
    """
    assemble_operator() on a 2-patch domain split must reproduce the same
    global matrix as integrate_operator() on the legacy-validated single
    patch -- same self-consistency check as
    test_multipatch_stiffness.py::test_assemble_domain_split_matches_single_patch,
    but for a custom Python operator instead of the built-in kernel.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/2_elts_C0_d2_rect')

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

    dofs_per_control_point = [2] * mgr.n_points
    global_dof_manager = GlobalDOFManager(dofs_per_control_point)

    su_left = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv_left = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    mapping_left = [0, 1, 2, 5, 6, 7, 10, 11, 12]
    dof_manager_left = PatchDOFManager(2, mapping_left, global_dof_manager)
    patch_left = Patch(BSplineSurface(su_left, sv_left), mgr, mapping_left, [3, 3],
                        dof_manager_left)

    su_right = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv_right = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    mapping_right = [2, 3, 4, 7, 8, 9, 12, 13, 14]
    dof_manager_right = PatchDOFManager(2, mapping_right, global_dof_manager)
    patch_right = Patch(BSplineSurface(su_right, sv_right), mgr, mapping_right, [3, 3],
                         dof_manager_right)

    assembly = PatchAssembly()
    assembly.add_patch(patch_left)
    assembly.add_patch(patch_right)
    assembly.detect_shared_control_points()

    op_left = PythonStiffnessOperator(210000, 0.3)
    op_right = PythonStiffnessOperator(210000, 0.3)
    operator_matrix = PatchIntegrator.assemble_operator(assembly, [op_left, op_right])

    assert operator_matrix.shape == stiff_legacy.shape
    np.testing.assert_allclose(
        operator_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


if __name__ == '__main__':
    test_integrate_operator_reproduces_stiffness()
    test_integrate_operator_lumped_mass_conserves_total_mass()
    test_assemble_operator_domain_split_matches_single_patch()
