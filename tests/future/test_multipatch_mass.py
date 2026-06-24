"""
Test multipatch mass matrix assembly (PatchIntegrator.assemble_mass)
"""

import os

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (BSpline, BSplineSurface, ControlPointManager,
    Patch, GlobalDOFManager, PatchDOFManager, PatchAssembly, PatchIntegrator,
    MaterialProperties)

from test_mass import mass_matrix_legagy, RHO


def test_assemble_mass_domain_split_matches_single_patch():
    """
    Splitting a single, legacy-validated patch into two conforming patches
    sharing the interior boundary must not change the assembled mass:
    PatchIntegrator.assemble_mass() on the 2-patch assembly must reproduce
    the same global mass matrix as the single-patch legacy reference for
    '2_elts_C0_d2_rect' (already validated in test_mass.py).
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    mass_legacy = mass_matrix_legagy(f'{script_dir}/2_elts_C0_d2_rect')

    # Same 15 control points as test_mass_2_elements_C0 (nu=5, nv=3,
    # u-fastest: iv outer, iu inner).
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])   # 0
    mgr.add_point([1.5, 0.0])   # 1
    mgr.add_point([3.0, 0.0])   # 2  (shared column)
    mgr.add_point([4.5, 0.0])   # 3
    mgr.add_point([6.0, 0.0])   # 4
    mgr.add_point([0.0, 0.5])   # 5
    mgr.add_point([1.5, 0.5])   # 6
    mgr.add_point([3.0, 0.5])   # 7  (shared column)
    mgr.add_point([4.5, 0.5])   # 8
    mgr.add_point([6.0, 0.5])   # 9
    mgr.add_point([0.0, 1.0])   # 10
    mgr.add_point([1.5, 1.0])   # 11
    mgr.add_point([3.0, 1.0])   # 12 (shared column)
    mgr.add_point([4.5, 1.0])   # 13
    mgr.add_point([6.0, 1.0])   # 14

    dofs_per_control_point = [2] * mgr.n_points
    global_dof_manager = GlobalDOFManager(dofs_per_control_point)

    # Left patch: columns iu=0,1,2
    su_left = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv_left = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    mapping_left = [0, 1, 2, 5, 6, 7, 10, 11, 12]
    dof_manager_left = PatchDOFManager(2, mapping_left, global_dof_manager)
    patch_left = Patch(BSplineSurface(su_left, sv_left), mgr, mapping_left, [3, 3],
                        dof_manager_left)

    # Right patch: columns iu=2,3,4 (column 2 shared with the left patch)
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

    shared_map = assembly.get_shared_control_points_map()
    assert set(shared_map.keys()) == {2, 7, 12}

    material = MaterialProperties(210000, 0.3, rho=RHO)
    mass_matrix = PatchIntegrator.assemble_mass(assembly, [material, material])

    assert mass_matrix.shape == mass_legacy.shape
    np.testing.assert_allclose(
        mass_matrix.toarray(), mass_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_assemble_mass_requires_rho_on_every_material():
    """
    assemble_mass() must fail loudly if any patch's material has no rho set.
    """
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.0, 1.0])
    mgr.add_point([2.0, 0.0])
    mgr.add_point([2.0, 1.0])

    dofs_per_control_point = [2] * mgr.n_points
    global_dof_manager = GlobalDOFManager(dofs_per_control_point)

    su_a = BSpline(1, np.array([0., 0., 1., 1.]))
    sv_a = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping_a = [0, 1, 2, 3]
    dof_manager_a = PatchDOFManager(2, mapping_a, global_dof_manager)
    patch_a = Patch(BSplineSurface(su_a, sv_a), mgr, mapping_a, [2, 2], dof_manager_a)

    su_b = BSpline(1, np.array([0., 0., 1., 1.]))
    sv_b = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping_b = [1, 4, 3, 5]
    dof_manager_b = PatchDOFManager(2, mapping_b, global_dof_manager)
    patch_b = Patch(BSplineSurface(su_b, sv_b), mgr, mapping_b, [2, 2], dof_manager_b)

    assembly = PatchAssembly()
    assembly.add_patch(patch_a)
    assembly.add_patch(patch_b)

    materials = [MaterialProperties(210000, 0.3, rho=RHO), MaterialProperties(210000, 0.3)]
    try:
        PatchIntegrator.assemble_mass(assembly, materials)
        assert False, "expected an exception when one material has no rho set"
    except Exception:
        pass


if __name__ == '__main__':
    test_assemble_mass_domain_split_matches_single_patch()
    test_assemble_mass_requires_rho_on_every_material()
