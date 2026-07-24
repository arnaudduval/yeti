"""
Test multipatch stiffness assembly (PatchIntegrator.assemble_stiffness)
"""

import os

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (BSpline, BSplineSurface, ControlPointManager,
    Patch, GlobalDOFManager, PatchDOFManager, PatchAssembly, PatchIntegrator,
    Material, PlaneStress, SubdivisionRefiner)

from test_stiffness import stiffness_matrix_legagy


def test_assemble_domain_split_matches_single_patch():
    """
    Splitting a single, legacy-validated patch into two conforming patches
    sharing the interior boundary must not change the assembled physics:
    PatchIntegrator.assemble_stiffness() on the 2-patch assembly must reproduce the
    same global stiffness matrix as the single-patch legacy reference for
    '2_elts_C0_d2_rect' (already validated in test_stiffness.py).
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    stiff_legacy = stiffness_matrix_legagy(f'{script_dir}/2_elts_C0_d2_rect')

    # Same 15 control points as test_integration_2_elements_C0 (nu=5, nv=3,
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

    material = PlaneStress(Material(E=210000, nu=0.3))
    stiffness_matrix = PatchIntegrator.assemble_stiffness(assembly, [material, material])

    assert stiffness_matrix.shape == stiff_legacy.shape
    np.testing.assert_allclose(
        stiffness_matrix.toarray(), stiff_legacy.toarray(), rtol=1.e-5, atol=1.e-8)


def test_assemble_crossed_reversed_interface():
    """
    Exercise assemble_stiffness() on a non-axis-aligned topology (u of one patch glued
    to v of the other, reversed traversal order -- see
    test_assembly.py::test_update_dof_managers_merges_shared_dofs for the
    same fixture) after a propagated refinement. Checks basic structural
    correctness: symmetry, expected size, and that the shared-edge dofs
    carry a stiffer (summed) contribution than either patch would produce
    alone.
    """
    cp_manager = ControlPointManager(dim=2)
    cp_manager.add_point([0.0, 0.0])   # CP 0
    cp_manager.add_point([1.0, 0.0])   # CP 1 (shared)
    cp_manager.add_point([0.0, 1.0])   # CP 2
    cp_manager.add_point([1.0, 1.0])   # CP 3 (shared)
    cp_manager.add_point([2.0, 0.0])   # CP 4
    cp_manager.add_point([2.0, 1.0])   # CP 5

    dofs_per_control_point = [2] * cp_manager.n_points
    global_dof_manager = GlobalDOFManager(dofs_per_control_point)

    su_a = BSpline(1, np.array([0., 0., 1., 1.]))
    sv_a = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping_a = [0, 1, 2, 3]
    dof_manager_a = PatchDOFManager(2, mapping_a, global_dof_manager)
    patch_a = Patch(BSplineSurface(su_a, sv_a), cp_manager, mapping_a, [2, 2], dof_manager_a)

    su_b = BSpline(1, np.array([0., 0., 1., 1.]))
    sv_b = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping_b = [3, 1, 4, 5]
    dof_manager_b = PatchDOFManager(2, mapping_b, global_dof_manager)
    patch_b = Patch(BSplineSurface(su_b, sv_b), cp_manager, mapping_b, [2, 2], dof_manager_b)

    assembly = PatchAssembly()
    assembly.add_patch(patch_a)
    assembly.add_patch(patch_b)
    assembly.detect_shared_control_points()
    assembly.detect_interfaces()

    def refine_1d_fn(patch, direction, protected_global_ids):
        SubdivisionRefiner(direction=direction, n_levels=1).refine_1d(patch, protected_global_ids)

    assembly.refine_with_propagation(patch_index=0, direction=1, refine_1d_fn=refine_1d_fn)
    assembly.update_dof_managers(global_dof_manager, dofs_per_cp=2)

    material = PlaneStress(Material(E=210000, nu=0.3))
    stiffness_matrix = PatchIntegrator.assemble_stiffness(assembly, [material, material]).toarray()

    n_dofs = global_dof_manager.get_dof_indices(cp_manager.n_points - 1)[-1] + 1
    assert stiffness_matrix.shape == (n_dofs, n_dofs)
    np.testing.assert_allclose(stiffness_matrix, stiffness_matrix.T, rtol=1.e-10, atol=1.e-6)

    # The merged midpoint on the shared edge must carry contributions from
    # BOTH patches: its diagonal stiffness must be strictly greater than
    # what patch_a alone would produce for that same dof.
    shared_map = assembly.get_shared_control_points_map()
    midpoint_cp = next(cp for cp, patches in shared_map.items()
                        if len(patches) == 2 and cp not in (1, 3))
    dof = global_dof_manager.get_dof_indices(midpoint_cp)[0]

    material_list_a_only = [material, PlaneStress(Material(E=0.0, nu=0.3))]
    stiffness_a_only = PatchIntegrator.assemble_stiffness(assembly, material_list_a_only).toarray()
    assert stiffness_matrix[dof, dof] > stiffness_a_only[dof, dof]


if __name__ == '__main__':
    test_assemble_domain_split_matches_single_patch()
    test_assemble_crossed_reversed_interface()
