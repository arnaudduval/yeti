"""
Test patch assembly
"""
import numpy as np

#pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (Patch, ControlPointManager,
    GlobalDOFManager, BSpline, BSplineSurface, PatchDOFManager,
    PatchAssembly)

def test_assembly_creation():
    """
    Test assembly of patchs
    """

    cp_manager = ControlPointManager(dim=2)

    # 6 control points for 2 patchs. 2 pointrs are shared
    cp_manager.add_point([0.0, 0.0])   # CP 0
    cp_manager.add_point([1.0, 0.0])   # CP 1 (shared)
    cp_manager.add_point([2.0, 0.0])   # CP 2
    cp_manager.add_point([0.0, 1.0])   # CP 3
    cp_manager.add_point([1.0, 1.0])   # CP 4 (shared)
    cp_manager.add_point([2.0, 1.0])   # CP 5

    dofs_per_control_point = [2] * cp_manager.n_points
    global_dof_manager = GlobalDOFManager(dofs_per_control_point)

    su1 = BSpline(1, np.array([0., 0., 1., 1.]))
    sv1 = BSpline(1, np.array([0., 0., 1., 1.]))
    surf1 = BSplineSurface(su1, sv1)
    # u-fastest mapping (nu=2, nv=2): (iu=0,iv=0)→0,(iu=1,iv=0)→1,(iu=0,iv=1)→3,(iu=1,iv=1)→4
    mapping1 = [0, 1, 3, 4]

    su2 = BSpline(1, np.array([0., 0., 1., 1.]))
    sv2 = BSpline(1, np.array([0., 0., 1., 1.]))
    surf2 = BSplineSurface(su2, sv2)
    # u-fastest mapping: (iu=0,iv=0)→1,(iu=1,iv=0)→2,(iu=0,iv=1)→4,(iu=1,iv=1)→5
    mapping2 = [1, 2, 4, 5]

    local_shape = [2, 2]

    patch_dof_manager1 = PatchDOFManager(
        dofs_per_control_point=2,
        control_points=mapping1,
        global_dof_manager=global_dof_manager
    )
    patch_dof_manager2 = PatchDOFManager(
        dofs_per_control_point=2,
        control_points=mapping2,
        global_dof_manager=global_dof_manager
    )

    # Create patchs
    patch1 = Patch(
        tensor=surf1,
        cp_manager=cp_manager,
        global_indices=mapping1,
        local_shape=local_shape,
        dof_manager=patch_dof_manager1
    )
    patch2 = Patch(
        tensor=surf2,
        cp_manager=cp_manager,
        global_indices=mapping2,
        local_shape=local_shape,
        dof_manager=patch_dof_manager2
    )

    # Create assembly
    assembly = PatchAssembly()
    assembly.add_patch(patch1)
    assembly.add_patch(patch2)

    # Detect shared CP
    assembly.detect_shared_control_points()

    shared_map = assembly.get_shared_control_points_map()
    print('Shared control points:')
    for cp_idx, patch_indices in shared_map.items():
        print(f"Global ID {cp_idx} shared by patchs {patch_indices}")

    # Detect compatible interfaces (Phase 1): patch1's u=last edge ([1, 4])
    # matches patch2's u=first edge ([1, 4]), same direction, same order.
    assembly.detect_interfaces()
    interfaces = assembly.get_interfaces()
    assert len(interfaces) == 1
    itf = interfaces[0]
    assert itf['patch_a'] == 0 and itf['direction_a'] == 0 and itf['side_a'] == 1
    assert itf['patch_b'] == 1 and itf['direction_b'] == 0 and itf['side_b'] == 0
    assert itf['varying_direction_a'] == 1 and itf['varying_direction_b'] == 1
    assert itf['reversed'] is False


def test_detect_interfaces_cross_direction_reversed():
    """
    Interface detection must handle a shared edge that runs along direction
    u for one patch and direction v for the other (e.g. patch A's u=last
    edge glued to patch B's v=first edge), with reversed traversal order.
    """
    cp_manager = ControlPointManager(dim=2)
    cp_manager.add_point([0.0, 0.0])   # CP 0
    cp_manager.add_point([1.0, 0.0])   # CP 1 (shared)
    cp_manager.add_point([0.0, 1.0])   # CP 2
    cp_manager.add_point([1.0, 1.0])   # CP 3 (shared)
    cp_manager.add_point([2.0, 0.0])   # CP 4
    cp_manager.add_point([2.0, 1.0])   # CP 5

    su_a = BSpline(1, np.array([0., 0., 1., 1.]))
    sv_a = BSpline(1, np.array([0., 0., 1., 1.]))
    surf_a = BSplineSurface(su_a, sv_a)
    # u-fastest mapping (nu=2,nv=2): patch A's u=last edge (varying v) -> [1, 3]
    mapping_a = [0, 1, 2, 3]

    su_b = BSpline(1, np.array([0., 0., 1., 1.]))
    sv_b = BSpline(1, np.array([0., 0., 1., 1.]))
    surf_b = BSplineSurface(su_b, sv_b)
    # patch B's v=first edge (varying u) -> [mapping_b[0], mapping_b[1]] = [3, 1]
    # i.e. patch A's [1, 3] in REVERSED order, and along direction v not u.
    mapping_b = [3, 1, 4, 5]

    local_shape = [2, 2]

    patch_a = Patch(surf_a, cp_manager, mapping_a, local_shape)
    patch_b = Patch(surf_b, cp_manager, mapping_b, local_shape)

    assembly = PatchAssembly()
    assembly.add_patch(patch_a)
    assembly.add_patch(patch_b)
    assembly.detect_interfaces()

    interfaces = assembly.get_interfaces()
    assert len(interfaces) == 1
    itf = interfaces[0]
    assert itf['patch_a'] == 0 and itf['direction_a'] == 0 and itf['side_a'] == 1
    assert itf['patch_b'] == 1 and itf['direction_b'] == 1 and itf['side_b'] == 0
    assert itf['varying_direction_a'] == 1   # patch A's edge runs along v
    assert itf['varying_direction_b'] == 0   # patch B's edge runs along u
    assert itf['reversed'] is True


def test_refine_with_propagation_cross_direction_reversed():
    """
    refine_with_propagation must correctly propagate refinement across a
    crossed AND reversed interface (the most demanding combination): the
    neighbor must be refined along its OWN matching direction (not the one
    requested for the first patch), and the independently-created boundary
    midpoints must be merged into a single shared id despite the reversed
    traversal order.
    """
    cp_manager = ControlPointManager(dim=2)
    cp_manager.add_point([0.0, 0.0])   # CP 0
    cp_manager.add_point([1.0, 0.0])   # CP 1 (shared)
    cp_manager.add_point([0.0, 1.0])   # CP 2
    cp_manager.add_point([1.0, 1.0])   # CP 3 (shared)
    cp_manager.add_point([2.0, 0.0])   # CP 4
    cp_manager.add_point([2.0, 1.0])   # CP 5

    su_a = BSpline(1, np.array([0., 0., 1., 1.]))
    sv_a = BSpline(1, np.array([0., 0., 1., 1.]))
    patch_a = Patch(BSplineSurface(su_a, sv_a), cp_manager, [0, 1, 2, 3], [2, 2])

    su_b = BSpline(1, np.array([0., 0., 1., 1.]))
    sv_b = BSpline(1, np.array([0., 0., 1., 1.]))
    patch_b = Patch(BSplineSurface(su_b, sv_b), cp_manager, [3, 1, 4, 5], [2, 2])

    assembly = PatchAssembly()
    assembly.add_patch(patch_a)
    assembly.add_patch(patch_b)
    assembly.detect_shared_control_points()
    assembly.detect_interfaces()

    def refine_1d_fn(patch, direction, protected_global_ids):
        # pylint: disable=import-outside-toplevel
        from yeti_iga.future.bspline import SubdivisionRefiner
        SubdivisionRefiner(direction=direction, n_levels=1).refine_1d(patch, protected_global_ids)

    assembly.refine_with_propagation(patch_index=0, direction=1, refine_1d_fn=refine_1d_fn)

    assert patch_a.n_cp == 6
    assert patch_b.n_cp == 6

    shared_map = assembly.get_shared_control_points_map()
    assert len(shared_map) == 3   # 2 original corners + 1 merged midpoint

    # patch_a's edge (u=last, varying v): x=1, y = 0, 0.5, 1
    edge_a = [patch_a.control_point(i) for i in (1, 3, 5)]
    np.testing.assert_allclose(edge_a, [[1, 0], [1, 0.5], [1, 1]])

    # patch_b's edge (v=first, varying u): x=1, y = 1, 0.5, 0 (reversed order)
    edge_b = [patch_b.control_point(i) for i in (0, 1, 2)]
    np.testing.assert_allclose(edge_b, [[1, 1], [1, 0.5], [1, 0]])


def test_update_dof_managers_merges_shared_dofs():
    """
    update_dof_managers must rebuild every patch's PatchDOFManager from its
    current global_indices, so that a boundary control point merged by
    refine_with_propagation (across a crossed AND reversed interface) ends
    up with the SAME global dofs on both patches -- without any explicit
    tracking of which control point ids were merged.
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
        # pylint: disable=import-outside-toplevel
        from yeti_iga.future.bspline import SubdivisionRefiner
        SubdivisionRefiner(direction=direction, n_levels=1).refine_1d(patch, protected_global_ids)

    assembly.refine_with_propagation(patch_index=0, direction=1, refine_1d_fn=refine_1d_fn)

    assembly.update_dof_managers(global_dof_manager, dofs_per_cp=2)

    shared_map = assembly.get_shared_control_points_map()
    assert len(shared_map) == 3   # 2 original corners + 1 merged midpoint

    for cp_id, patch_indices in shared_map.items():
        expected = global_dof_manager.get_dof_indices(cp_id)
        for pidx in patch_indices:
            patch = patch_a if pidx == 0 else patch_b
            local_pos = patch.global_indices.index(cp_id)
            got = patch.dof_manager.get_global_dof_indices(local_pos)
            assert got == expected


if __name__ == '__main__':
    test_assembly_creation()
    test_detect_interfaces_cross_direction_reversed()
    test_refine_with_propagation_cross_direction_reversed()
    test_update_dof_managers_merges_shared_dofs()