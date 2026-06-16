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


if __name__ == '__main__':
    test_assembly_creation()