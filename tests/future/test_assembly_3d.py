"""
Tests for 3D (BSplineVolume) multi-patch interface detection and refinement
propagation: PatchAssembly.detect_interfaces()/refine_with_propagation() for
patches sharing a FACE (ndim == 3), generalizing the existing 2D (shared
edge) tests in test_assembly.py / test_h_refinement.py.

A 3D face has TWO varying directions, so two matching faces can correspond
in up to 8 ways (the symmetry group of a rectangle): axis transpose x
independent per-axis reversal -- the real new complexity relative to a 2D
edge's single varying direction (only "same order"/"reversed").
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineVolume, Patch, ControlPointManager, PatchAssembly,
    GlobalDOFManager, PatchDOFManager, SubdivisionRefiner, IGABasis1D,
    Material, Elasticity3D, PatchIntegrator,
)


def _refine_1d_fn(patch, direction, protected_ids):
    SubdivisionRefiner(direction=direction, n_levels=1).refine_1d(patch, protected_ids)


def _build_shared_face_assembly(transpose=False, with_dofs=False):
    """
    Two unit-cube BSplineVolume patches (degree 1, 1 element/direction)
    sharing the x=1 face: left x in [0,1], right x in [1,2]. If
    transpose=True, right's own (v0, v1) axes at that face are deliberately
    swapped relative to left's (its v0, the fast/direction-1 axis, holds
    physical z instead of y) to exercise the 8-way matching logic instead
    of only the trivial identity case.
    """
    mgr = ControlPointManager(dim=3)
    coords = {}

    def add(x, y, z):
        key = (x, y, z)
        if key not in coords:
            coords[key] = mgr.add_point([x, y, z])
        return coords[key]

    left_ids = []
    for iw in range(2):
        for iv in range(2):
            for iu in range(2):
                left_ids.append(add(float(iu), float(iv), float(iw)))

    right_ids = []
    if not transpose:
        for iw in range(2):
            for iv in range(2):
                for iu in range(2):
                    right_ids.append(add(float(iu + 1), float(iv), float(iw)))
    else:
        for a in range(2):      # v1_b (slow) -> physical y
            for b in range(2):  # v0_b (fast) -> physical z
                for iu in range(2):
                    right_ids.append(add(float(iu + 1), float(a), float(b)))

    kv = np.array([0., 0., 1., 1.])

    dof_manager = GlobalDOFManager([3] * mgr.n_points) if with_dofs else None

    def make_patch(ids):
        tensor = BSplineVolume(BSpline(1, kv), BSpline(1, kv), BSpline(1, kv))
        if with_dofs:
            pdm = PatchDOFManager(3, ids, dof_manager)
            return Patch(tensor, mgr, ids, [2, 2, 2], pdm)
        return Patch(tensor, mgr, ids, [2, 2, 2])

    patch_left = make_patch(left_ids)
    patch_right = make_patch(right_ids)

    assembly = PatchAssembly()
    assembly.add_patch(patch_left)
    assembly.add_patch(patch_right)
    assembly.detect_shared_control_points()
    assembly.detect_interfaces()
    if with_dofs:
        return assembly, patch_left, patch_right, dof_manager
    return assembly, patch_left, patch_right


def test_detect_interfaces_3d_identity():
    assembly, _, _ = _build_shared_face_assembly(transpose=False)
    itfs = assembly.get_interfaces()
    assert len(itfs) == 1
    itf = itfs[0]
    assert itf['ndim'] == 3
    assert itf['varying_directions_a'] == (1, 2)
    assert itf['varying_directions_b'] == (1, 2)
    assert itf['transpose'] is False
    assert itf['reverse_0'] is False
    assert itf['reverse_1'] is False


def test_detect_interfaces_3d_transpose():
    assembly, _, _ = _build_shared_face_assembly(transpose=True)
    itfs = assembly.get_interfaces()
    assert len(itfs) == 1
    assert itfs[0]['transpose'] is True


def test_refine_with_propagation_3d_both_axes_identity():
    assembly, patch_left, patch_right = _build_shared_face_assembly(transpose=False)

    # A 3D face has two varying directions (1 and 2 here) -- propagate each
    # separately, the stored orientation is reused unchanged by the second call.
    assembly.refine_with_propagation(patch_index=0, direction=1, refine_1d_fn=_refine_1d_fn)
    assembly.refine_with_propagation(patch_index=0, direction=2, refine_1d_fn=_refine_1d_fn)

    assert patch_left.local_shape == [2, 3, 3]
    assert patch_right.local_shape == [2, 3, 3]

    face_left = patch_left.boundary_control_points(0, 1)
    face_right = patch_right.boundary_control_points(0, 0)
    ids_left = [patch_left.global_indices[i] for i in face_left]
    ids_right = [patch_right.global_indices[i] for i in face_right]
    assert ids_left == ids_right  # identity case: same traversal order too

    coords_left = [patch_left.control_point(i) for i in face_left]
    coords_right = [patch_right.control_point(i) for i in face_right]
    for a, b in zip(coords_left, coords_right):
        np.testing.assert_allclose(a, b)


def test_refine_with_propagation_3d_both_axes_transpose():
    assembly, patch_left, patch_right = _build_shared_face_assembly(transpose=True)

    assembly.refine_with_propagation(patch_index=0, direction=1, refine_1d_fn=_refine_1d_fn)
    assembly.refine_with_propagation(patch_index=0, direction=2, refine_1d_fn=_refine_1d_fn)

    face_left = patch_left.boundary_control_points(0, 1)
    face_right = patch_right.boundary_control_points(0, 0)
    ids_left = {patch_left.global_indices[i] for i in face_left}
    ids_right = {patch_right.global_indices[i] for i in face_right}
    # Transposed case: each patch keeps its OWN local traversal convention
    # after the merge, so the two id lists aren't expected to match
    # element-wise -- only as sets (the shared-id merge is what matters).
    assert ids_left == ids_right
    assert len(face_left) == len(face_right) == 9


def test_detect_interfaces_3d_incompatible_ordering_throws():
    """
    Shared ids match as a SET but their arrangement is not any of the 8
    rectangle symmetries (mirrors the 2D "incompatible ordering" check,
    generalized from 2 possible orderings to 8).
    """
    mgr = ControlPointManager(dim=3)
    left_ids = [mgr.add_point([float(iu), float(iv), float(iw)])
                for iw in range(2) for iv in range(2) for iu in range(2)]

    private = [mgr.add_point([2.0, y, z]) for y, z in [(0, 0), (1, 0), (0, 1), (1, 1)]]
    right_ids = [0] * 8
    right_ids[0], right_ids[1] = 1, private[0]
    right_ids[2], right_ids[3] = 7, private[1]   # scrambled (left has 3 here)
    right_ids[4], right_ids[5] = 5, private[2]
    right_ids[6], right_ids[7] = 3, private[3]   # scrambled (left has 7 here)

    kv = np.array([0., 0., 1., 1.])

    def make_patch(ids):
        return Patch(BSplineVolume(BSpline(1, kv), BSpline(1, kv), BSpline(1, kv)),
                     mgr, ids, [2, 2, 2])

    patch_left = make_patch(left_ids)
    patch_right = make_patch(right_ids)

    assembly = PatchAssembly()
    assembly.add_patch(patch_left)
    assembly.add_patch(patch_right)
    assembly.detect_shared_control_points()
    try:
        assembly.detect_interfaces()
        assert False, "expected a RuntimeError for incompatible 3D boundary ordering"
    except RuntimeError:
        pass


def test_merged_3d_interface_dof_managers_and_stiffness():
    """
    End to end: detect + propagate-refine a 3D face interface on both of
    its varying directions, update DOF managers so the merged boundary CPs
    get shared dofs, then assemble the stiffness matrix of each patch and
    check the globally-scattered matrix is symmetric -- a topological/
    assembly-integrity patch test (not a new physics check; the
    single-patch 3D stiffness kernel itself is already validated in
    test_stiffness_3d.py).
    """
    assembly, patch_left, patch_right, global_dof_manager = _build_shared_face_assembly(
        transpose=False, with_dofs=True)

    assembly.refine_with_propagation(patch_index=0, direction=1, refine_1d_fn=_refine_1d_fn)
    assembly.refine_with_propagation(patch_index=0, direction=2, refine_1d_fn=_refine_1d_fn)
    assembly.update_dof_managers(global_dof_manager, dofs_per_cp=3)

    law = Elasticity3D(Material(E=210000.0, nu=0.3))
    n_dof = global_dof_manager.n_control_points() * 3
    K = np.zeros((n_dof, n_dof))

    for patch in (patch_left, patch_right):
        su, sv, sw = patch.tensor.components
        basis_u = IGABasis1D.build(su, su.degree + 1)
        basis_v = IGABasis1D.build(sv, sv.degree + 1)
        basis_w = IGABasis1D.build(sw, sw.degree + 1)
        integrator = PatchIntegrator(patch, basis_u, basis_v, basis_w, law)
        K_local = integrator.integrate_stiffness().toarray()

        for a in range(patch.n_cp):
            dofs_a = patch.dof_manager.get_global_dof_indices(a)
            for b in range(patch.n_cp):
                dofs_b = patch.dof_manager.get_global_dof_indices(b)
                for ia, ga in enumerate(dofs_a):
                    for ib, gb in enumerate(dofs_b):
                        K[ga, gb] += K_local[3*a + ia, 3*b + ib]

    np.testing.assert_allclose(K, K.T, atol=1e-6)
    # Merged interface CPs must carry the SAME global dofs on both patches.
    shared_map = assembly.get_shared_control_points_map()
    assert len(shared_map) == 9  # the fully-merged 3x3 shared face
    for cp_id, patch_indices in shared_map.items():
        expected = global_dof_manager.get_dof_indices(cp_id)
        for pidx in patch_indices:
            patch = patch_left if pidx == 0 else patch_right
            local_pos = patch.global_indices.index(cp_id)
            assert patch.dof_manager.get_global_dof_indices(local_pos) == expected


if __name__ == '__main__':
    test_detect_interfaces_3d_identity()
    test_detect_interfaces_3d_transpose()
    test_refine_with_propagation_3d_both_axes_identity()
    test_refine_with_propagation_3d_both_axes_transpose()
    test_detect_interfaces_3d_incompatible_ordering_throws()
    test_merged_3d_interface_dof_managers_and_stiffness()
