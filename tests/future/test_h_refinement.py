"""
Tests and examples for HRefiner — B-spline knot insertion (h-refinement).

Run directly:   python tests/future/test_h_refinement.py
Via pytest:     pytest tests/future/test_h_refinement.py
"""
import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager, Patch, HRefiner,
    GlobalDOFManager, PatchDOFManager, PatchAssembly, SubdivisionRefiner,
    nd_transition_from_1d,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_surface_patch(p_u, p_v, kv_u, kv_v, cp_coords):
    """
    Build a geometry-only B-spline surface patch (no DOF manager).

    Parameters
    ----------
    p_u, p_v : int
        Polynomial degrees in u and v.
    kv_u, kv_v : array_like
        Knot vectors.
    cp_coords : list of [x, y], length nu*nv
        Control points in u-fastest order: index = iu + iv*nu
        (direction 0 / u varies fastest, direction 1 / v varies slowest).

    Returns
    -------
    patch : Patch
    nu, nv : int
        Number of CPs in each direction.
    """
    kv_u = np.asarray(kv_u, dtype=float)
    kv_v = np.asarray(kv_v, dtype=float)
    nu = len(kv_u) - p_u - 1
    nv = len(kv_v) - p_v - 1
    assert len(cp_coords) == nu * nv, (
        f"Expected {nu * nv} CPs, got {len(cp_coords)}"
    )

    cp_manager = ControlPointManager(dim=2)
    for xy in cp_coords:
        cp_manager.add_point(list(xy))

    su = BSpline(p_u, kv_u)
    sv = BSpline(p_v, kv_v)
    surf = BSplineSurface(su, sv)
    patch = Patch(surf, cp_manager, list(range(nu * nv)), [nu, nv])
    return patch, nu, nv


def get_cp_coords(patch, n_cp):
    """Return (n_cp, 2) array of a patch's control point coordinates."""
    return np.array([patch.control_point(i) for i in range(n_cp)])


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_hrefiner_transition_matrix_shape():
    """Transition matrix must have shape (nb_new_cp, nb_old_cp)."""
    patch, nu, nv = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        # u-fastest: (u=0,v=0),(u=1,v=0),(u=0,v=1),(u=1,v=1)
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )

    # Insert knot in u: nu 2 → 3, nv unchanged
    T = HRefiner(direction=0, knot=0.5).refine(patch)
    assert T.shape == (3 * nv, nu * nv), (
        f"Expected T shape ({3 * nv}, {nu * nv}), got {T.shape}"
    )
    print(f"[PASS] test_hrefiner_transition_matrix_shape  T.shape={T.shape}")


def test_hrefiner_transition_matrix_properties():
    """Each row of T must sum to 1 and be non-negative (convex combination)."""
    patch, _, _ = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )

    T = HRefiner(direction=0, knot=0.5).refine(patch)

    assert np.all(T >= -1e-14), f"T has negative entries:\n{T}"
    row_sums = T.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-14), (
        f"Row sums not equal to 1: {row_sums}"
    )
    print("[PASS] test_hrefiner_transition_matrix_properties")


def test_hrefiner_geometry_preservation_d1():
    """
    Degree-1 bilinear patch, knot insertion at u=0.5.

    New CPs in u-fastest order after inserting u=0.5 (nu: 2→3, nv=2):
      iv=0: (0,0),(0.5,0),(1,0)
      iv=1: (0,1),(0.5,1),(1,1)  ← inserted at iu=1
    """
    patch, _, _ = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )

    old_coords = get_cp_coords(patch, 4)
    T = HRefiner(direction=0, knot=0.5).refine(patch)
    new_coords = get_cp_coords(patch, 6)

    # Geometry preservation: new_cp == T @ old_cp
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"new_cp != T @ old_cp\n  T@old:\n{T @ old_coords}\n  actual:\n{new_coords}"
    )

    # u-fastest expected: iv=0,1; for each iv: iu=0,1,2
    expected = np.array([
        [0., 0.], [0.5, 0.], [1., 0.],   # iv=0
        [0., 1.], [0.5, 1.], [1., 1.]    # iv=1 (inserted at iu=1)
    ])
    assert np.allclose(new_coords, expected, atol=1e-12), (
        f"Unexpected CP positions\n  expected:\n{expected}\n  got:\n{new_coords}"
    )
    print("[PASS] test_hrefiner_geometry_preservation_d1")
    print(f"  old CPs:\n{old_coords}")
    print(f"  new CPs:\n{new_coords}")
    print(f"  T:\n{T}")


def test_hrefiner_direction_v():
    """Insert knot in the v direction: nv 2 → 3, nu unchanged.

    New CPs in u-fastest order (nu=2 unchanged, nv: 2→3):
      iv=0: (0,0),(1,0)
      iv=1: (0,0.5),(1,0.5)  ← inserted
      iv=2: (0,1),(1,1)
    """
    patch, nu, nv = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )

    old_coords = get_cp_coords(patch, 4)
    T = HRefiner(direction=1, knot=0.5).refine(patch)
    new_coords = get_cp_coords(patch, 6)

    assert T.shape == (nu * 3, nu * nv), (
        f"Expected T shape ({nu * 3}, {nu * nv}), got {T.shape}"
    )
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12)

    # u-fastest expected: for each iv: iu=0,1
    expected = np.array([
        [0., 0.], [1., 0.],      # iv=0
        [0., 0.5], [1., 0.5],   # iv=1 (inserted)
        [0., 1.], [1., 1.]      # iv=2
    ])
    assert np.allclose(new_coords, expected, atol=1e-12), (
        f"Unexpected CP positions\n  expected:\n{expected}\n  got:\n{new_coords}"
    )
    print("[PASS] test_hrefiner_direction_v")
    print(f"  new CPs:\n{new_coords}")


def test_hrefiner_geometry_preservation_d2():
    """
    Degree-2 single-element surface, insert knot at u=0.5.
    Verifies both geometry preservation and transition matrix properties.

    3x3 CP grid in u-fastest order (iv outer, iu inner):
      iv=0: (0,0),(1.5,0),(3,0)
      iv=1: (0,0.5),(1.5,0.5),(3,0.5)
      iv=2: (0,1),(1.5,1),(3,1)
    """
    cp_coords = [
        [0., 0.], [1.5, 0.], [3., 0.],
        [0., 0.5], [1.5, 0.5], [3., 0.5],
        [0., 1.], [1.5, 1.], [3., 1.]
    ]
    patch, _, _ = build_surface_patch(
        2, 2,
        kv_u=[0., 0., 0., 1., 1., 1.],
        kv_v=[0., 0., 0., 1., 1., 1.],
        cp_coords=cp_coords
    )

    old_coords = get_cp_coords(patch, 9)
    T = HRefiner(direction=0, knot=0.5).refine(patch)
    # Insert in u: nu 3 → 4, nv stays 3 → 12 new CPs
    new_coords = get_cp_coords(patch, 12)

    assert T.shape == (12, 9), f"Expected (12, 9), got {T.shape}"
    assert np.all(T >= -1e-14)
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"Geometry not preserved\n  T@old:\n{T @ old_coords}\n  actual:\n{new_coords}"
    )
    print("[PASS] test_hrefiner_geometry_preservation_d2")
    print(f"  T ({T.shape[0]}x{T.shape[1]}):\n{T}")
    print(f"  new CPs:\n{new_coords}")


def test_hrefiner_refine_1d_consistency():
    """
    refine_1d() must yield the same geometry as refine() but return only the
    1D transition matrix.  The full nD matrix reconstructed via
    nd_transition_from_1d() must match the one returned by refine().
    """
    # Reference patch: bilinear, 2x2 CPs
    def make_patch():
        return build_surface_patch(
            1, 1,
            kv_u=[0., 0., 1., 1.],
            kv_v=[0., 0., 1., 1.],
            cp_coords=[[0., 0.], [1., 0.],
                       [0., 1.], [1., 1.]]
        )

    # --- full refine ---
    patch_full, nu, nv = make_patch()
    shape_before = [nu, nv]
    T_nd = HRefiner(direction=0, knot=0.5).refine(patch_full)
    coords_full = get_cp_coords(patch_full, patch_full.n_cp)
    shape_after_full = list(patch_full.local_shape)

    # --- refine_1d ---
    patch_1d, nu2, nv2 = make_patch()
    T_1d = HRefiner(direction=0, knot=0.5).refine_1d(patch_1d)
    coords_1d = get_cp_coords(patch_1d, patch_1d.n_cp)
    shape_after_1d = list(patch_1d.local_shape)

    # 1D matrix must be smaller than nD
    assert T_1d.shape == (nu + 1, nu), (
        f"Expected T_1d shape ({nu+1}, {nu}), got {T_1d.shape}"
    )
    assert T_nd.shape == ((nu + 1) * nv, nu * nv), (
        f"Expected T_nd shape ({(nu+1)*nv}, {nu*nv}), got {T_nd.shape}"
    )

    # Geometry must be identical
    assert np.allclose(coords_full, coords_1d, atol=1e-12), (
        f"Geometry mismatch between refine() and refine_1d()\n"
        f"  refine:    {coords_full}\n  refine_1d: {coords_1d}"
    )
    assert shape_after_full == shape_after_1d

    # Reconstructed nD matrix must match the full one
    T_nd_reconstructed = nd_transition_from_1d(T_1d, 0, shape_before, shape_after_1d)
    assert np.allclose(T_nd, T_nd_reconstructed, atol=1e-12), (
        f"nd_transition_from_1d(T_1d) != T from refine()"
    )
    print("[PASS] test_hrefiner_refine_1d_consistency")


def test_hrefiner_refine_1d_with_propagation():
    """
    HRefiner.refine_1d() must work as refine_1d_fn in
    PatchAssembly.refine_with_propagation(), propagating a targeted knot
    insertion across a shared interface.

    Two unit-square patches glued along their right/left edge (u=1 / u=0).
    Insert knot v=0.5 in patch 0; verify it propagates to patch 1.
    """
    def make_unit_square(x_offset=0.0):
        cp_coords = [
            [x_offset + 0., 0.], [x_offset + 1., 0.],
            [x_offset + 0., 1.], [x_offset + 1., 1.],
        ]
        return build_surface_patch(
            1, 1,
            kv_u=[0., 0., 1., 1.],
            kv_v=[0., 0., 1., 1.],
            cp_coords=cp_coords
        )

    patch0, nu0, nv0 = make_unit_square(0.0)
    patch1, nu1, nv1 = make_unit_square(1.0)

    # Share the right boundary of patch0 with the left boundary of patch1
    # Right boundary of patch0 (direction=0, side=1): local indices {1, 3}  (u=1)
    # Left  boundary of patch1 (direction=0, side=0): local indices {0, 2}  (u=0)
    # Re-use global ids from patch0 in patch1 to establish sharing
    shared_right = patch0.boundary_control_points(direction=0, side=1)   # [1, 3]
    left_of_p1   = patch1.boundary_control_points(direction=0, side=0)   # [0, 2]

    # Remap patch1's left-edge local positions to use patch0's global ids
    new_gi = list(patch1.global_indices)
    cp_mgr = patch1.cp_manager
    for loc1, loc0 in zip(left_of_p1, shared_right):
        gid0 = patch0.global_indices[loc0]
        old_gid1 = new_gi[loc1]
        new_gi[loc1] = gid0
        # Patch0's cp_manager holds the shared coords; link patch1 to it
    # Use a single shared cp_manager for both patches
    shared_mgr = patch0.cp_manager
    # Add patch1's private CPs (right boundary) to shared_mgr
    private_coords = [
        [2., 0.],  # local 1 (u=1, v=0) of patch1
        [2., 1.],  # local 3 (u=1, v=1) of patch1
    ]
    private_locs = patch1.boundary_control_points(direction=0, side=1)  # [1, 3]
    for loc, xy in zip(private_locs, private_coords):
        new_gi[loc] = shared_mgr.add_point(list(xy))

    surf1 = BSplineSurface(
        BSpline(1, np.array([0., 0., 1., 1.])),
        BSpline(1, np.array([0., 0., 1., 1.]))
    )
    patch1 = Patch(surf1, shared_mgr, new_gi, [nu1, nv1])

    assembly = PatchAssembly()
    assembly.add_patch(patch0)
    assembly.add_patch(patch1)
    assembly.detect_shared_control_points()
    assembly.detect_interfaces()

    knot = 0.5

    def refine_1d_fn(patch, direction, protected_global_ids):
        HRefiner(direction, knot).refine_1d(patch, protected_global_ids)

    assembly.refine_with_propagation(
        patch_index=0, direction=1, refine_1d_fn=refine_1d_fn
    )

    patches = assembly.get_patchs()
    # Both patches must now have nv = 3 (one knot inserted in v)
    for i, p in enumerate(patches):
        assert p.local_shape[1] == 3, (
            f"Patch {i}: expected nv=3 after propagation, got {p.local_shape[1]}"
        )
    # The shared boundary CPs (in v) must be identical between the two patches
    right_boundary = patches[0].boundary_control_points(direction=0, side=1)
    left_boundary  = patches[1].boundary_control_points(direction=0, side=0)
    for loc0, loc1 in zip(right_boundary, left_boundary):
        gid0 = patches[0].global_indices[loc0]
        gid1 = patches[1].global_indices[loc1]
        assert gid0 == gid1, (
            f"Shared boundary CP mismatch: patch0 gid={gid0}, patch1 gid={gid1}"
        )
    print("[PASS] test_hrefiner_refine_1d_with_propagation")


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    test_hrefiner_transition_matrix_shape()
    test_hrefiner_transition_matrix_properties()
    test_hrefiner_geometry_preservation_d1()
    test_hrefiner_direction_v()
    test_hrefiner_geometry_preservation_d2()
    test_hrefiner_refine_1d_consistency()
    test_hrefiner_refine_1d_with_propagation()
    print("\nAll HRefiner tests passed!")
