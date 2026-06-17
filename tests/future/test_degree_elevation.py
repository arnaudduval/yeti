"""
Tests for PRefiner — B-spline degree elevation (Piegl & Tiller A5.9).

Run directly:   python tests/future/test_degree_elevation.py
Via pytest:     pytest tests/future/test_degree_elevation.py
"""
import numpy as np

from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager, Patch, PRefiner
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_surface_patch(p_u, p_v, kv_u, kv_v, cp_coords):
    """cp_coords in u-fastest order: index = iu + iv*nu."""
    kv_u = np.asarray(kv_u, dtype=float)
    kv_v = np.asarray(kv_v, dtype=float)
    nu = len(kv_u) - p_u - 1
    nv = len(kv_v) - p_v - 1
    assert len(cp_coords) == nu * nv, f"Expected {nu*nv} CPs, got {len(cp_coords)}"
    cpm = ControlPointManager(dim=2)
    for xy in cp_coords:
        cpm.add_point(list(xy))
    surf = BSplineSurface(BSpline(p_u, kv_u), BSpline(p_v, kv_v))
    return Patch(surf, cpm, list(range(nu * nv)), [nu, nv]), nu, nv


def get_cp_coords(patch, n):
    return np.array([patch.control_point(i) for i in range(n)])


def get_degree(patch, direction):
    sp = patch.tensor.components[direction]
    return sp.degree


def get_kv(patch, direction):
    sp = patch.tensor.components[direction]
    return np.array(sp.knot_vector)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_prefiner_degree_increases():
    """Degree in the elevation direction must increase by exactly 1."""
    patch, _, _ = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        # u-fastest: (u=0,v=0),(u=1,v=0),(u=0,v=1),(u=1,v=1)
        cp_coords=[[0.,0.],[1.,0.],[0.,1.],[1.,1.]]
    )
    assert get_degree(patch, 0) == 1
    PRefiner(direction=0).refine(patch)
    assert get_degree(patch, 0) == 2, f"Expected degree 2, got {get_degree(patch, 0)}"
    assert get_degree(patch, 1) == 1, "Direction 1 should be unchanged"
    print(f"[PASS] test_prefiner_degree_increases  kv_u={get_kv(patch,0)}")


def test_prefiner_single_element_max_continuity():
    """
    Single-element patch: no interior knots → elevation preserves max continuity.
    Degree 1 → 2: kv [0,0,1,1] → [0,0,0,1,1,1], no interior knot.
    Degree 3 → 4: kv [0,0,0,0,1,1,1,1] → [0,0,0,0,0,1,1,1,1,1].
    """
    # Degree 1 → 2
    patch, _, _ = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0.,0.],[1.,0.],[0.,1.],[1.,1.]]
    )
    PRefiner(direction=0).refine(patch)
    kv = get_kv(patch, 0)
    interior = kv[(kv > kv[0]) & (kv < kv[-1])]
    assert len(interior) == 0, f"Single-element: no interior knots expected, got {interior}"
    print(f"[PASS] test_prefiner_single_element_max_continuity (d=1→2)  kv={kv}")

    # Degree 3 → 4 (nu=4, nv=2, u-fastest: iv outer, iu inner)
    patch2, _, _ = build_surface_patch(
        3, 1,
        kv_u=[0., 0., 0., 0., 1., 1., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[i, j] for j in range(2) for i in range(4)]
    )
    PRefiner(direction=0).refine(patch2)
    assert get_degree(patch2, 0) == 4
    kv2 = get_kv(patch2, 0)
    interior2 = kv2[(kv2 > kv2[0]) & (kv2 < kv2[-1])]
    assert len(interior2) == 0, f"Expected no interior knots, got {interior2}"
    print(f"[PASS] test_prefiner_single_element_max_continuity (d=3→4)  kv={kv2}")


def test_prefiner_transition_matrix_properties():
    """T rows sum to 1, all entries non-negative."""
    patch, _, _ = build_surface_patch(
        2, 1,
        kv_u=[0., 0., 0., 1., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        # u-fastest: iv=0,1; iu=0,1,2 for each iv
        cp_coords=[[i, j] for j in range(2) for i in range(3)]
    )
    T = PRefiner(direction=0).refine(patch)
    assert np.all(T >= -1e-14), "T has negative entries"
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14), f"Row sums: {T.sum(axis=1)}"
    print(f"[PASS] test_prefiner_transition_matrix_properties  T.shape={T.shape}")


def test_prefiner_geometry_preservation_linear():
    """
    Degree 1 → 2 on a single-element bilinear surface.
    Geometry must be exactly preserved: new_coords == T @ old_coords.

    Old CPs u-fastest (nu=2,nv=2): (0,0),(2,0),(0,3),(2,3)
    New CPs u-fastest (nu=3,nv=2): iv=0:(0,0),(1,0),(2,0)  iv=1:(0,3),(1,3),(2,3)
    """
    patch, _, _ = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        # u-fastest: (u=0,v=0),(u=1,v=0),(u=0,v=1),(u=1,v=1)
        cp_coords=[[0.,0.],[2.,0.],[0.,3.],[2.,3.]]
    )
    old_coords = get_cp_coords(patch, 4)
    T = PRefiner(direction=0).refine(patch)
    # New nu = 3 (p=1→2 adds 1 CP per element), nv=2 → 6 total
    new_coords = get_cp_coords(patch, 6)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"Geometry not preserved\n  T@old:\n{T @ old_coords}\n  actual:\n{new_coords}"
    )
    # u-fastest: iv=0,1; iu=0,1,2 for each iv
    expected = np.array([
        [0.,0.], [1.,0.], [2.,0.],  # iv=0
        [0.,3.], [1.,3.], [2.,3.],  # iv=1 (midpoint at iu=1)
    ])
    assert np.allclose(new_coords, expected, atol=1e-12), (
        f"Wrong CP positions\n  expected:\n{expected}\n  got:\n{new_coords}"
    )
    print("[PASS] test_prefiner_geometry_preservation_linear")
    print(f"  T:\n{T}")
    print(f"  new CPs:\n{new_coords}")


def test_prefiner_geometry_preservation_multi_element():
    """
    Degree 1 → 2 on a 2-element patch: geometry must be preserved.
    Interior knot gets multiplicity 2 (C^0 preserved).

    u-fastest cp_coords for nu=3, nv=2.
    """
    patch, nu, nv = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 0.5, 1., 1.],
        kv_v=[0., 0., 1., 1.],
        # u-fastest: iv=0,1; iu=0,1,2 for each iv
        cp_coords=[[0.,0.],[0.5,0.],[1.,0.],
                   [0.,1.],[0.5,1.],[1.,1.]]
    )
    old_coords = get_cp_coords(patch, nu * nv)
    T = PRefiner(direction=0).refine(patch)

    # nu 3 → 5 (2 elements × 1 extra CP each, minus 1 shared = 3+2=5), nv=2 → 10 CPs
    new_nu = 5
    new_coords = get_cp_coords(patch, new_nu * nv)

    assert T.shape == (new_nu * nv, nu * nv), f"T.shape={T.shape}"
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"Geometry not preserved\n  T@old:\n{T @ old_coords}\n  actual:\n{new_coords}"
    )

    # New kv should have 0.5 with multiplicity 2
    kv = get_kv(patch, 0)
    n05 = np.sum(np.abs(kv - 0.5) < 1e-14)
    assert n05 == 2, f"Interior knot 0.5 should have mult 2, got {n05}"

    print("[PASS] test_prefiner_geometry_preservation_multi_element")
    print(f"  T ({T.shape}):\n{T}")
    print(f"  new kv: {kv}")
    print(f"  new CPs:\n{new_coords}")


def test_prefiner_direction_v():
    """Elevate in v direction: nv increases, nu unchanged.

    New CPs u-fastest (nu=2 unchanged, nv: 2→3):
      iv=0: (0,0),(1,0)
      iv=1: (0,?),(1,?)  ← inserted
      iv=2: (0,1),(1,1)
    """
    patch, nu, nv = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0.,0.],[1.,0.],[0.,1.],[1.,1.]]
    )
    old_coords = get_cp_coords(patch, nu * nv)
    T = PRefiner(direction=1).refine(patch)

    # nv 2 → 3, nu=2 unchanged → 6 CPs
    new_coords = get_cp_coords(patch, 6)
    assert T.shape == (6, 4), f"T.shape={T.shape}"
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12)
    assert get_degree(patch, 0) == 1, "Direction 0 unchanged"
    assert get_degree(patch, 1) == 2, "Direction 1 elevated"
    print(f"[PASS] test_prefiner_direction_v  new CPs:\n{new_coords}")


def test_prefiner_degree3_single_element():
    """Degree 3 → 4: 4 CPs → 5 CPs, geometry preserved, C^3 continuity (no interior knots).

    u-fastest cp_coords for nu=4, nv=2.
    """
    patch, nu, nv = build_surface_patch(
        3, 1,
        kv_u=[0., 0., 0., 0., 1., 1., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        # u-fastest: iv=0,1; iu=0,1,2,3 for each iv
        cp_coords=[[0.,0.],[1.,0.],[2.,0.],[3.,0.],
                   [0.,1.],[1.,1.],[2.,1.],[3.,1.]]
    )
    old_coords = get_cp_coords(patch, nu * nv)
    T = PRefiner(direction=0).refine(patch)
    # nu 4 → 5, nv=2 → 10 CPs
    new_coords = get_cp_coords(patch, 10)

    assert T.shape == (10, 8), f"T.shape={T.shape}"
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"Geometry not preserved\n  T@old:\n{T@old_coords}\n  actual:\n{new_coords}"
    )
    assert get_degree(patch, 0) == 4
    kv = get_kv(patch, 0)
    assert np.allclose(kv, [0,0,0,0,0,1,1,1,1,1]), f"kv={kv}"
    print(f"[PASS] test_prefiner_degree3_single_element  T.shape={T.shape}")
    print(f"  T:\n{T}")


def test_prefiner_multiple_elevations():
    """
    n_elevations=3 on a single-element patch: degree 1 → 4.
    Equivalent to applying PRefiner(dir, 1) three times.
    Geometry and composed T must match.
    """
    def fresh_patch():
        return build_surface_patch(
            1, 1,
            kv_u=[0., 0., 1., 1.],
            kv_v=[0., 0., 1., 1.],
            cp_coords=[[0.,0.],[1.,0.],[0.,1.],[1.,1.]]
        )[0]

    # Three separate elevations
    p3 = fresh_patch()
    old_coords = get_cp_coords(p3, 4)
    T1 = PRefiner(direction=0, n_elevations=1).refine(p3)
    T2 = PRefiner(direction=0, n_elevations=1).refine(p3)
    T3 = PRefiner(direction=0, n_elevations=1).refine(p3)
    T_manual = T3 @ T2 @ T1

    # Single call with n_elevations=3
    p1 = fresh_patch()
    T_combined = PRefiner(direction=0, n_elevations=3).refine(p1)
    coords_combined = get_cp_coords(p1, T_combined.shape[0])

    assert get_degree(p1, 0) == 4, f"Expected degree 4, got {get_degree(p1, 0)}"
    assert T_combined.shape == T_manual.shape, (
        f"Shape mismatch: {T_combined.shape} vs {T_manual.shape}"
    )
    assert np.allclose(T_combined, T_manual, atol=1e-12), "T mismatch"
    assert np.allclose(T_combined.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(T_combined @ old_coords, coords_combined, atol=1e-12)
    print(f"[PASS] test_prefiner_multiple_elevations  T.shape={T_combined.shape}")
    print(f"  new CPs:\n{coords_combined}")


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    test_prefiner_degree_increases()
    test_prefiner_single_element_max_continuity()
    test_prefiner_transition_matrix_properties()
    test_prefiner_geometry_preservation_linear()
    test_prefiner_geometry_preservation_multi_element()
    test_prefiner_direction_v()
    test_prefiner_degree3_single_element()
    test_prefiner_multiple_elevations()
    print("\nAll PRefiner tests passed!")
