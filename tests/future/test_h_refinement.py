"""
Tests and examples for HRefiner — B-spline knot insertion (h-refinement).

Run directly:   python tests/future/test_h_refinement.py
Via pytest:     pytest tests/future/test_h_refinement.py
"""
import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager, Patch, HRefiner
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
        Control points, row-major order: index = iv * nu + iu.

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
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )

    # Insert knot in u: nu 2 → 3, nv unchanged
    _, T = HRefiner(direction=0, knot=0.5).refine(patch)
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

    _, T = HRefiner(direction=0, knot=0.5).refine(patch)

    assert np.all(T >= -1e-14), f"T has negative entries:\n{T}"
    row_sums = T.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-14), (
        f"Row sums not equal to 1: {row_sums}"
    )
    print("[PASS] test_hrefiner_transition_matrix_properties")


def test_hrefiner_geometry_preservation_d1():
    """
    Degree-1 bilinear patch, knot insertion at u=0.5.

    Expected result: inserting 0.5 in [0,0,1,1] with degree 1 gives
    new CPs that interpolate midpoints between old ones.
    """
    patch, _, _ = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )

    old_coords = get_cp_coords(patch, 4)
    new_patch, T = HRefiner(direction=0, knot=0.5).refine(patch)
    new_coords = get_cp_coords(new_patch, 6)

    # Geometry preservation: new_cp == T @ old_cp
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"new_cp != T @ old_cp\n  T@old:\n{T @ old_coords}\n  actual:\n{new_coords}"
    )

    # Analytical check
    expected = np.array([
        [0., 0.], [0.5, 0.], [1., 0.],  # v=0 row
        [0., 1.], [0.5, 1.], [1., 1.]   # v=1 row
    ])
    assert np.allclose(new_coords, expected, atol=1e-12), (
        f"Unexpected CP positions\n  expected:\n{expected}\n  got:\n{new_coords}"
    )
    print("[PASS] test_hrefiner_geometry_preservation_d1")
    print(f"  old CPs:\n{old_coords}")
    print(f"  new CPs:\n{new_coords}")
    print(f"  T:\n{T}")


def test_hrefiner_direction_v():
    """Insert knot in the v direction: nv 2 → 3, nu unchanged."""
    patch, nu, nv = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )

    old_coords = get_cp_coords(patch, 4)
    new_patch, T = HRefiner(direction=1, knot=0.5).refine(patch)
    new_coords = get_cp_coords(new_patch, 6)

    assert T.shape == (nu * 3, nu * nv), (
        f"Expected T shape ({nu * 3}, {nu * nv}), got {T.shape}"
    )
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12)

    expected = np.array([
        [0., 0.], [1., 0.],    # v=0 row
        [0., 0.5], [1., 0.5],  # v=0.5 row (inserted)
        [0., 1.], [1., 1.]     # v=1 row
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
    """
    # 3x3 CPs on a [0,3]x[0,1] rectangle
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
    new_patch, T = HRefiner(direction=0, knot=0.5).refine(patch)
    # Insert in u: nu 3 → 4, nv stays 3 → 12 new CPs
    new_coords = get_cp_coords(new_patch, 12)

    assert T.shape == (12, 9), f"Expected (12, 9), got {T.shape}"
    assert np.all(T >= -1e-14)
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"Geometry not preserved\n  T@old:\n{T @ old_coords}\n  actual:\n{new_coords}"
    )
    print("[PASS] test_hrefiner_geometry_preservation_d2")
    print(f"  T ({T.shape[0]}x{T.shape[1]}):\n{T}")
    print(f"  new CPs:\n{new_coords}")


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    test_hrefiner_transition_matrix_shape()
    test_hrefiner_transition_matrix_properties()
    test_hrefiner_geometry_preservation_d1()
    test_hrefiner_direction_v()
    test_hrefiner_geometry_preservation_d2()
    print("\nAll HRefiner tests passed!")
