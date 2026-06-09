"""
Tests for SubdivisionRefiner — bisection of all knot spans (uniform h-refinement).

Run directly:   python tests/future/test_subdivision.py
Via pytest:     pytest tests/future/test_subdivision.py
"""
import numpy as np

from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager, Patch, SubdivisionRefiner
)


# ---------------------------------------------------------------------------
# Helpers (shared with test_h_refinement)
# ---------------------------------------------------------------------------

def build_surface_patch(p_u, p_v, kv_u, kv_v, cp_coords):
    kv_u = np.asarray(kv_u, dtype=float)
    kv_v = np.asarray(kv_v, dtype=float)
    nu = len(kv_u) - p_u - 1
    nv = len(kv_v) - p_v - 1
    assert len(cp_coords) == nu * nv

    cp_manager = ControlPointManager(dim=2)
    for xy in cp_coords:
        cp_manager.add_point(list(xy))

    surf = BSplineSurface(BSpline(p_u, kv_u), BSpline(p_v, kv_v))
    return Patch(surf, cp_manager, list(range(nu * nv)), [nu, nv]), nu, nv


def get_cp_coords(patch, n_cp):
    return np.array([patch.control_point(i) for i in range(n_cp)])


def count_elements(patch, direction):
    """Number of non-zero knot spans in the given direction."""
    kv = patch.tensor.components[direction].knot_vector
    return sum(1 for i in range(len(kv) - 1) if kv[i + 1] - kv[i] > 1e-14)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_subdivision_doubles_elements():
    """One level of subdivision halves each span → element count doubles."""
    patch, nu, nv = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )
    assert count_elements(patch, 0) == 1

    T = SubdivisionRefiner(direction=0, n_levels=1).refine(patch)

    assert count_elements(patch, 0) == 2, (
        f"Expected 2 elements after subdivision, got {count_elements(patch, 0)}"
    )
    print(f"[PASS] test_subdivision_doubles_elements  T.shape={T.shape}")


def test_subdivision_transition_matrix_properties():
    """T rows must sum to 1 and be non-negative."""
    patch, _, _ = build_surface_patch(
        2, 1,
        kv_u=[0., 0., 0., 1., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[i, j] for j in range(2) for i in range(3)]
    )

    T = SubdivisionRefiner(direction=0, n_levels=1).refine(patch)

    assert np.all(T >= -1e-14), f"T has negative entries"
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14), (
        f"Row sums: {T.sum(axis=1)}"
    )
    print(f"[PASS] test_subdivision_transition_matrix_properties  T.shape={T.shape}")


def test_subdivision_geometry_preservation():
    """new_coords == T @ old_coords for a degree-1 bilinear patch."""
    patch, _, _ = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [2., 0.],
                   [0., 1.], [2., 1.]]
    )

    old_coords = get_cp_coords(patch, 4)
    T = SubdivisionRefiner(direction=0, n_levels=1).refine(patch)
    new_coords = get_cp_coords(patch, patch.tensor.components[0].knot_vector.shape[0]
                                      - patch.tensor.components[0].degree - 1
                                      # nu_new
                                      )

    # nu_new * nv = 3*2 = 6
    new_coords = get_cp_coords(patch, 6)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"Geometry not preserved\n  T@old:\n{T @ old_coords}\n  actual:\n{new_coords}"
    )
    expected = np.array([
        [0., 0.], [1., 0.], [2., 0.],
        [0., 1.], [1., 1.], [2., 1.]
    ])
    assert np.allclose(new_coords, expected, atol=1e-12), (
        f"Wrong CP positions\n  expected:\n{expected}\n  got:\n{new_coords}"
    )
    print(f"[PASS] test_subdivision_geometry_preservation")
    print(f"  T:\n{T}")
    print(f"  new CPs:\n{new_coords}")


def test_subdivision_two_levels():
    """Two levels of subdivision → 4 elements from 1."""
    patch, _, _ = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [1., 0.],
                   [0., 1.], [1., 1.]]
    )

    old_coords = get_cp_coords(patch, 4)
    T = SubdivisionRefiner(direction=0, n_levels=2).refine(patch)

    assert count_elements(patch, 0) == 4, (
        f"Expected 4 elements after 2 levels, got {count_elements(patch, 0)}"
    )
    new_n = count_elements(patch, 0) + 1  # degree 1: nb_cp = nb_elems + 1
    new_coords = get_cp_coords(patch, new_n * 2)  # *nv=2
    assert T.shape == (new_n * 2, 4), f"T.shape={T.shape}"
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12)
    print(f"[PASS] test_subdivision_two_levels  T.shape={T.shape}")


def test_subdivision_multi_element_patch():
    """Patch with 2 initial elements: subdivision must give 4 elements."""
    patch, nu, nv = build_surface_patch(
        1, 1,
        kv_u=[0., 0., 0.5, 1., 1.],
        kv_v=[0., 0., 1., 1.],
        cp_coords=[[0., 0.], [0.5, 0.], [1., 0.],
                   [0., 1.], [0.5, 1.], [1., 1.]]
    )
    assert count_elements(patch, 0) == 2

    old_coords = get_cp_coords(patch, nu * nv)
    T = SubdivisionRefiner(direction=0, n_levels=1).refine(patch)

    assert count_elements(patch, 0) == 4, (
        f"Expected 4 elements, got {count_elements(patch, 0)}"
    )
    assert np.allclose(T.sum(axis=1), 1.0, atol=1e-14)
    new_nu = count_elements(patch, 0) + 1  # degree 1
    new_coords = get_cp_coords(patch, new_nu * nv)
    assert np.allclose(new_coords, T @ old_coords, atol=1e-12), (
        f"Geometry not preserved\n  T@old:\n{T @ old_coords}\n  actual:\n{new_coords}"
    )
    print(f"[PASS] test_subdivision_multi_element_patch  T.shape={T.shape}")
    print(f"  new CPs:\n{new_coords}")


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    test_subdivision_doubles_elements()
    test_subdivision_transition_matrix_properties()
    test_subdivision_geometry_preservation()
    test_subdivision_two_levels()
    test_subdivision_multi_element_patch()
    print("\nAll SubdivisionRefiner tests passed!")
