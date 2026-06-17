"""
Tests for BezierExtractor — Bézier extraction operator (Borden et al. 2011).

For a B-spline of degree p with ne elements, BezierExtractor returns ne matrices
C^e (shape p+1 × p+1) such that:

    N_active_e(xi) = C^e @ B^p(xi_hat)

where N_active_e are the active B-spline basis functions on element e and B^p are
the Bernstein polynomials on [0,1].

Run directly:   python tests/future/test_bezier_extraction.py
Via pytest:     pytest tests/future/test_bezier_extraction.py
"""
import numpy as np
import pytest

from yeti_iga.future.bspline import BSpline, BSplineSurface, ControlPointManager, Patch, BezierExtractor


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def bernstein(p, j, xi_hat):
    """Bernstein polynomial B^p_j evaluated at xi_hat ∈ [0,1]."""
    from math import comb
    return comb(p, j) * (xi_hat ** j) * ((1 - xi_hat) ** (p - j))


def bernstein_vec(p, xi_hat):
    """Full Bernstein vector [B^p_0(xi_hat), ..., B^p_p(xi_hat)]."""
    return np.array([bernstein(p, j, xi_hat) for j in range(p + 1)])


def eval_bspline_all(spline, xi):
    """Evaluate all B-spline basis functions at xi using one_basis_fun."""
    n = len(spline.knot_vector) - spline.degree - 1
    return np.array([spline.one_basis_fun(xi, i) for i in range(n)])


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_single_element_identity():
    """Single-element patch: extraction operator must be the identity."""
    for p in [1, 2, 3]:
        kv = np.array([0.] * (p+1) + [1.] * (p+1))
        spline = BSpline(p, kv)
        Ce = BezierExtractor.extract_1d(spline)
        assert len(Ce) == 1, f"p={p}: expected 1 element, got {len(Ce)}"
        assert np.allclose(Ce[0], np.eye(p+1), atol=1e-14), f"p={p}: expected identity, got\n{Ce[0]}"
    print("[PASS] test_single_element_identity")


def test_element_count():
    """Number of extraction operators matches number of non-zero knot spans."""
    # degree 2, 3 elements
    kv = np.array([0., 0., 0., 1./3, 2./3, 1., 1., 1.])
    spline = BSpline(2, kv)
    Ce = BezierExtractor.extract_1d(spline)
    assert len(Ce) == 3, f"Expected 3 elements, got {len(Ce)}"
    print(f"[PASS] test_element_count  ne={len(Ce)}")


def test_column_sum_unity():
    """
    Column sums of C^e must equal 1.

    sum_i N_i(xi) = 1  (B-spline partition of unity)
    = sum_i (sum_j C^e[i,j] B^p_j) = sum_j (sum_i C^e[i,j]) B^p_j
    Since sum_j B^p_j = 1, we need: sum_i C^e[i,j] = 1 for every j.
    i.e., C.sum(axis=0) == [1, 1, ..., 1].
    """
    kv = np.array([0., 0., 0., 0.5, 1., 1., 1.])
    spline = BSpline(2, kv)
    Ce = BezierExtractor.extract_1d(spline)
    for e, C in enumerate(Ce):
        col_sums = C.sum(axis=0)
        assert np.allclose(col_sums, 1.0, atol=1e-14), (
            f"Element {e}: column sums = {col_sums}, expected all 1"
        )
    print("[PASS] test_column_sum_unity")


def test_nonnegative():
    """C^e must have non-negative entries (convex combination)."""
    kv = np.array([0., 0., 0., 0.5, 1., 1., 1.])
    spline = BSpline(2, kv)
    Ce = BezierExtractor.extract_1d(spline)
    for e, C in enumerate(Ce):
        assert np.all(C >= -1e-14), f"Element {e}: negative entry\n{C}"
    print("[PASS] test_nonnegative")


def test_reconstruction_p2_two_elements():
    """
    Verify N_active_e(xi) = C^e @ B^p(xi_hat) pointwise.

    kv = [0,0,0, 0.5, 1,1,1], p=2, 4 basis functions, 2 elements.
    """
    kv = np.array([0., 0., 0., 0.5, 1., 1., 1.])
    spline = BSpline(2, kv)
    Ce = BezierExtractor.extract_1d(spline)
    spans = BezierExtractor.element_spans(spline)
    p = spline.degree

    assert len(Ce) == 2

    # Active basis function indices: [span-p, ..., span] for each element
    for e, (C, span) in enumerate(zip(Ce, spans)):
        xi_a = spline.knot_vector[span]
        xi_b = spline.knot_vector[span + 1]
        active = list(range(span - p, span + 1))

        for xi_hat in np.linspace(0.0, 1.0, 11):
            xi = xi_a + xi_hat * (xi_b - xi_a)
            # Direct B-spline evaluation
            N_all = eval_bspline_all(spline, xi)
            N_active = N_all[active]
            # Bézier reconstruction
            B = bernstein_vec(p, xi_hat)
            N_reconstructed = C @ B
            assert np.allclose(N_active, N_reconstructed, atol=1e-12), (
                f"Element {e}, xi={xi:.4f} (xi_hat={xi_hat:.2f}):\n"
                f"  N_active     = {N_active}\n"
                f"  C^e @ B^p    = {N_reconstructed}\n"
                f"  diff         = {N_active - N_reconstructed}"
            )

    print("[PASS] test_reconstruction_p2_two_elements")


def test_reconstruction_p2_three_elements():
    """Same check with 3 elements of unequal size."""
    kv = np.array([0., 0., 0., 0.3, 0.7, 1., 1., 1.])
    spline = BSpline(2, kv)
    Ce = BezierExtractor.extract_1d(spline)
    spans = BezierExtractor.element_spans(spline)
    p = spline.degree
    assert len(Ce) == 3

    for e, (C, span) in enumerate(zip(Ce, spans)):
        xi_a = spline.knot_vector[span]
        xi_b = spline.knot_vector[span + 1]
        active = list(range(span - p, span + 1))

        for xi_hat in np.linspace(0.0, 1.0, 11):
            xi = xi_a + xi_hat * (xi_b - xi_a)
            N_active = eval_bspline_all(spline, xi)[active]
            N_reconstructed = C @ bernstein_vec(p, xi_hat)
            assert np.allclose(N_active, N_reconstructed, atol=1e-12), (
                f"Element {e}, xi={xi:.4f}: mismatch"
            )

    print("[PASS] test_reconstruction_p2_three_elements")


def test_reconstruction_p3():
    """Degree 3, 2 elements."""
    kv = np.array([0., 0., 0., 0., 0.5, 1., 1., 1., 1.])
    spline = BSpline(3, kv)
    Ce = BezierExtractor.extract_1d(spline)
    spans = BezierExtractor.element_spans(spline)
    p = spline.degree
    assert len(Ce) == 2

    for e, (C, span) in enumerate(zip(Ce, spans)):
        xi_a = spline.knot_vector[span]
        xi_b = spline.knot_vector[span + 1]
        active = list(range(span - p, span + 1))

        for xi_hat in np.linspace(0.0, 1.0, 11):
            xi = xi_a + xi_hat * (xi_b - xi_a)
            N_active = eval_bspline_all(spline, xi)[active]
            N_reconstructed = C @ bernstein_vec(p, xi_hat)
            assert np.allclose(N_active, N_reconstructed, atol=1e-12), (
                f"p=3 Element {e}, xi={xi:.4f}: mismatch"
            )

    print("[PASS] test_reconstruction_p3")


def test_reconstruction_full_continuity():
    """C^0 interior knot (multiplicity p-1): more blending needed."""
    # degree 2, interior knot at 0.5 with multiplicity 1 → C^1 (standard)
    kv = np.array([0., 0., 0., 0.5, 0.5, 1., 1., 1.])  # mult=2=p → C^0
    spline = BSpline(2, kv)
    Ce = BezierExtractor.extract_1d(spline)
    spans = BezierExtractor.element_spans(spline)
    p = spline.degree
    assert len(Ce) == 2

    for e, (C, span) in enumerate(zip(Ce, spans)):
        xi_a = spline.knot_vector[span]
        xi_b = spline.knot_vector[span + 1]
        active = list(range(span - p, span + 1))

        for xi_hat in np.linspace(0.0, 1.0, 11):
            xi = xi_a + xi_hat * (xi_b - xi_a)
            N_active = eval_bspline_all(spline, xi)[active]
            N_reconstructed = C @ bernstein_vec(p, xi_hat)
            assert np.allclose(N_active, N_reconstructed, atol=1e-12), (
                f"C^0 Element {e}, xi={xi:.4f}: mismatch"
            )

    print("[PASS] test_reconstruction_full_continuity")


def test_patch_extract():
    """BezierExtractor.extract(patch) delegates to the correct direction."""
    su = BSpline(2, np.array([0., 0., 0., 0.5, 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mgr = ControlPointManager(dim=2)
    nu, nv = 4, 2
    for iv in range(nv):
        for iu in range(nu):
            mgr.add_point([float(iu), float(iv)])
    patch = Patch(surf, mgr, list(range(nu * nv)), [nu, nv])

    Ce_u = BezierExtractor(direction=0).extract(patch)
    Ce_v = BezierExtractor(direction=1).extract(patch)

    assert len(Ce_u) == 2, f"Expected 2 elements in u, got {len(Ce_u)}"
    assert len(Ce_v) == 1, f"Expected 1 element in v, got {len(Ce_v)}"
    assert Ce_u[0].shape == (3, 3)  # p_u+1 = 3
    assert Ce_v[0].shape == (2, 2)  # p_v+1 = 2

    print(f"[PASS] test_patch_extract  |Ce_u|={len(Ce_u)}, |Ce_v|={len(Ce_v)}")


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    test_single_element_identity()
    test_element_count()
    test_row_sum_unity()
    test_nonnegative()
    test_reconstruction_p2_two_elements()
    test_reconstruction_p2_three_elements()
    test_reconstruction_p3()
    test_reconstruction_full_continuity()
    test_patch_extract()
    print("\nAll BezierExtractor tests passed!")
