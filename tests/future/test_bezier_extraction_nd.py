"""
Tests for BezierExtractor.extract_nd() — ND tensor-product Bézier extraction.

For each ND element e = (e_0, ..., e_{n-1}) the extraction operator C satisfies:

    N_active(xi) = C @ B_nd(xi_hat)

where:
  - N_active  are the prod_d(p_d+1) active B-spline basis functions on this element
  - B_nd      is the tensor-product Bernstein basis (u-fastest on [0,1]^ndim)
  - C         has shape (n_local, n_local) with n_local = prod_d(p_d+1)

Run directly:   python tests/future/test_bezier_extraction_nd.py
Via pytest:     pytest tests/future/test_bezier_extraction_nd.py
"""
import numpy as np
import pytest

from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager, Patch, BezierExtractor
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def bernstein(p, j, t):
    from math import comb
    return comb(p, j) * (t**j) * ((1 - t)**(p - j))


def bernstein_nd(degrees, xi_hats):
    """Tensor-product Bernstein vector (u-fastest) at xi_hats ∈ [0,1]^ndim."""
    vecs = [np.array([bernstein(p, j, t) for j in range(p + 1)])
            for p, t in zip(degrees, xi_hats)]
    # Kronecker product in u-fastest order: kron(v_{n-1}, ..., kron(v_1, v_0))
    result = vecs[0]
    for v in vecs[1:]:
        result = np.kron(v, result)
    return result


def build_patch(degrees, kvs, local_shape):
    """Build a unit-square/cube patch with identity CP mapping."""
    splines = [BSpline(p, np.asarray(kv, float)) for p, kv in zip(degrees, kvs)]
    if len(splines) == 2:
        from yeti_iga.future.bspline import BSplineSurface
        tensor = BSplineSurface(splines[0], splines[1])
    else:
        from yeti_iga.future.bspline import BSplineVolume
        tensor = BSplineVolume(splines[0], splines[1], splines[2])

    n_cp = 1
    for n in local_shape:
        n_cp *= n

    mgr = ControlPointManager(dim=len(degrees))
    for i in range(n_cp):
        mgr.add_point([0.0] * len(degrees))   # coords don't matter for basis tests

    return Patch(tensor, mgr, list(range(n_cp)), local_shape)


def eval_bspline_nd_active(patch, xi, active_indices):
    """
    Evaluate each active B-spline basis function by evaluating each 1D factor
    using one_basis_fun and taking the tensor product.
    """
    ndim = len(patch.tensor.components)
    n_local = len(active_indices)

    # Convert flat active indices back to multi-indices, then evaluate each
    local_shapes = patch.local_shape
    strides = [1]
    for d in range(ndim - 1):
        strides.append(strides[-1] * local_shapes[d])

    values = []
    for flat in active_indices:
        # Decode u-fastest flat index into per-direction CP index
        rem = flat
        cp_idx = []
        for d in range(ndim):
            cp_idx.append(rem % local_shapes[d])
            rem //= local_shapes[d]

        # Evaluate each 1D basis function at xi[d] using one_basis_fun
        val = 1.0
        for d, spline in enumerate(patch.tensor.components):
            val *= spline.one_basis_fun(xi[d], cp_idx[d])
        values.append(val)

    return np.array(values)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_nd_element_count_2d():
    """2D: total elements = ne_u * ne_v."""
    patch = build_patch(
        degrees=[2, 1],
        kvs=[[0,0,0, 0.5, 1,1,1], [0,0, 1,1]],
        local_shape=[4, 2]
    )
    elems = BezierExtractor.extract_nd(patch)
    assert len(elems) == 2 * 1, f"Expected 2 elements, got {len(elems)}"
    print(f"[PASS] test_nd_element_count_2d  ne={len(elems)}")


def test_nd_element_count_2d_multi():
    """2D: 3 × 2 = 6 elements."""
    patch = build_patch(
        degrees=[2, 1],
        kvs=[[0,0,0, 1/3, 2/3, 1,1,1], [0,0, 0.5, 1,1]],
        local_shape=[5, 3]
    )
    elems = BezierExtractor.extract_nd(patch)
    assert len(elems) == 3 * 2, f"Expected 6 elements, got {len(elems)}"
    print(f"[PASS] test_nd_element_count_2d_multi  ne={len(elems)}")


def test_nd_matrix_shape():
    """C must be square of size prod(p_d + 1)."""
    patch = build_patch(
        degrees=[2, 1],
        kvs=[[0,0,0, 0.5, 1,1,1], [0,0, 1,1]],
        local_shape=[4, 2]
    )
    elems = BezierExtractor.extract_nd(patch)
    n_local = (2 + 1) * (1 + 1)   # (p_u+1)*(p_v+1)
    for elem in elems:
        assert elem.C.shape == (n_local, n_local), (
            f"Expected ({n_local},{n_local}), got {elem.C.shape}"
        )
    print(f"[PASS] test_nd_matrix_shape  n_local={n_local}")


def test_nd_active_indices_count():
    """Each element has exactly prod(p_d+1) active indices."""
    patch = build_patch(
        degrees=[2, 1],
        kvs=[[0,0,0, 0.5, 1,1,1], [0,0, 1,1]],
        local_shape=[4, 2]
    )
    elems = BezierExtractor.extract_nd(patch)
    n_local = (2 + 1) * (1 + 1)
    for elem in elems:
        assert len(elem.active_indices) == n_local, (
            f"Expected {n_local} active indices, got {len(elem.active_indices)}"
        )
    print(f"[PASS] test_nd_active_indices_count")


def test_nd_column_sum_unity():
    """Column sums of C^e must equal 1 (ND partition of unity)."""
    patch = build_patch(
        degrees=[2, 1],
        kvs=[[0,0,0, 0.5, 1,1,1], [0,0, 1,1]],
        local_shape=[4, 2]
    )
    elems = BezierExtractor.extract_nd(patch)
    for e, elem in enumerate(elems):
        col_sums = elem.C.sum(axis=0)
        assert np.allclose(col_sums, 1.0, atol=1e-14), (
            f"Element {elem.elem_index}: column sums = {col_sums}"
        )
    print("[PASS] test_nd_column_sum_unity")


def test_nd_nonnegative():
    """C^e must have non-negative entries."""
    patch = build_patch(
        degrees=[2, 2],
        kvs=[[0,0,0, 0.5, 1,1,1], [0,0,0, 0.5, 1,1,1]],
        local_shape=[4, 4]
    )
    elems = BezierExtractor.extract_nd(patch)
    for elem in elems:
        assert np.all(elem.C >= -1e-14), f"Negative entry in C for elem {elem.elem_index}"
    print("[PASS] test_nd_nonnegative")


def test_nd_elem_index_ufastest():
    """Elements are returned in u-fastest order."""
    patch = build_patch(
        degrees=[2, 1],
        kvs=[[0,0,0, 1/3, 2/3, 1,1,1], [0,0, 0.5, 1,1]],
        local_shape=[5, 3]
    )
    elems = BezierExtractor.extract_nd(patch)
    # u-fastest: e_0 (u) varies fastest
    expected = [(eu, ev) for ev in range(2) for eu in range(3)]
    for elem, (eu, ev) in zip(elems, expected):
        assert tuple(elem.elem_index) == (eu, ev), (
            f"Expected elem_index {(eu,ev)}, got {elem.elem_index}"
        )
    print("[PASS] test_nd_elem_index_ufastest")


def patch_local_shape(patch):
    """Return [n_0, n_1, ...] from the patch's spline components."""
    return [len(s.knot_vector) - s.degree - 1
            for s in patch.tensor.components]


def eval_active_nd(patch, xi, active_indices):
    """
    Evaluate each active B-spline function at xi by directly calling one_basis_fun
    per direction and taking the tensor product.

    This is the correct reference: it uses the element's own active indices,
    so it works at knot boundaries too (unlike basis_funs_nd which can assign
    a boundary point to the next element's span).
    """
    local_shapes = patch_local_shape(patch)
    values = []
    for flat in active_indices:
        rem = flat
        cp_idx = []
        for n in local_shapes:
            cp_idx.append(rem % n)
            rem //= n
        val = 1.0
        for d, spline in enumerate(patch.tensor.components):
            val *= spline.one_basis_fun(xi[d], cp_idx[d])
        values.append(val)
    return np.array(values)


def test_nd_reconstruction_2d():
    """
    N_active(xi) = C^e @ B_nd(xi_hat) for every point on a 2D patch.

    Reference: evaluate each active B-spline function via one_basis_fun (product
    per direction), which is correct even at knot boundaries.
    """
    kv_u = np.array([0., 0., 0., 0.5, 1., 1., 1.])
    kv_v = np.array([0., 0., 1., 1.])
    su = BSpline(2, kv_u)
    sv = BSpline(1, kv_v)
    from yeti_iga.future.bspline import BSplineSurface
    surf = BSplineSurface(su, sv)
    nu, nv = 4, 2
    mgr = ControlPointManager(dim=2)
    for iv in range(nv):
        for iu in range(nu):
            mgr.add_point([float(iu), float(iv)])
    patch = Patch(surf, mgr, list(range(nu * nv)), [nu, nv])

    degrees = [su.degree, sv.degree]
    elems = BezierExtractor.extract_nd(patch)
    spans_u = BezierExtractor.element_spans(su)
    spans_v = BezierExtractor.element_spans(sv)

    for elem in elems:
        eu, ev = elem.elem_index
        xi_u_a, xi_u_b = kv_u[spans_u[eu]], kv_u[spans_u[eu] + 1]
        xi_v_a, xi_v_b = kv_v[spans_v[ev]], kv_v[spans_v[ev] + 1]

        for t_u in np.linspace(0.0, 1.0, 7):
            for t_v in np.linspace(0.0, 1.0, 7):
                xi = np.array([xi_u_a + t_u * (xi_u_b - xi_u_a),
                               xi_v_a + t_v * (xi_v_b - xi_v_a)])

                N_active_ref = eval_active_nd(patch, xi, elem.active_indices)
                B_nd = bernstein_nd(degrees, [t_u, t_v])
                N_reconstructed = elem.C @ B_nd

                assert np.allclose(N_active_ref, N_reconstructed, atol=1e-11), (
                    f"Elem {elem.elem_index}, xi={xi}:\n"
                    f"  ref      = {N_active_ref}\n"
                    f"  C @ B_nd = {N_reconstructed}"
                )

    print("[PASS] test_nd_reconstruction_2d")


def test_nd_single_element_identity():
    """Single-element 2D patch: C must be the identity."""
    patch = build_patch(
        degrees=[1, 2],
        kvs=[[0,0, 1,1], [0,0,0, 1,1,1]],
        local_shape=[2, 3]
    )
    elems = BezierExtractor.extract_nd(patch)
    assert len(elems) == 1
    n_local = 2 * 3
    assert np.allclose(elems[0].C, np.eye(n_local), atol=1e-14)
    assert elems[0].elem_index == [0, 0]
    print("[PASS] test_nd_single_element_identity")


# ---------------------------------------------------------------------------

if __name__ == '__main__':
    test_nd_element_count_2d()
    test_nd_element_count_2d_multi()
    test_nd_matrix_shape()
    test_nd_active_indices_count()
    test_nd_column_sum_unity()
    test_nd_nonnegative()
    test_nd_elem_index_ufastest()
    test_nd_reconstruction_2d()
    test_nd_single_element_identity()
    print("\nAll ND BezierExtractor tests passed!")
