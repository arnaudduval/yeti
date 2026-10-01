"""
matrix_free_apply_2d: 2D port of pymfiga's MatrixFree.apply (Kronecker
sandwich product without ever forming the Kronecker product matrix). Pure
linear-algebra identity, cross-checked against scipy.sparse.kron -- plus one
check using real WeightedQuadrature1D matrices, the actual intended use
case (Phase C: matrix-free elasticity).
"""

import numpy as np
import pytest
import scipy.sparse as sp

from yeti_iga.future.bspline import BSpline, WeightedQuadrature1D, matrix_free_apply_2d


@pytest.mark.parametrize("sparse_format", ["csc", "csr"])
def test_matrix_free_apply_2d_matches_kronecker(sparse_format):
    rng = np.random.default_rng(0)
    nu_out, nu_in = 5, 4
    nv_out, nv_in = 3, 6
    Mu = sp.random(nu_out, nu_in, density=0.6, random_state=rng, format=sparse_format)
    Mv = sp.random(nv_out, nv_in, density=0.6, random_state=rng, format=sparse_format)
    v_in = rng.standard_normal(nu_in * nv_in)

    result = matrix_free_apply_2d(Mu, Mv, v_in, False)
    expected = sp.kron(Mv, Mu) @ v_in

    assert np.allclose(result, expected, atol=1e-10)


def test_matrix_free_apply_2d_transpose_matches_kronecker():
    rng = np.random.default_rng(1)
    nu_out, nu_in = 5, 4
    nv_out, nv_in = 3, 6
    Mu = sp.random(nu_out, nu_in, density=0.6, random_state=rng, format="csc")
    Mv = sp.random(nv_out, nv_in, density=0.6, random_state=rng, format="csc")
    v_out = rng.standard_normal(nu_out * nv_out)

    result = matrix_free_apply_2d(Mu, Mv, v_out, True)
    expected = sp.kron(Mv, Mu).T @ v_out

    assert np.allclose(result, expected, atol=1e-10)


def test_matrix_free_apply_2d_rejects_size_mismatch():
    Mu = sp.eye(3, format="csc")
    Mv = sp.eye(4, format="csc")
    with pytest.raises(Exception):
        matrix_free_apply_2d(Mu, Mv, np.zeros(5), False)


def test_matrix_free_apply_2d_with_weighted_quadrature_matrices():
    """
    The actual intended use case: gather a 2D tensor-product field to WQ
    points (B0_u (x) B0_v), scatter back (W00_u (x) W00_v), and reuse B0 as
    its own scatter via is_transpose -- each checked against explicit
    Kronecker assembly.
    """
    degree = 2
    kv = np.array([0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1], dtype=float)
    bsp = BSpline(degree, kv)
    wq = WeightedQuadrature1D.build(bsp)

    nbctrlpts = wq.nbctrlpts
    nq = len(wq.quadpts)
    rng = np.random.default_rng(2)

    # Gather: control-point field -> values at WQ points (u and v directions
    # using the same 1D rule here, for simplicity).
    dofs = rng.standard_normal(nbctrlpts * nbctrlpts)
    gathered = matrix_free_apply_2d(wq.B0, wq.B0, dofs, False)
    expected_gathered = sp.kron(wq.B0, wq.B0) @ dofs
    assert np.allclose(gathered, expected_gathered, atol=1e-10)
    assert gathered.shape == (nq * nq,)

    # Scatter: values at WQ points -> control-point field, using W00
    # directly (W00 already has shape nbctrlpts x nq, the opposite
    # convention from B0, so no transpose needed here).
    point_values = rng.standard_normal(nq * nq)
    scattered = matrix_free_apply_2d(wq.W00, wq.W00, point_values, False)
    expected_scattered = sp.kron(wq.W00, wq.W00) @ point_values
    assert np.allclose(scattered, expected_scattered, atol=1e-10)
    assert scattered.shape == (nbctrlpts * nbctrlpts,)

    # is_transpose=True lets the SAME gather matrix (B0, shape nq x
    # nbctrlpts) double as a scatter, without storing a separately
    # transposed copy -- the common B0^T @ (...) Galerkin-assembly pattern.
    scattered_via_transpose = matrix_free_apply_2d(wq.B0, wq.B0, point_values, True)
    expected_via_transpose = sp.kron(wq.B0, wq.B0).T @ point_values
    assert np.allclose(scattered_via_transpose, expected_via_transpose, atol=1e-10)
    assert scattered_via_transpose.shape == (nbctrlpts * nbctrlpts,)
