from typing import List, Optional
import numpy as np


def kron(arrays: List[np.ndarray]) -> np.ndarray:
    a = arrays[0]
    for curr in arrays[1:]:
        a = np.kron(curr, a)
    return a


def combine_arrays(
    arr: List[np.ndarray], coefs: np.ndarray, brr: Optional[List[np.ndarray]] = None
) -> np.ndarray:
    ndim = len(arr)
    assert ndim > 0 and len(coefs) == ndim
    if brr is None:
        brr = [np.ones_like(a) for a in arr]

    mixed = [[brr[i].copy() for i in range(ndim)].copy() for _ in range(ndim)]
    for i in range(ndim):
        mixed[i][i] = arr[i] * coefs[i]

    v_out = kron(mixed[0])
    for i in range(1, ndim):
        v_out += kron(mixed[i])
    return v_out


def compute_special_schur(
    H: np.ndarray,
    bfac: np.ndarray,
) -> np.ndarray:
    """
    H     : (k+1, m) diagonals of H_i, i = 1, ..., k+1
    bfac  : (k,)  scalars b_i
    """

    H_top = H[:-1, :]  # shape k, m
    H_last = H[-1, :]  # shape m

    # Schur complement
    schur = H_last + np.sum((np.abs(bfac) ** 2)[:, None] / H_top, axis=0)

    return schur


def solve_special_arrowhead(
    H: np.ndarray,
    bfac: np.ndarray,
    rhs: np.ndarray,
    schur: Optional[np.ndarray] = None,
):
    """
    H     : (k+1, m) diagonals of H_i, i = 1, ..., k+1
    bfac  : (k,)  scalars b_i
    rhs   : (r, (k+1)*m) RHS
    """

    H_top = H[:-1, :]  # shape k, m
    r = rhs.shape[0]

    rhs = np.reshape(rhs, (r, *H.shape))  # shape r, k+1, m
    rhs_top = rhs[:, :-1, :]  # shape r, k, m
    rhs_last = rhs[:, -1, :]  # shape r, m

    # y_i = H_i^{-1} rhs_i
    y = rhs_top / H_top[None, :, :]  # shape r, k, m

    # Schur complement
    if schur is None:
        schur = compute_special_schur(H, bfac)  # shape m

    # RHS of last block: shape r, m
    rhs_s = rhs_last + np.sum(np.conj(bfac)[None, :, None] * y, axis=1)

    x_last = rhs_s / np.asarray(schur)[None, :]  # shape r, m

    # back substitution
    multp = bfac[:, None] / H_top  # shape k, m
    x_top = y - multp[None, :, :] * x_last[:, None, :]  # shape r, k, m

    x = np.concatenate([x_top, x_last[:, None, :]], axis=1)

    return np.reshape(x, (r, -1))  # shape r, (k+1)*m
