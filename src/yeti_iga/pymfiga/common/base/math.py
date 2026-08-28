from .enum import Constants
from typing import Tuple, List
from functools import reduce
from operator import mul
import numpy as np
import itertools


def eval_inverse_and_determinant(matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Computes the determinant and the inverse of a square matrix (n-array)

    Args:
        matrix (np.ndarray): matrix of size m x m x ..., with m <= 3

    Returns:
        A tuple with the determinant and inverse of the matrix

    Raises:
        NotImplementedError: if m > 3 with m being the size of the matrix
        ZeroDivisionError: if any value of the determinant is close to zero
    """
    assert isinstance(matrix, np.ndarray) and matrix.ndim > 2
    m1, m2 = np.shape(matrix)[:2]
    assert m1 == m2
    other_dim = matrix.shape[2:]
    inv = np.zeros(shape=(m2, m1) + other_dim)
    if m1 == 1:
        det = matrix[0, 0, ...].copy()
        inv = np.ones((1, 1, *np.shape(det)))
    elif m1 == 2:
        det = (
            matrix[0, 0, ...] * matrix[1, 1, ...]
            - matrix[0, 1, ...] * matrix[1, 0, ...]
        )
        inv[0, 0, ...] = matrix[1, 1, ...]
        inv[1, 1, ...] = matrix[0, 0, ...]
        inv[0, 1, ...] = -matrix[0, 1, ...]
        inv[1, 0, ...] = -matrix[1, 0, ...]
    elif m1 == 3:
        det = (
            matrix[0, 1, ...] * matrix[1, 2, ...] * matrix[2, 0, ...]
            - matrix[0, 2, ...] * matrix[1, 1, ...] * matrix[2, 0, ...]
            + matrix[0, 2, ...] * matrix[1, 0, ...] * matrix[2, 1, ...]
            - matrix[0, 0, ...] * matrix[1, 2, ...] * matrix[2, 1, ...]
            + matrix[0, 0, ...] * matrix[1, 1, ...] * matrix[2, 2, ...]
            - matrix[0, 1, ...] * matrix[1, 0, ...] * matrix[2, 2, ...]
        )
        inv[0, 0, ...] = (
            matrix[1, 1, ...] * matrix[2, 2, ...]
            - matrix[1, 2, ...] * matrix[2, 1, ...]
        )
        inv[0, 1, ...] = (
            matrix[0, 2, ...] * matrix[2, 1, ...]
            - matrix[0, 1, ...] * matrix[2, 2, ...]
        )
        inv[0, 2, ...] = (
            matrix[0, 1, ...] * matrix[1, 2, ...]
            - matrix[0, 2, ...] * matrix[1, 1, ...]
        )
        inv[1, 0, ...] = (
            matrix[1, 2, ...] * matrix[2, 0, ...]
            - matrix[1, 0, ...] * matrix[2, 2, ...]
        )
        inv[1, 1, ...] = (
            matrix[0, 0, ...] * matrix[2, 2, ...]
            - matrix[0, 2, ...] * matrix[2, 0, ...]
        )
        inv[1, 2, ...] = (
            matrix[0, 2, ...] * matrix[1, 0, ...]
            - matrix[0, 0, ...] * matrix[1, 2, ...]
        )
        inv[2, 0, ...] = (
            matrix[1, 0, ...] * matrix[2, 1, ...]
            - matrix[1, 1, ...] * matrix[2, 0, ...]
        )
        inv[2, 1, ...] = (
            matrix[0, 1, ...] * matrix[2, 0, ...]
            - matrix[0, 0, ...] * matrix[2, 1, ...]
        )
        inv[2, 2, ...] = (
            matrix[0, 0, ...] * matrix[1, 1, ...]
            - matrix[0, 1, ...] * matrix[1, 0, ...]
        )
    else:
        raise NotImplementedError("Only 1x1, 2x2 and 3x3 matrices are supported")

    if np.any(np.abs(det) < Constants.SAFEGUARD):
        raise ZeroDivisionError("There are near to zero determinants")

    inv = inv / det
    return det, inv


def kron_nonzero_indices(indices_list: List[list], nnz_list: List[int]) -> list:
    """
    Finds the nonzero indices of a kronecker product of sparse arrays.
    For example, let say A_1, A_2, ..., A_n are sparse arrays, then
    the result A = A_n x ... x A_2 x A_1 is also sparse. It computes then
    the nonzero values of A knowing the nonzero values of A_1, ..., A_n.

    Args:
        indices_list (List[list]): list that constains the nonzero
            indices of the different arrays
        nnz_list (List[int]): list that contains the arrays' size

    Returns:
        A list of nonzero indices of the resulting array
    """
    strides = [reduce(mul, nnz_list[i + 1 :], 1) for i in range(len(nnz_list))]
    flat_indices = []
    for multi_idx in itertools.product(*indices_list):
        flat = sum(i * s for i, s in zip(multi_idx, strides))
        flat_indices.append(flat)
    return flat_indices
