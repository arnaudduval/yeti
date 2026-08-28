from typing import Union, Sequence
from scipy import sparse as sp
import numpy as np


class MatrixFree:

    # NOTE: We define the max number of matrices as a class property
    # Since we are in mechanics, it is equal to 4: X, Y, Z, T (space and time)
    # WARNING: high dimensionality could lead to expensive CPU time in reshaping

    tensor_ndim = 4

    @staticmethod
    def _tensor_matrix_product(
        tensor: np.ndarray,
        matrix: Union[sp.csr_array, np.ndarray],
        mode: int,
        is_transpose: bool = False,
    ) -> np.ndarray:

        mat = matrix.T if is_transpose else matrix

        if sp.issparse(mat):
            tensor_ndim = MatrixFree.tensor_ndim
            old_shape = tensor.shape
            tensor_perm = np.moveaxis(tensor, mode, 0)
            tensor_reshaped = np.reshape(tensor_perm, (old_shape[mode], -1))
            new_tensor = mat @ tensor_reshaped
            new_shape = [mat.shape[0]] + [
                old_shape[i] for i in range(tensor_ndim) if i != mode
            ]
            new_tensor_reshaped = np.reshape(new_tensor, new_shape)
            return np.moveaxis(new_tensor_reshaped, 0, mode)

        elif isinstance(mat, np.ndarray):
            new_tensor = np.tensordot(tensor, mat, axes=(mode, 1))
            return np.moveaxis(new_tensor, -1, mode)

        else:
            raise NotImplementedError("Method only for scipy.sparse and numpy.ndarray")

    @staticmethod
    def apply(
        matrix_list: Sequence[Union[sp.csr_array, np.ndarray]],
        array_in: np.ndarray,
        is_transpose: bool = False,
    ) -> np.ndarray:
        """
        Computes the matrix-free product M @ v, with M being the result of
        (M_n x ... x M_2 x M_1), where 'x' represents Kronecker product and
        M_i are 2-dimensional matrices.
        Note: It adapts the tensor-matrix product for dense and sparse matrices.
        Property: (M_n x ... x M_2 x M_1) . v = V x M_n x_n M_(n-1) x_(n-1) ... x_1 M_1
        where v is a ravel of V and x_i is the i-mode tensor product

        Args:
            matrix_list (Sequence[Union[np.ndarray, sp.spmatrix]]):
                a list of sparse of dense matrices [M_1, M_2, ..., M_n]
            array_in (np.ndarray): the vector to be multiplied with the resulting matrix
            is_transpose (bool): if true it performs transpose(M) @ v, if false M @ v.
                By default is set to false.
        Returns:
            np.ndarray, the result of the matrix-vector product
        """
        tensor_dim = MatrixFree.tensor_ndim
        nmodes = len(matrix_list)

        # We reverse the order to match Mn x ... x M1
        matrices = matrix_list[::-1]
        nbcols = [matrix.shape[0 if is_transpose else 1] for matrix in matrices]
        original_shape = [nbcols[i] if i < nmodes else 1 for i in range(tensor_dim)]

        tensor = np.reshape(array_in, original_shape)
        for i, matrix in zip(range(tensor_dim), matrices):
            # NOTE: the mode is i and not (nmodes - i - 1) due to C-ordering
            # resulting on a match with Kronecker product. Otherwise, we should
            # reshape the tensor using F-ordering which is not native in python
            # leading to unnecesary copies to access memory.
            tensor = MatrixFree._tensor_matrix_product(
                tensor, matrix, mode=i, is_transpose=is_transpose
            )
        return np.ravel(tensor)
