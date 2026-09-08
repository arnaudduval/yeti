from yeti_iga.pymfiga.common.base.enum import Constants
from typing import Tuple, List
from scipy import sparse as sp
import numpy as np
from geomdl import helpers


class Operations:
    @staticmethod
    def find_multiplicity(
        knotvector: np.ndarray, knot: float, threshold: float = Constants.TINY
    ) -> int:
        """
        Finds the multiplicity of a knot in the knot-vector within the given threshold.
        Args:
            knotvector (np.ndarray):
            knot (float):
        """
        return np.add.reduce(np.abs(knotvector - knot) <= threshold)

    @staticmethod
    def increase_multiplicity_to_knotvector(
        repeat: int, degree: int, knotvector: np.ndarray
    ) -> np.ndarray:
        """
        Returns a knot-vector with higher multiplicity than the reference.
        If a knot has a multiplicity m in current knotvector, it will have
        a multiplicity of m + r (repeat) in output
        Args:
            repeat (int): the number to add to current multiplicity
            degree (int): degree of the spline
            knotvector (np.ndarray): the reference knot-vector
        """
        if len(knotvector) < 2 * (degree + 1):
            raise ValueError("Knot vector must contain at least 2*(degree + 1) knots")
        kv_out = list(knotvector[: degree + 1]) + list(knotvector[-(degree + 1) :])
        kv_unique = np.unique(knotvector[degree + 1 : -(degree + 1)])

        for knot in kv_unique:
            m = Operations.find_multiplicity(knotvector, knot) + repeat
            if m > degree + 1:
                raise ValueError("Introduces numerical instability")
            kv_out.extend([knot] * m)

        return np.sort(np.array(kv_out))

    @staticmethod
    def eval_ders_basis_COO_format(
        degree: int,
        knotvector: np.ndarray,
        knots: np.ndarray,
        nders: int = 1,
        is_periodic: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluates B-spline functions and its first (nders)-derivatives at given knots.
        It returns the COO format of matrices. Let say that there are n spline functions
        and m knots, each matrix should have a size of n x m.
        Args:
            degree (int): degree of the spline
            knotvector (np.ndarray): the knot-vector that defines the spline
            knots (np.ndarray): the list of knots where the spline is evaluated
            nders (int): number of derivatives to compute. By default 1
            is_periodic (bool): True if the spline is periodic, False otherwise.
                By default False
        """
        assert degree >= 0, "Degree must be a positive integer"

        nbctrlpts = len(knotvector) - degree - 1
        basis, indices_i, indices_j = [], [], []
        ring_size = nbctrlpts - degree if is_periodic else nbctrlpts

        for j, knot in enumerate(knots):
            knot_span = helpers.find_span_linear(degree, knotvector, nbctrlpts, knot)
            basis_ders = helpers.basis_function_ders(
                degree, knotvector, knot_span, knot, nders
            )

            for i, output in enumerate(zip(*basis_ders)):
                basis.append(output)
                indices_i.append((knot_span - degree + i) % ring_size)
                indices_j.append(j)

        return np.array(basis), np.array(indices_i), np.array(indices_j)

    @staticmethod
    def eval_ders_basis_sparse(
        degree: int,
        knotvector: np.ndarray,
        knots: np.ndarray,
        nders: int = 1,
        is_periodic: bool = False,
    ) -> List[sp.csr_array]:
        """
        Evaluates B-spline functions and its first (nders)-derivatives at given knots.
        It returns the list of CSR matrices. Let say that there are n spline functions
        and m knots, each matrix should have a size of n x m.
        Args:
            degree (int): degree of the spline
            knotvector (np.ndarray): the knot-vector that defines the spline
            knots (np.ndarray): the list of knots where the spline is evaluated
            nders (int): number of derivatives to compute. By default 1
            is_periodic (bool): True if the spline is periodic, False otherwise.
                By default False
        """
        basis_coo, indi_coo, indj_coo = Operations.eval_ders_basis_COO_format(
            degree,
            knotvector,
            knots,
            nders=nders,
            is_periodic=is_periodic,
        )

        # Loop through the basis functions and derivatives
        nbctrlpts = len(knotvector) - degree - 1
        nbctrlpts = nbctrlpts - degree if is_periodic else nbctrlpts
        basis_list = []
        for i in range(basis_coo.shape[1]):
            b = sp.coo_array(
                (basis_coo[:, i], (indj_coo, indi_coo)), shape=(len(knots), nbctrlpts)
            ).tocsr()
            b.eliminate_zeros()
            basis_list.append(sp.csr_array(b))
        return basis_list
