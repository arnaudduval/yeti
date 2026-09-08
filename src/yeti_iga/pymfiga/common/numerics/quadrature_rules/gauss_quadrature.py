from .operations import Operations
from .core import IGAQuadratureRule
from .data import legendre_table, lobatto_table
from dataclasses import dataclass
from typing import Tuple, Literal, Dict
import numpy as np
import logging

logger = logging.getLogger("SRC.QUADRATURE")


@dataclass
class StandardGaussArgs:
    quadperel: int

    def __post_init__(self):
        if self.quadperel <= 0:
            raise ValueError(
                f"Number of quadrature points must be positive, got {self.quadperel}"
            )


class StandardGauss(IGAQuadratureRule):
    def __init__(
        self,
        degree: int,
        knotvector: np.ndarray,
        quadtype: Literal["legendre", "lobatto"],
        is_periodic: bool = False,
        quad_args: Dict[str, int] = {},
    ):
        """
        Parameters:

            degree: int, polynomial degree of spline
            knotvector: np.ndarray, knot-vector that defines the spline
            quadtype: Literal["legendre", "lobatto"], for Gauss-Legendre or Gauss-Lobatto
            is_periodic: bool, True if spline is periodic, False otherwise.
                by default is set to False
            quad_args: dict, it has the information needed to
                build the quadrature rule.

        Notes:
        ------
                quad_args may containt the parameter
                - quadperel: int, the number of quadrature points per element, default degree + 1
        """
        super().__init__(degree, knotvector, is_periodic=is_periodic)
        self._quadrature_class = "GAUSS"
        self._set_quadrature_type(quadtype)

        # Public variables
        self._parametric_weights: np.ndarray = np.array([])

        # Private variables
        self._isoparametric_positions = None
        self._isoparametric_weights = None
        self._set_isoparametric_info(quad_args)

    @property
    def parametric_weights(self):
        "Weights in the reference space [0, 1]"
        return self._parametric_weights

    def _set_quadrature_type(self, value):
        assert value in ["lobatto", "legendre"]
        self._quadrature_type = value

    def _set_isoparametric_info(self, quad_args: dict):
        val = 1 if self.quadrature_type == "legendre" else 2
        input_args = {"quadperel": self.degree + val}
        if "quadperel" in quad_args:
            input_args["quadperel"] = quad_args["quadperel"]
        args = StandardGaussArgs(**input_args)
        table = {"legendre": legendre_table, "lobatto": lobatto_table}[
            self.quadrature_type
        ]
        self._isoparametric_positions, self._isoparametric_weights = table(
            args.quadperel
        )

    def _set_quadrature_points(self):
        knots = self.unique_kv
        quadpts = np.concatenate(
            [
                0.5
                * (
                    (knots[i + 1] - knots[i]) * self._isoparametric_positions
                    + knots[i]
                    + knots[i + 1]
                )
                for i in range(self.nbelem)
            ]
        )
        self.quadpts = quadpts
        if self.quadrature_type == "lobatto" and self.degree == 1:
            logger.warning(
                "Becarefull, Gauss-Lobatto quadrature is not"
                "supposed to work out with linear polynomials."
            )

    def compute_parametric_weights(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        "Computes the basis and weights of the quadrature class (COO format)"
        knots = self.unique_kv
        self._parametric_weights = np.concatenate(
            [
                0.5 * (knots[i + 1] - knots[i]) * self._isoparametric_weights
                for i in range(self.nbelem)
            ]
        )
        basis, indi, indj = Operations.eval_ders_basis_COO_format(
            self.degree,
            self.knotvector,
            self.quadpts,
            is_periodic=self.is_periodic,
        )
        nnz = np.shape(basis)[0]
        weights = np.zeros((nnz, 4))
        weights[:, 0] = basis[:, 0] * self.parametric_weights[indj]
        weights[:, 3] = basis[:, 1] * self.parametric_weights[indj]
        weights[:, 1] = weights[:, 0]
        weights[:, 2] = weights[:, 3]
        return basis, weights, indi, indj

    def export_quadrature_rules(self):
        self._set_quadrature_points()
        basis, weights, indi, indj = self.compute_parametric_weights()
        super()._set_coo_basis_weights(basis, weights, [indi, indj])
        super()._assemble_csr_basis_weights()
        logger.info(repr(self))
        return self
