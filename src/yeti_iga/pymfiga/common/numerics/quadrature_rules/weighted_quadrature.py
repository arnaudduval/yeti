from .operations import Operations
from .core import IGAQuadratureRule
from .gauss_quadrature import StandardGauss
from typing import Tuple, Callable, Literal
from dataclasses import dataclass
from scipy import sparse as sp
import numpy as np
import logging

logger = logging.getLogger("SRC.QUADRATURE")


def solve_optimization_problem(
    Z: np.ndarray,
    A: np.ndarray,
    B: np.ndarray,
) -> np.ndarray:
    """
    Solves the optimization problem:
        minimize ||diag(Z_i)^{-1} @  w||^2
        subject to A @ w = B_i, for i = 1, ..., n

    Parameters:
        Z (np.ndarray): Coefficient matrix of size (n, m).
        A (np.ndarray): Constraint matrix of size (p, m).
        b (np.ndarray): Right-hand side matrix of size (n, p).
        threshold (float): Threshold below which numbers are treated as zero.

    Returns:
        np.ndarray: Solution vector w of size (n, m,).
    """
    assert B.shape[0] == Z.shape[0], "No. columns in b must match the no. rows in z"
    assert B.shape[1] == A.shape[0], "b and A must have compatible dimensions"
    assert A.shape[1] == Z.shape[1], "A and Z must have compatible dimensions"
    solution_all = np.zeros_like(Z)

    for ii, (zz, bb) in enumerate(zip(Z, B)):

        solution_all[ii] = (
            np.diag(zz) @ np.linalg.lstsq(A @ np.diag(zz), bb, rcond=None)[0]
        )

    return solution_all


@dataclass
class WeightedQuadratureArgs:
    position_rule: str
    rule_parameters: dict

    def __post_init__(self):
        if self.position_rule not in ["midpoint", "internal", "external_source"]:
            raise ValueError(
                f"Position rule must be midpoint, internal or external_source, got {self.position_rule}"
            )
        if self.position_rule in ["midpoint", "internal"]:
            keys_types = [("s", int), ("r", int), ("include_boundaries", bool)]
            for ky, tp in keys_types:
                assert (
                    ky in self.rule_parameters
                ), f"{ky} must be include in rule_parameters"
                assert isinstance(self.rule_parameters.get(ky), tp)
        elif self.position_rule == "external_source":
            ky = "quadrature_points"
            assert (
                ky in self.rule_parameters
            ), f"{ky} must be include in rule_parameters"
            assert isinstance(self.rule_parameters.get(ky), np.ndarray)


class WeightedQuadrature(IGAQuadratureRule):
    def __init__(
        self,
        degree: int,
        knotvector: np.ndarray,
        quadtype: Literal["1", "2"],
        is_periodic: bool = False,
        quad_args: dict = {},
    ):
        """
        Parameters:
            degree: int, polynomial degree of spline
            knotvector: np.ndarray, knot-vector that defines the spline
            quadtype: Literal["1", "2"], for first and second implementation
            is_periodic: bool, True if spline is periodic, False otherwise.
                by default is set to False
            quad_args: dict, it has the information needed to
                build the quadrature rule

        Notes:
        ------
                quad_args may containt the parameters:
                - position_rule: Literal["midpoint", "internal", "external_source"],
                        it defines the algorithm to compute the position of quadrature points
                - rule_parameters: dict, it has information for position rule algorithm.
        """
        super().__init__(degree, knotvector, is_periodic=is_periodic)
        self._quadrature_class = "WEIGHTED"
        self._set_quadrature_type(quadtype)

        # Private variables
        self._use_other: bool = False
        self._read_data_from_args(quad_args)

    def _set_quadrature_type(self, value: str):
        assert value in ["1", "2"]
        self._quadrature_type = value

    def _read_data_from_args(self, quad_args: dict):
        output = self._get_position_rule_and_defaults(self.quadrature_type)
        default_position_rule, default_rule_param = output
        self._args = WeightedQuadratureArgs(
            position_rule=quad_args.get("position_rule", default_position_rule),
            rule_parameters=quad_args.get("rule_parameters", default_rule_param),
        )
        self._use_other = False if self.degree > 1 else True
        if self._use_other:
            logger.warning(
                "Becarefull, weighted quadrature is not"
                "supposed to work out with linear polynomials."
                "By default, Gauss quadrature will be used."
            )

    def _get_position_rule_and_defaults(self, quadrature_type: str):
        if quadrature_type == "1":
            return "midpoint", {
                "s": 1,
                "r": self.degree + 2,
                "include_boundaries": True,
            }
        elif quadrature_type == "2":
            return "midpoint", {
                "s": 2,
                "r": self.degree + 3,
                "include_boundaries": True,
            }
        # NOTE: we can add other strategies
        else:
            raise ValueError(f"Unknown quadrature type: {quadrature_type}")

    def _set_quadrature_points(self):
        rule_param = self._args.rule_parameters
        s = rule_param.get("s", 1)
        r = rule_param.get("r", self.degree + 2)
        include_boundaries = rule_param.get("include_boundaries", True)
        position_rule = self._args.position_rule
        if position_rule == "midpoint":
            quadpts = self._midpoint_rule(
                s=s, r=r, include_boundaries=include_boundaries
            )
        elif position_rule == "internal":
            quadpts = self._internal_rule(
                s=s, r=r, include_boundaries=include_boundaries
            )
        elif position_rule == "external_source":
            quadpts = rule_param.get("quadrature_points", np.array([]))
            assert quadpts is not None
        else:
            raise NotImplementedError(
                f"Unknown position rule algorithm: {position_rule}"
            )
        self.quadpts = quadpts

    def _midpoint_rule(self, s: int, r: int, include_boundaries: bool) -> np.ndarray:
        return self._generate_quadrature_points(
            s, r, np.linspace, include_boundaries=include_boundaries
        )

    def _internal_rule(self, s: int, r: int, include_boundaries: bool) -> np.ndarray:
        def get_quadrature_points_elementwise(start, end, n):
            return np.array(
                [
                    (1 - (2 * k - 1) / (2 * n)) * start + ((2 * k - 1) / (2 * n)) * end
                    for k in range(1, n + 1)
                ]
            )

        return self._generate_quadrature_points(
            s,
            r,
            get_quadrature_points_elementwise,
            include_boundaries=include_boundaries,
        )

    def _generate_quadrature_points(
        self, s: int, r: int, algorithm: Callable, include_boundaries: bool
    ) -> np.ndarray:
        quadpts = []
        knots = self.unique_kv

        # First span
        tmp = algorithm(knots[0], knots[1], r)
        if not include_boundaries:
            tmp = tmp[1:-1]
        quadpts.extend(tmp)

        # Last span
        tmp = algorithm(knots[-2], knots[-1], r)
        if not include_boundaries:
            tmp = tmp[1:-1]
        quadpts.extend(tmp)

        # Inner spans
        for i in range(1, self.nbelem - 1):
            tmp = algorithm(knots[i], knots[i + 1], 2 + s)
            if not include_boundaries:
                tmp = tmp[1:-1]
            quadpts.extend(tmp)

        return np.sort(np.unique(np.array(quadpts)))

    def _compute_knot_support(self, quadpts: np.ndarray) -> np.ndarray:
        quadpts_extended = np.concatenate(
            [np.array([-quadpts[0]]), quadpts, np.array([2 - quadpts[-1]])]
        )
        mean_quapts_extendend = (quadpts_extended[:-1] + quadpts_extended[1:]) / 2.0
        return np.diff(mean_quapts_extendend)

    def compute_parametric_weights(
        self, method: str
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Computes the basis and weights of the quadrature class (COO format)
        Args:
            method (Literal["1", "2"]): the algorithm name to compute the weights
        """

        # Test space
        gauss_test = StandardGauss(
            self.degree,
            self.knotvector,
            quadtype="legendre",
            is_periodic=self.is_periodic,
        )
        gauss_test.export_quadrature_rules()
        B0cgg_test = gauss_test.basis[0]
        W0cgg_test, _, W1cgg_test, _ = gauss_test.weights
        basis_coo, indi_coo, indj_coo = Operations.eval_ders_basis_COO_format(
            self.degree,
            self.knotvector,
            self.quadpts,
            nders=1,
            is_periodic=self.is_periodic,
        )
        B0wq_test = sp.coo_array((basis_coo[:, 0], (indi_coo, indj_coo))).tocsr()
        B1wq_test = sp.coo_array((basis_coo[:, 1], (indi_coo, indj_coo))).tocsr()

        # Target space
        if method == "1":
            # Space S^[p-1]_[r-1]
            degree_target = self.degree - 1
            knotvector_target = self.knotvector[1:-1]
        elif method == "2":
            # Space S^[p]_[r-1]
            degree_target = self.degree
            knotvector_target = Operations.increase_multiplicity_to_knotvector(
                1, degree_target, self.knotvector
            )
        else:
            raise NotImplementedError()

        B0cgg_target = Operations.eval_ders_basis_sparse(
            degree_target,
            knotvector_target,
            gauss_test.quadpts,
            nders=1,
            is_periodic=self.is_periodic,
        )[0]
        B0wq_target = (
            Operations.eval_ders_basis_sparse(
                degree_target,
                knotvector_target,
                self.quadpts,
                nders=1,
                is_periodic=self.is_periodic,
            )[0]
            .toarray()
            .T
        )

        # Compute quadrature points support
        quadpts_support: np.ndarray = self._compute_knot_support(self.quadpts)

        # Compute the weights
        list_weights: list = []

        regularization = sp.csr_array(B0wq_test @ sp.diags(quadpts_support))
        if method == "1":
            # Computation of W00
            integral = sp.csr_array(W0cgg_test @ B0cgg_test)
            list_weights.append(
                solve_optimization_problem(
                    Z=regularization.toarray(),
                    A=B0wq_test.toarray(),
                    B=integral.toarray(),
                )
            )

        # Computation of W01 for method 1 or W0 for method 2
        integral = sp.csr_array(W0cgg_test @ B0cgg_target)
        list_weights.append(
            solve_optimization_problem(
                Z=regularization.toarray(),
                A=B0wq_target,
                B=integral.toarray(),
            )
        )

        regularization = sp.csr_array(B1wq_test @ sp.diags(quadpts_support))
        if method == "1":
            # Computation of W10
            integral = sp.csr_array(W1cgg_test @ B0cgg_test)
            list_weights.append(
                solve_optimization_problem(
                    Z=regularization.toarray(),
                    A=B0wq_test.toarray(),
                    B=integral.toarray(),
                )
            )

        # Computation of W11 for method 1 or W1 for method 2
        integral = sp.csr_array(W1cgg_test @ B0cgg_target)
        list_weights.append(
            solve_optimization_problem(
                Z=regularization.toarray(),
                A=B0wq_target,
                B=integral.toarray(),
            )
        )

        weights_coo = np.zeros((np.size(basis_coo, axis=0), 4))
        list_of_indices = {"1": [0, 1, 2, 3], "2": [0, 0, 1, 1]}[method]
        for idx, (i, j) in enumerate(zip(indi_coo, indj_coo)):
            weights_coo[idx] = [
                list_weights[list_of_indices[k]][i, j]
                for k in range(len(list_of_indices))
            ]

        return basis_coo, weights_coo, indi_coo, indj_coo

    def export_quadrature_rules(self):
        self._set_quadrature_points()
        quadtype = self.quadrature_type
        if self._use_other:
            # Overwrite the quadrature points, basis and weights with other quadrature
            value = 1 if quadtype == "1" else 2
            quadrule = StandardGauss(
                degree=self.degree,
                knotvector=self.knotvector,
                quadtype="legendre",
                is_periodic=self.is_periodic,
                quad_args={
                    "quadperel": (self.degree + value),
                },
            )
            quadrule.export_quadrature_rules()
            self.quadpts = quadrule.quadpts
            basis, weights, indi, indj = quadrule.compute_parametric_weights()

        else:
            basis, weights, indi, indj = self.compute_parametric_weights(
                method=quadtype
            )
        super()._set_coo_basis_weights(basis, weights, [indi, indj])
        super()._assemble_csr_basis_weights()
        logger.info(repr(self))
        return self
