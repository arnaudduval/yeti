from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.numerics.operations import BsplineOperations, NurbsOperations
from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from yeti_iga.pymfiga.common.base.math import eval_inverse_and_determinant
from typing import List, Tuple
import numpy as np


class MathMixin:
    @staticmethod
    def inverse_rectangular_matrix(jac: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the Moore-Penrose pseudoinverse of a rectangular matrix A (here calle jac).

        Args:
            jac (np.ndarray): Input matrix of shape (m, n, ...)

        Returns:
            det_jac (np.ndarray): Determinant-like of jac, shape (...)
            inv_jac (np.ndarray): Pseudoinverse of jac, shape (n, m, ...)
        """
        # Metric tensor G = J^T J over param space
        G = np.einsum("li...,lj...->ij...", jac, jac, optimize=True)
        det_G, inv_G = eval_inverse_and_determinant(G)
        det_jac = np.sqrt(det_G)
        # J^{-1}_phys = inv(J^T J) J^T
        inv_jac = np.einsum("il...,jl...->ij...", inv_G, jac, optimize=True)
        return det_jac, inv_jac

    @staticmethod
    def eval_transformation(
        quadrule_list: List[IGAQuadratureRule],
        ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute |J|, J, J^{-1} and physical coordinates at quadrature points.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                for each parametric direction.
            ctrlpts (np.ndarray): control points array.
            nurbs_weights (np.ndarray): optional array of weights for NURBS.

        Returns:
            (det_jac, jac, inv_jac, knots_phy)

        Exceptions:
            ValueError: if any quadrature point has non-positive Jacobian determinant.
        """
        operator_engine = (
            NurbsOperations if nurbs_weights.size > 0 else BsplineOperations
        )

        jac = operator_engine.eval_jacobien(
            quadrule_list, ctrlpts, nurbs_weights=nurbs_weights
        )  # shape: (spatial_dim, param_dim, ...)

        det_jac, inv_jac = MathMixin.inverse_rectangular_matrix(jac)

        knots_phy = operator_engine.interpolate_meshgrid(
            quadrule_list, ctrlpts, nurbs_weights=nurbs_weights
        )

        if not np.all(det_jac > 0.0):
            nb_bad = int(np.sum(det_jac <= 0.0))
            frac = 100.0 * nb_bad / det_jac.size
            msg = (
                f"Potential geometry issue: {frac:.2e}% quadrature points have detJ <= 0 "
                f"({nb_bad}/{det_jac.size})."
            )
            raise ValueError(msg)

        return det_jac, jac, inv_jac, knots_phy

    @staticmethod
    def detect_true_dimension(ctrlpts: np.ndarray, tol: float = Constants.TINY) -> int:
        """
        Detects the true geometric dimension of control points by checking
        which coordinates actually vary.
        """
        ndim_detected = 0
        for i in range(ctrlpts.shape[0]):  # loop over x, y, z
            if np.max(ctrlpts[i]) - np.min(ctrlpts[i]) > tol:
                ndim_detected += 1
        return ndim_detected
