from yeti_iga.pymfiga.common.base.enum import Constants, BoundarySide, ParametricDirection
from yeti_iga.pymfiga.common.base.math import eval_inverse_and_determinant
from yeti_iga.pymfiga.common.numerics.quadrature_rules import (
    IGAQuadratureRule,
    StandardGauss,
    WeightedQuadrature,
)
from yeti_iga.pymfiga.common.numerics.solvers import NonLinearSolver
from .base import BasePatchIGA
from .math import MathMixin
from typing import List, Tuple, Optional, Sequence, Literal, Union
from geomdl import BSpline, NURBS
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.GEOMETRY")


class ProjectionMixin(BasePatchIGA):
    def __init__(
        self,
        obj: Union[
            BSpline.Curve,
            BSpline.Surface,
            BSpline.Volume,
            NURBS.Curve,
            NURBS.Surface,
            NURBS.Volume,
        ],
        quadclass: Literal["gs", "wq"],
        **kwargs
    ):
        super().__init__(obj)
        self.nbqp: np.ndarray = np.ones(self._ndim_default, dtype=int)
        self.quadrule_list: List[IGAQuadratureRule] = []
        self._set_quadrature_rules(quadclass, **kwargs)

    @property
    def nbqp_total(self) -> int:
        "Total number of quadrature points (product over all parametric directions)."
        return int(np.prod(self.nbqp[: self.ndim]))

    def _set_quadrature_rules(self, quadclass: Literal["gs", "wq"], **kwargs):
        self.quadrule_list.clear()
        assert quadclass in ["gs", "wq"], NotImplementedError()
        quadrule_class = {"gs": StandardGauss, "wq": WeightedQuadrature}[quadclass]
        if "quadtype" in kwargs:
            quadrule_type = kwargs["quadtype"]
            kwargs.pop("quadtype")
        else:
            quadrule_type = {"gs": "legendre", "wq": "2"}[quadclass]
        if "is_periodic" in kwargs:
            is_periodic = kwargs["is_periodic"]
            kwargs.pop("is_periodic")
            assert isinstance(is_periodic, (list, tuple)) and all(
                isinstance(x, bool) for x in is_periodic
            )
        else:
            is_periodic = [False for _ in range(self.ndim)]

        for i in range(self.ndim):
            q: IGAQuadratureRule = quadrule_class(
                int(self.degree[i]),
                self.knotvector[i],
                quadtype=quadrule_type,
                is_periodic=is_periodic[i],
                quad_args=kwargs,
            )
            q.export_quadrature_rules()
            self.nbqp[i] = int(q.nbquadpts)
            self.quadrule_list.append(q)

    def compute_global_mesh_parameter(self) -> Tuple[float, float]:
        """Return the maximum unique-knot spacing over all parametric directions."""
        max_distance = []
        for i in range(self.ndim):
            q = self.quadrule_list[i]
            max_distance.append(q.max_h_size)
        max_h = float(np.max(max_distance)) if max_distance else 0.0
        min_p = float(np.min(self.degree[: self.ndim]))
        return max_h, min_p

    def interpolate_field(
        self,
        knots_list: Optional[List[np.ndarray]] = None,
        u_ctrlpts: Optional[np.ndarray] = None,
        eval_geometry: bool = True,
    ) -> Sequence[Optional[np.ndarray]]:
        """Interpolates a scalar/vector field defined by control point values.

        Args:
            knots_list: custom parametric grid(s) to evaluate on (or None -> quadrature grid).
            u_ctrlpts: array of control point values, shape (ncomp, nbctrlpts_total) or (nbctrlpts_total,)
            eval_geometry: if True, also evaluate detJ and physical coordinates.

        Returns:
            (u_interp, pts_phy, det_jac)
        """
        u_interp: Optional[np.ndarray] = None
        pts_phy: Optional[np.ndarray] = None
        det_jac: Optional[np.ndarray] = None
        jac: Optional[np.ndarray] = None

        if knots_list is None:
            knots_list = [quadrule.knots_to_sample for quadrule in self.quadrule_list]

        for q, knots in zip(self.quadrule_list, knots_list):
            q.knots_to_sample = np.asarray(knots)

        if isinstance(u_ctrlpts, np.ndarray):
            u_interp = self.operator_engine.interpolate_meshgrid(
                self.quadrule_list,
                u_ctrlpts,
                nurbs_weights=self.nurbs_weights,
            )

        if eval_geometry:
            det_jac, jac, _, pts_phy = MathMixin.eval_transformation(
                self.quadrule_list,
                self.ctrlpts,
                nurbs_weights=self.nurbs_weights,
            )

        for q in self.quadrule_list:
            q.clear_sample()

        return u_interp, pts_phy, det_jac, jac

    def inverse_projection_on_boundary(
        self,
        X_target: np.ndarray,
        boundary_id: Tuple[ParametricDirection, BoundarySide],
        tolerance: float = Constants.TINY,
        method: Literal["picard", "newton"] = "newton",
        warping: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Find the parametric coordinates corresponding to a target physical point on a specified boundary.

        Args:
            X_target (np.ndarray): Target physical coordinates, shape (spatial_dim, n_points).
            boundary_id (Tuple[ParametricDirection, BoundarySide]):
                Tuple specifying the parametric direction and boundary side.
            tolerance (float): threshold for the convergence criteria in nonlinear solver.
            method (Literal["picard", "newton"]): Picard is a fixed point algorithm whereas Newton computes
                consistent tangent matrix (it requires that the second derivative of splines are defined).
            warping (Optional[np.ndarray]): displacement offset to add to control points. Default None.

        Exceptions:
            ValueError: if the nonlinear solver fails to converge or if the solution is out of bounds.
        """
        assert self.ndim >= 2
        assert (
            isinstance(X_target, np.ndarray)
            and X_target.ndim == 2
            and X_target.shape[0] >= 2
        )
        n_ctrlpts = X_target.shape[1]

        logger.info("Computing projection")
        nonlinearsolver = NonLinearSolver(
            tolerance=tolerance,
            maxiters=20,  # TODO: let the user defined it ?
            allow_acceleration=True if method == "picard" else False,
            allow_line_search=True,
            verbose=True,
        )

        solution = (
            np.zeros((self.ndim, n_ctrlpts))
            if boundary_id[1].value == 0
            else np.ones((self.ndim, n_ctrlpts))
        )
        active_params = [i for i in range(self.ndim) if i != boundary_id[0].value]
        solution[active_params] = np.ones((self.ndim - 1, n_ctrlpts)) * 0.5
        projector = ClosestPointProjector(boundary_id[0], self.ndim, method=method)
        nonlinearsolver.solve(
            solution,
            compute_residual=projector.compute_residual,
            compute_increment=projector.compute_increment,
            residual_args={
                "X_target": X_target,
                "patch": self,
                "warping": warping,
            },
        )
        if not (
            np.all((solution >= -Constants.TINY) & (solution <= 1 + Constants.TINY))
        ):
            raise ValueError("Point out of bounds")
        return solution


class ClosestPointProjector:

    def __init__(
        self,
        parame_dir: ParametricDirection,
        ndim: int,
        method: Literal["picard", "newton"] = "newton",
    ):
        self.ndim = ndim
        self.parame_dir = parame_dir.value
        self.active_params = [i for i in range(ndim) if i != self.parame_dir]
        self.order_params = self.active_params + [self.parame_dir]
        self.method = method

    def _set_knots_to_evaluate(self, solution: np.ndarray):
        Xi_k = np.clip(solution, -Constants.TINY, 1 + Constants.TINY)
        extra_knot = float(np.min(Xi_k[self.parame_dir]))
        knot_list = (Xi_k[self.active_params]).tolist() + [[extra_knot]]
        knot_list = [knot_list[i] for i in self.order_params]
        return knot_list

    def compute_residual(
        self, solution: np.ndarray, **res_args
    ) -> Tuple[np.ndarray, dict]:
        # Read inputs
        X_target: np.ndarray = res_args["X_target"]
        patch: ProjectionMixin = res_args["patch"]
        warping: Optional[np.ndarray] = res_args["warping"]

        # Preprocessing
        knot_list = self._set_knots_to_evaluate(solution)
        ctrlpts = patch.ctrlpts.copy()
        if isinstance(warping, np.ndarray):
            ctrlpts += warping

        # Computations
        for q, knots in zip(patch.quadrule_list, knot_list):
            q.knots_to_sample = np.asarray(knots)

        X_k = patch.operator_engine.interpolate_meshgrid(
            patch.quadrule_list,
            ctrlpts,
            nurbs_weights=patch.nurbs_weights,
        )  # (x_i)
        J_k = patch.operator_engine.eval_jacobien(
            patch.quadrule_list,
            ctrlpts,
            nurbs_weights=patch.nurbs_weights,
        )  # (d x_i / d xi_j)

        res_k = np.einsum("ik,ijk->jk", X_target - X_k, J_k, optimize=True)
        tan_k = np.einsum("lik,ljk->ijk", J_k, J_k, optimize=True)

        if self.method == "newton" and np.min(patch.degree) > 1:
            H_k = patch.operator_engine.eval_hessian(
                patch.quadrule_list,
                ctrlpts,
                nurbs_weights=patch.nurbs_weights,
            )  # (d^2 x_i / d xi_j /d xi_k)
            tan_k -= np.einsum("il,ijkl->jkl", X_target - X_k, H_k, optimize=True)

        for q in patch.quadrule_list:
            q.clear_sample()

        return res_k[self.active_params], {
            "tangent": tan_k[np.ix_(self.active_params, self.active_params)],
        }

    def compute_increment(self, residual: np.ndarray, **kwargs) -> np.ndarray:
        tangent: np.ndarray = kwargs["tangent"]
        invtang = eval_inverse_and_determinant(tangent)[-1]
        incr = np.zeros((self.ndim, residual.shape[1]))
        incr[self.active_params] = np.einsum(
            "ijk,jk->ik", invtang, residual, optimize=True
        )
        return incr
