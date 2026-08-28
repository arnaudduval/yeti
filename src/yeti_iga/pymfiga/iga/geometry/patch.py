from yeti_iga.pymfiga.common.base.enum import ParametricDirection
from yeti_iga.pymfiga.common.numerics.quadrature_rules import StandardGauss
from .mixin.math import MathMixin
from .mixin.projection import ProjectionMixin
from .mixin.transformation import TransformationMixin
from typing import Union, List, Literal, Dict
from geomdl import BSpline, NURBS
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.GEOMETRY")


class SinglePatch(MathMixin, ProjectionMixin, TransformationMixin):
    """Extracts essential info from a geomdl object (Curve/Surface/Volume)
    and prepares quadrature + geometric mappings for IGA.

    Notes / key design choices:
    - `_ndim_default` stays 3 (XYZ container) but results are sliced to `ndim` **always**.
        This was previously gated by `_is_trimmed`, which was misleading.
    - Control points & weights are read via `ctrlpts2d/3d` (and weights analogues) when available;
        otherwise we fall back to a robust flattening with documented u-fastest ordering.
    """

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
        **kwargs,
    ):
        MathMixin.__init__(self)
        ProjectionMixin.__init__(self, obj, quadclass, **kwargs)
        TransformationMixin.__init__(self, obj)

    @property
    def dir_degrees(self):
        return {ParametricDirection(i): self.degree[i] for i in range(self.ndim)}

    @property
    def dir_knotvectors(self):
        return {ParametricDirection(i): self.knotvector[i] for i in range(self.ndim)}

    @property
    def dir_is_periodic(self):
        return {
            ParametricDirection(i): self.quadrule_list[i].is_periodic
            for i in range(self.ndim)
        }

    def _get_quadrature_for_boundary(
        self,
        boundary_dir: ParametricDirection,
        p_master_degrees: Dict[ParametricDirection, int],
    ) -> List[StandardGauss]:
        integration_dirs = ParametricDirection.get_integration_dirs(
            boundary_dir, self.ndim
        )
        all_quadrules = []
        for direction in integration_dirs:
            p_S = self.dir_degrees.get(direction, 1)
            p_M = p_master_degrees.get(direction, 1)
            nb_pts_per_el = int(max(p_S, p_M)) + 1
            quadrature = StandardGauss(
                self.dir_degrees[direction],
                self.dir_knotvectors[direction],
                quadtype="legendre",
                is_periodic=self.dir_is_periodic[direction],
                quad_args={"quadperel": nb_pts_per_el},
            )
            quadrature.export_quadrature_rules()
            all_quadrules.append(quadrature)
        return all_quadrules

    def get_data_on_surface(
        self,
        indices: List[int],
        boundary_dir: ParametricDirection,
        p_degrees: Dict[ParametricDirection, int] = {},
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, List[StandardGauss]]:
        """
        Computes geometric data on a boundary surface for projection purposes.
        This includes:
        Physical coordinates of quadrature points on the boundary
        Jacobian of the transformation at quadrature points
        NURBS weights at quadrature points (if applicable)
        Quadrature rules used for integration on the boundary

        Args:
            indices (List[int]):
                List of control point indices that lie on the selected boundary.
            boundary_dir (ParametricDirection):
                Parametric direction on boundary.
            p_degrees (Dict[ParametricDirection, int]):
                Optional dict specifying the polynomial degrees in each parametric direction for the master patch.
                This is used to determine the appropriate quadrature rules for projection.
        """
        # 1. Get control points and weights on boundary
        ctrlpts = self.ctrlpts[:, indices]
        nurbs_weights = (
            self.nurbs_weights[indices] if self.nurbs_weights.size > 0 else np.array([])
        )

        # 2. Get quadrature on boundary
        p_degrees = {d: 0 for d in ParametricDirection} if not p_degrees else p_degrees
        quadrules = self._get_quadrature_for_boundary(
            boundary_dir, p_master_degrees=p_degrees
        )

        # 3. Mapping: jacobian and determinant
        jacobien, _, _, quadpts_phys = SinglePatch.eval_transformation(
            quadrules, ctrlpts, nurbs_weights=nurbs_weights
        )
        return (
            quadpts_phys,
            jacobien,
            nurbs_weights,
            quadrules,
        )

    def generate(self):
        "Initialize the patch by computing geometric transformations at quadrature points."
        det_jac, _, inv_jac, qp_phy = SinglePatch.eval_transformation(
            self.quadrule_list, self.ctrlpts, nurbs_weights=self.nurbs_weights
        )
        ndim_geom = MathMixin.detect_true_dimension(self.ctrlpts)
        self.ctrlpts = self.ctrlpts[:ndim_geom, ...]
        self.inv_jac = inv_jac[:ndim_geom, :ndim_geom, ...]
        self.qp_phy = qp_phy[:ndim_geom, ...]
        self.det_jac = det_jac
        logger.info(repr(self))
        return self

    def __repr__(self) -> str:
        message = f""""
            GEOMETRY:
            dimensionality: {self.ndim}
            polynomial degrees: {self.degree.tolist()}
            nb of control points: {self.nbctrlpts.tolist()}
            total nb of c.p. : {self.nbctrlpts_total}
            nb of quadrature points: {self.nbqp.tolist()}
            total nb of q.p. : {self.nbqp_total}
        """
        return message
