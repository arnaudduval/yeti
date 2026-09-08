from yeti_iga.pymfiga.common.base.enum import DOF, ParametricDirection, BoundarySide
from yeti_iga.pymfiga.common.material import Material
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from .core import SingleModel
from typing import Callable, Tuple, Dict
from abc import ABC, abstractmethod
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MODEL")


class SingleSpatialModel(SingleModel, ABC):
    def __init__(
        self,
        material: Material,
        patch: SinglePatch,
        boundary: BoundaryCondition,
    ):
        super().__init__(material, patch, boundary)

    def verify_fun_args(self, args: dict):
        assert isinstance(
            args, dict
        ), "allowed extra arguments should be in dictionnary"
        args.update(
            {
                "position": self.part.qp_phy,
                "shape_quadpts": (self.part.nbqp_total,),
            }
        )

    def _assemble_scalar_volume_force(
        self, fun: Callable, **external_args
    ) -> np.ndarray:
        start = time()
        fun_evaluated = np.atleast_2d(fun(external_args))
        assert np.ndim(fun_evaluated) == 2
        prop: np.ndarray = np.einsum("i,ji->ji", self.part.det_jac, fun_evaluated)
        array_out = self.operator_engine.assemble_scalar_u_force(
            self.part.quadrule_list, prop, nurbs_weights=self.part.nurbs_weights
        )
        logger.info(f"Computing volume force in {time() - start:.2e} seconds")
        return array_out

    def _assemble_scalar_surface_force(
        self, fun: Callable, info: dict, **external_args
    ) -> np.ndarray:
        start = time()
        assert self.ndim > 1, "Method only valid for multivariate geometries"
        assert callable(fun) and isinstance(info, dict)

        # Boundary definition
        idx_direction: ParametricDirection = info["direction"]
        idx_side: BoundarySide = info["face"]
        assert isinstance(idx_direction, ParametricDirection)
        assert isinstance(idx_side, BoundarySide)

        # Get data on boundary
        sltd_indices = self.boundary.boundary_dofs[(idx_direction, idx_side)]
        output = self.part.get_data_on_surface(sltd_indices, idx_direction)
        (
            sltd_quadpts_phys,
            sltd_dsurf,
            sltd_nurbs_weights,
            sltd_quadrules,
        ) = output

        # Evaluate the function on quadrature points
        args = {**external_args, "position": sltd_quadpts_phys}
        # fun_evaluated: (N_components, N_quadpts)
        fun_evaluated = np.atleast_2d(fun(args))
        assert fun_evaluated.ndim == 2

        # Prepare data for integration
        # Prop: (N_gauss, N_components)
        prop: np.ndarray = np.einsum("i,ji->ji", sltd_dsurf, fun_evaluated)

        # Assemble the vector
        array_out = np.zeros((prop.shape[0], self.part.nbctrlpts_total))
        array_out[:, sltd_indices] = self.operator_engine.assemble_scalar_u_force(
            sltd_quadrules, prop, nurbs_weights=sltd_nurbs_weights
        )

        logger.info(f"Computing surface force in {time() - start:.2e} seconds")
        return array_out

    def assemble_surface_force(
        self, fun_info: Dict[DOF, Tuple[dict, Callable]], **external_args
    ) -> np.ndarray:
        assert isinstance(fun_info, dict)
        assert all(
            isinstance(dof, DOF) and isinstance(cond, dict) and callable(fun)
            for dof, (cond, fun) in fun_info.items()
        )
        self.verify_fun_args(external_args)
        array_out = {}
        for dof, (info, fun) in fun_info.items():
            output = self._assemble_scalar_surface_force(fun, info, **external_args)
            array_out.update({dof: output})
        return super().export_force(array_out)

    def assemble_volumetric_force(
        self, fun_info: Dict[DOF, Callable], **external_args
    ) -> np.ndarray:
        assert isinstance(fun_info, dict)
        assert all(
            callable(fun) for _, fun in fun_info.items()
        ), "Insert a list of functions"
        self.verify_fun_args(external_args)
        array_out = {}
        for dof, fun in fun_info.items():
            output = self._assemble_scalar_volume_force(fun, **external_args)
            array_out.update({dof: output})
        return super().export_force(array_out)

    @abstractmethod
    def compute_residual(self, array_in, **kwargs):
        raise NotImplementedError()

    @abstractmethod
    def solve_linearized_system(self, array_in, **kwargs):
        raise NotImplementedError()
