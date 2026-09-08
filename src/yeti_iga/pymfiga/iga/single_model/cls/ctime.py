from yeti_iga.pymfiga.common.base.math import kron_nonzero_indices
from yeti_iga.pymfiga.common.numerics.operations import BsplineOperations, NurbsOperations
from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from yeti_iga.pymfiga.common.material import Material
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from .core import SingleModel
from typing import Callable, List, Dict
from abc import ABC, abstractmethod
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MODEL")


class SingleSpaceTimeModel(SingleModel, ABC):

    _sptm_free_nodes: List[int] = []
    _sptm_constraint_nodes: List[int] = []
    _sptm_nurbs_weights: np.ndarray = np.array([])
    _sptm_quadrature_list: List[IGAQuadratureRule] = []

    def __init__(
        self,
        material: Material,
        space_patch: SinglePatch,
        time_patch: SinglePatch,
        boundary: BoundaryCondition,
    ):
        if time_patch.ndim != 1:
            raise ValueError("Time is only one dimension")
        super().__init__(material, space_patch, boundary)
        self._set_time(time_patch)

        # Propagate quadrature list in time
        self._sptm_quadrature_list: List[IGAQuadratureRule] = []
        self._sptm_free_nodes: List[int] = []
        self._sptm_constraint_nodes: List[int] = []
        self._sptm_nurbs_weights: np.ndarray = np.array([])
        self._propagate_quadrule_in_time()
        self._propagate_free_and_constraint_nodes_in_time()
        self._propagate_nurbs_weights_in_time()

    @property
    def time(self):
        return self._time

    @property
    def sptm_constraint_nodes(self):
        return self._sptm_constraint_nodes

    @property
    def sptm_free_nodes(self):
        return self._sptm_free_nodes

    @property
    def sptm_nurbs_weights(self):
        return self._sptm_nurbs_weights

    @property
    def sptm_quadrature_list(self):
        return self._sptm_quadrature_list

    def _set_time(self, value):
        assert isinstance(value, SinglePatch)
        self._time = value

    def _propagate_quadrule_in_time(self):
        self._sptm_quadrature_list = self.part.quadrule_list + self.time.quadrule_list

    def _propagate_free_and_constraint_nodes_in_time(self):
        nbctrlpts_sp = self.nbvars * self.part.nbctrlpts_total
        nbctrlpts_tm = self.time.nbctrlpts_total
        nnz_sptm_list = [nbctrlpts_tm, nbctrlpts_sp]

        indices_list = [
            np.arange(1, nbctrlpts_tm, dtype=int).tolist(),
            self.free_nodes,
        ]
        global_indices = kron_nonzero_indices(indices_list, nnz_sptm_list)
        self._sptm_free_nodes = global_indices

        indices_list = [
            np.arange(1, nbctrlpts_tm, dtype=int).tolist(),
            self.constraint_nodes,
        ]
        global_indices = np.arange(nbctrlpts_sp).tolist() + kron_nonzero_indices(
            indices_list, nnz_sptm_list
        )
        self._sptm_constraint_nodes = global_indices

    def _propagate_nurbs_weights_in_time(self):
        space_weights = self.part.nurbs_weights
        time_weights = self.time.nurbs_weights
        if space_weights.size > 0:
            if time_weights.size > 0:
                nurbs_weights = np.kron(
                    self.time.nurbs_weights, self.part.nurbs_weights
                )
            else:
                nurbs_weights = np.kron(
                    np.ones(self.time.nbctrlpts_total), self.part.nurbs_weights
                )
        else:
            if time_weights.size > 0:
                self.operator_engine = NurbsOperations
                nurbs_weights = np.kron(
                    self.time.nurbs_weights, np.ones(self.part.nbctrlpts_total)
                )
            else:
                self.operator_engine = BsplineOperations
                nurbs_weights = np.array([])

        self._sptm_nurbs_weights = nurbs_weights

    def _assemble_scalar_volume_force(
        self, fun: Callable, **external_args
    ) -> np.ndarray:
        start = time()
        fun_evaluated = np.atleast_2d(fun(external_args))
        assert np.ndim(fun_evaluated) == 2
        geometry_prop = np.kron(self.time.det_jac, self.part.det_jac)
        prop: np.ndarray = np.einsum("i,ji->ji", geometry_prop, fun_evaluated)
        array_out = self.operator_engine.assemble_scalar_u_force(
            self.sptm_quadrature_list, prop, nurbs_weights=self.sptm_nurbs_weights
        )
        logger.info(f"Computing volume force in {time() - start:.2e} seconds")
        return array_out

    def verify_fun_args(self, args: dict):
        # Overwrite SingleModel function
        assert isinstance(
            args, dict
        ), "allowed extra arguments should be in dictionnary"
        args.update(
            {
                "position": self.part.qp_phy,
                "shape_quadpts": (self.part.nbqp_total * self.time.nbqp_total,),
                "time": np.ravel(self.time.qp_phy),
            }
        )

    def export_force(self, array_in: dict) -> np.ndarray:
        # Overwrite SingleModel function
        boundary = self.boundary
        nnz_time = self.time.nbctrlpts_total
        dof_idx_list = {dof: idx for dof, (idx, _) in boundary.dofs_index.items()}
        array_out = [
            np.zeros(boundary.nbofdofs[dof] * nnz_time) for dof in dof_idx_list.keys()
        ]

        for dof, array in array_in.items():
            if dof in dof_idx_list.keys():
                idx = dof_idx_list[dof]
                array_out[idx] += np.ravel(array)
            elif dof is DOF.ALL:
                for idx in range(len(dof_idx_list)):
                    array_out[idx] += np.ravel(array[idx])
        return np.hstack(array_out)

    def assemble_volumetric_force(
        self, fun_info: Dict[DOF, Callable], **external_args
    ) -> np.ndarray:
        # Overwrite SingleModel function
        assert isinstance(fun_info, dict)
        assert all(
            callable(fun) for _, fun in fun_info.items()
        ), "Insert a list of functions"
        self.verify_fun_args(external_args)
        array_out = {}
        for dof, fun in fun_info.items():
            output = self._assemble_scalar_volume_force(fun, **external_args)
            array_out.update({dof: output})
        return self.export_force(array_out)

    def clear_bcs(self, array_in: np.ndarray):
        # Overwrite SingleModel function
        assert isinstance(self.sptm_constraint_nodes, list)
        array_in[self.sptm_constraint_nodes] = 0.0

    def get_free_and_constraint_nodes(self):
        # Overwrite SingleModel function
        return self.sptm_free_nodes, self.sptm_constraint_nodes

    @abstractmethod
    def compute_residual(self, array_in, **kwargs):
        raise NotImplementedError()

    @abstractmethod
    def solve_linearized_system(self, array_in, **kwargs):
        raise NotImplementedError()
