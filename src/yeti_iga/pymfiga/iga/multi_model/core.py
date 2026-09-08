from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.numerics.solvers import LagrangeSolver, UpdateManager
from yeti_iga.pymfiga.iga.single_model.cls.cspace import SingleSpatialModel
from yeti_iga.pymfiga.iga.fastdiagonalization import MultiFD, LagrangeFD
from yeti_iga.pymfiga.iga.model_manager import ModelManager
from .dual_mortar import MortarCoupling
from typing import Dict, Union, Literal
from abc import ABC, abstractmethod
import numpy as np


class MultiModel(ModelManager, ABC):
    _lagrange_type: Literal["augmented", "standard"] = "augmented"
    _preconditioner = _lagrange_solver = None

    def __init__(
        self,
        patch_models_dict: Dict[
            Union[int, str],
            SingleSpatialModel,
        ],
        lagrange_type: Literal["augmented", "standard"] = "augmented",
    ):
        super().__init__(patch_models_dict)
        self._lagrange_type = lagrange_type
        self._prepare_lagrange_solver()

    @property
    def lagrange_type(self):
        return self._lagrange_type

    @property
    def constraint_matrix(self):
        return self._constraint_matrix

    @property
    def constraint_vector(self):
        return self._constraint_vector

    @property
    def update_manager(self):
        if not isinstance(self._update_manager, UpdateManager):
            self._update_manager = UpdateManager()
        return self._update_manager

    @property
    def preconditioner(self):
        if isinstance(self._preconditioner, LagrangeFD):
            return self._preconditioner

        lagfd = LagrangeFD(
            preconditioner=MultiFD(self),
            constraint_matrix=self.constraint_matrix,
            lagrange_type=self.lagrange_type,
        )
        self._preconditioner = lagfd
        return lagfd

    @property
    def lagrange_solver(self):
        if isinstance(self._lagrange_solver, LagrangeSolver):
            return self._lagrange_solver

        # TODO: how modify values (tolerance and maxiters) for Lagrange solver ?
        solver = LagrangeSolver(
            constraint_matrix=self.constraint_matrix,
            constraint_vector=self.constraint_vector,
            tolerance=Constants.TINY,
            maxiters=100,
            cleandod=self.constraint_nodes,
            lagrange_type=self.lagrange_type,
        )
        self._lagrange_solver = solver
        return solver

    def _prepare_lagrange_solver(self):
        coupling = MortarCoupling(self)
        mat, vec = coupling.build_constraint_matrix_vector()
        self._constraint_matrix = mat
        self._constraint_vector = vec

    def set_update_manager(self, manager):
        if isinstance(manager, UpdateManager):
            self._update_manager = manager
            for model in self._patch_models:
                model.set_update_manager(manager)

    def clear_bcs(self, array_in: np.ndarray):
        array_in[self.constraint_nodes] = 0.0

    @abstractmethod
    def compute_residual(self, array_in, **kwargs):
        raise NotImplementedError()

    @abstractmethod
    def solve_linearized_system(self, array_in, **kwargs):
        raise NotImplementedError()
