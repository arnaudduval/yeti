from .enum import TargetPriority, BoundarySide, ParametricDirection, DOF
from typing import Sequence, Tuple, Dict, Set, List, Optional, Literal, Union, Any
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np


class BaseBoundaryCondition(ABC):
    @property
    @abstractmethod
    def dofs_index(self) -> dict:
        pass

    @property
    @abstractmethod
    def nbofdofs(self) -> dict:
        "Number of Degree of Freedom (DOF)"
        pass

    @abstractmethod
    def select_nodes_for_contact(self) -> Tuple[Set[int], Set[tuple]]:
        pass

    @abstractmethod
    def select_nodes_for_solving(
        self,
    ) -> Tuple[List[int], List[int]]:
        pass

    @abstractmethod
    def select_nodes_for_gluing(
        self,
    ) -> Tuple[
        Dict[
            Union[int, str],
            Tuple[TargetPriority, Tuple[ParametricDirection, BoundarySide, List]],
        ],
        List,
    ]:
        pass

    @abstractmethod
    def add_constraint(self, constraint_info, constraint_type) -> None:
        pass


class BaseMaterial(ABC):
    @property
    @abstractmethod
    def hasnonlinearmass(self) -> bool:
        pass

    @property
    @abstractmethod
    def hasnonlinearstiffness(self) -> bool:
        pass


class BasePatch(ABC):

    ndim: int

    @abstractmethod
    def compute_global_mesh_parameter(self) -> Tuple[float, float]:
        pass


class BasePreconditioner(ABC):
    @property
    @abstractmethod
    def global_space_eigenvalues_inverse(self) -> np.ndarray:
        pass

    @abstractmethod
    def matvec_space_eigenvectors(
        self, array_in: np.ndarray, is_transpose: bool
    ) -> np.ndarray:
        pass

    @abstractmethod
    def apply_spatial_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def apply_spacetime_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        pass


class BaseSingleModel(ABC):

    TypeProblemToBeApplied: Literal["FEA", "IGA"]

    @property
    @abstractmethod
    def update_manager(self) -> Any:
        pass

    @property
    @abstractmethod
    def part(self) -> BasePatch:
        pass

    @property
    @abstractmethod
    def boundary(self) -> BaseBoundaryCondition:
        pass

    @property
    @abstractmethod
    def material(self) -> BaseMaterial:
        pass

    @property
    def time(self) -> Optional[BasePatch]:
        pass

    @property
    def preconditioner(self) -> Optional[BasePreconditioner]:
        pass

    @property
    def nbvars(self):
        "Number of DoFs per node"
        return len(self.boundary.dofs_index)  # DoF per node

    @abstractmethod
    def set_update_manager(self, manager):
        pass

    @abstractmethod
    def get_size_of_arrays(self) -> int:
        "Get number of DoFs in space"
        pass

    @abstractmethod
    def get_free_and_constraint_nodes(self) -> Tuple[list, list]:
        pass

    @abstractmethod
    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, Dict]:
        pass

    @abstractmethod
    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        pass

    def export_force(self, array_in: dict) -> np.ndarray:
        assert self.boundary is not None
        boundary = self.boundary
        dof_idx_list = {dof: idx for dof, (idx, _) in boundary.dofs_index.items()}
        array_out = [np.zeros(boundary.nbofdofs[dof]) for dof in dof_idx_list.keys()]

        for dof, array in array_in.items():
            if dof in dof_idx_list.keys():
                idx = dof_idx_list[dof]
                array_out[idx] += np.ravel(array)
            elif dof is DOF.ALL:
                for idx in range(len(dof_idx_list)):
                    array_out[idx] += np.ravel(array[idx])
        return np.hstack(array_out)


class BaseMultiModel(ABC):

    TypeProblemToBeApplied: Literal["IGA"] = "IGA"

    @property
    @abstractmethod
    def update_manager(self) -> Any:
        pass

    @property
    def nbvars(self):
        "Number of DoFs per node"
        return np.max([m.nbvars for m in self.recover_models()])  # DoF per node

    @abstractmethod
    def set_update_manager(self, manager):
        pass

    @abstractmethod
    def recover_models(self) -> Sequence[BaseSingleModel]:
        pass

    @abstractmethod
    def get_size_of_arrays(self) -> int:
        "Get number of DoFs in space"
        pass

    @abstractmethod
    def get_free_and_constraint_nodes(self) -> Tuple[list, list]:
        pass

    @abstractmethod
    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        pass

    @abstractmethod
    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        pass


@dataclass
class MortarInterface:

    master_patch: Union[int, str]
    master_side: Tuple[ParametricDirection, BoundarySide, List[int]]
    slave_patch: Union[int, str]
    slave_side: Tuple[ParametricDirection, BoundarySide, List[int]]

    def __post_init__(self):
        if not isinstance(self.master_side, tuple):
            assert ValueError("master_side should be a dictionnary")
        if not isinstance(self.master_side[0], ParametricDirection) or not isinstance(
            self.master_side[1], BoundarySide
        ):
            assert ValueError("master_side should containt 'direction' and 'face'")
        if not isinstance(self.slave_side, tuple):
            assert ValueError("slave_side should be a dictionnary")
        if not isinstance(self.slave_side[0], ParametricDirection) or not isinstance(
            self.slave_side[1], BoundarySide
        ):
            assert ValueError("slave_side should containt 'direction' and 'face'")
