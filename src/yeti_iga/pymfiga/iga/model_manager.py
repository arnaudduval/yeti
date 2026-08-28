from yeti_iga.pymfiga.common.base.enum import TargetPriority, BoundarySide, ParametricDirection
from yeti_iga.pymfiga.common.base.cls import (
    MortarInterface,
    BaseMultiModel,
    BaseSingleModel,
)
from typing import Dict, Union, List, Tuple
from abc import ABC, abstractmethod
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MULTI_MODEL")


class ModelManager(BaseMultiModel, ABC):
    def __init__(
        self,
        patch_models_dict: Dict[
            Union[int, str],
            BaseSingleModel,
        ],
    ):
        logger.info("Initialize multi-model")
        self._patch_id_model_dict = patch_models_dict

        self._patch_ids: List[Union[str, int]] = []
        self._patch_models: List[BaseSingleModel] = []
        self._patch_offsets: Dict[Union[str, int], int] = {}
        self._interfaces: List[MortarInterface] = []
        self._slices: List[slice] = []
        self._size_of_arrays: int = 0
        self._identify_models()

        self._free_nodes: List[int] = []
        self._constraint_nodes: List[int] = []
        self._add_free_constraint_nodes()

    @property
    def patch_id_model_dict(self):
        return self._patch_id_model_dict

    @property
    def patch_ids(self):
        return self._patch_ids

    @property
    def patch_models(self):
        return self._patch_models

    @property
    def patch_offsets(self):
        return self._patch_offsets

    @property
    def interfaces(self):
        return self._interfaces

    @property
    def slices(self):
        return self._slices

    @property
    def constraint_nodes(self):
        return self._constraint_nodes

    @property
    def free_nodes(self):
        return self._free_nodes

    @property
    def num_nodes_patch(self):
        # NOTE: This only works if the number of DOFs is the same for all nodes
        return {
            pid: m.get_size_of_arrays() // m.nbvars
            for pid, m in self.patch_id_model_dict.items()
        }

    def _identify_models(self):

        patch_ids: List[Union[str, int]] = []
        patch_models: List[BaseSingleModel] = []
        size_arrays_list: List[int] = [0]
        interfaces: list[MortarInterface] = []

        nbvars = 0
        for pid, md in self.patch_id_model_dict.items():
            # NOTE: we assume that nbvars for all models is the same
            nbvars = max(nbvars, md.nbvars)
            patch_ids.append(pid)
            patch_models.append(md)
            output = md.boundary.select_nodes_for_gluing()[0]
            interfaces.extend(self._create_interfaces(pid, output))
            size_arrays_list.append(md.get_size_of_arrays())

        offsets = np.cumsum(size_arrays_list)
        slices = [slice(offsets[i], offsets[i + 1]) for i in range(len(offsets) - 1)]

        patch_offsets: Dict[Union[str, int], int] = {
            pid: offset for pid, offset in zip(patch_ids, offsets)
        }

        self._patch_ids = patch_ids
        self._patch_models = patch_models
        self._patch_offsets = patch_offsets
        self._interfaces = interfaces
        self._slices = slices
        self._size_of_arrays = offsets[-1]

    def _create_interfaces(
        self,
        master_id: Union[int, str],
        master_glued_data: Dict[
            Union[int, str],
            Tuple[TargetPriority, Tuple[ParametricDirection, BoundarySide, List]],
        ],
    ) -> List[MortarInterface]:
        interfaces = []
        for slave_id, (slave_priority, master_side) in master_glued_data.items():
            if slave_priority == TargetPriority.MASTER:
                # NOTE: slave can not be master
                continue
            slave = self.patch_id_model_dict[slave_id]
            slave_glued_data = slave.boundary.select_nodes_for_gluing()[0]
            # NOTE: by definition, slave should share with current master
            assert master_id in slave_glued_data.keys()
            master_priority, slave_side = slave_glued_data[master_id]
            assert master_priority == TargetPriority.MASTER
            interfaces.append(
                MortarInterface(
                    master_patch=master_id,
                    master_side=master_side,
                    slave_patch=slave_id,
                    slave_side=slave_side,
                )
            )
        return interfaces

    def _add_free_constraint_nodes(self):
        free_nodes: List[int] = []
        constraint_nodes: List[int] = []
        pid_boundary_list = [
            (id, self.patch_models[k].boundary) for k, id in enumerate(self.patch_ids)
        ]
        for pid, boundary in pid_boundary_list:
            offset_global = self._patch_offsets[pid]
            output = boundary.select_nodes_for_solving()
            tmp = offset_global + np.asarray(output[0])
            free_nodes.extend(tmp.tolist())
            tmp = offset_global + np.asarray(output[1])
            constraint_nodes.extend(tmp.tolist())

        self._free_nodes = free_nodes
        self._constraint_nodes = constraint_nodes

    def recover_models(self) -> List[BaseSingleModel]:
        return self.patch_models

    def get_size_of_arrays(self):
        return self._size_of_arrays

    def get_free_and_constraint_nodes(self):
        return self.free_nodes, self.constraint_nodes

    def cut(self, array_in: np.ndarray) -> Dict[Union[int, str], np.ndarray]:
        return {pid: array_in[s] for pid, s in zip(self.patch_ids, self.slices)}

    def glue(self, id_array_in: Dict[Union[int, str], np.ndarray]) -> np.ndarray:
        if len(id_array_in) == 0:
            return np.zeros(self.get_size_of_arrays())

        dtype = id_array_in[self.patch_ids[0]].dtype
        out = np.zeros(self.get_size_of_arrays(), dtype=dtype)
        for pid, s in zip(self.patch_ids, self.slices):
            block = id_array_in.get(pid)
            if block is None:
                out[s] = 0.0
            else:
                out[s] = block
        return out

    @abstractmethod
    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        pass

    @abstractmethod
    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        pass
