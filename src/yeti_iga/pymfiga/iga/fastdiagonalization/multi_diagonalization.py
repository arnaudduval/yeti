from yeti_iga.pymfiga.common.base.cls import BasePreconditioner
from yeti_iga.pymfiga.iga.model_manager import ModelManager
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.FASTDIAG")


class MultiFastDiagonalization(BasePreconditioner):
    def __init__(
        self,
        multipatch: ModelManager,
    ):
        self._multipatch = multipatch
        self._global_eigenvectors = None

    @property
    def global_space_eigenvalues_inverse(self) -> np.ndarray:
        local_eig_list = []
        for model in self._multipatch.patch_models:
            fd = model.preconditioner
            assert isinstance(fd, BasePreconditioner)
            local_eig = fd.global_space_eigenvalues_inverse
            local_eig_list.append(local_eig)
        eigenvalues = np.hstack(local_eig_list)
        return eigenvalues

    def matvec_space_eigenvectors(self, array_in: np.ndarray, is_transpose: bool):
        array_in_cutted = self._multipatch.cut(array_in)
        array_out_dict = {}
        for pid, model in zip(
            self._multipatch.patch_ids, self._multipatch.patch_models
        ):
            fd = model.preconditioner
            assert isinstance(fd, BasePreconditioner)
            array_out_cutted = fd.matvec_space_eigenvectors(
                array_in=array_in_cutted[pid],
                is_transpose=is_transpose,
            )
            array_out_dict.update({pid: array_out_cutted})
        array_out = self._multipatch.glue(array_out_dict)
        return array_out

    def apply_spatial_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        start = time()
        array_in_cutted = self._multipatch.cut(array_in)
        array_out_dict = {}
        for pid, model in zip(
            self._multipatch.patch_ids, self._multipatch.patch_models
        ):
            fd = model.preconditioner
            assert isinstance(fd, BasePreconditioner)
            array_out_cutted = fd.apply_spatial_preconditioner(
                array_in=array_in_cutted[pid]
            )
            array_out_dict.update({pid: array_out_cutted})
        array_out = self._multipatch.glue(array_out_dict)
        logger.debug(
            f"Multi-patch fast-diagonalization in {time() - start:.2e} seconds"
        )
        return array_out

    def apply_spacetime_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        raise NotImplementedError()
