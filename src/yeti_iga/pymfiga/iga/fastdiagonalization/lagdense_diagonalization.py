from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.base.cls import BasePreconditioner
from typing import Optional, Literal
from scipy import linalg as sclin
from scipy import sparse as sp
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.FASTDIAG")


class LagrangeFastDiagonalization:

    _VALID_TYPES = ["standard", "augmented"]
    _SMW_SS_FACT = None

    def __init__(
        self,
        preconditioner: BasePreconditioner,
        constraint_matrix: sp.csr_array,
        lagrange_type: Literal["standard", "augmented"] = "augmented",
    ):
        # Public variables
        self._lagrange_type = self._validate_lagrange_type(lagrange_type)
        self._constraint_matrix = sp.csr_array(constraint_matrix)
        self._lagrange_penalty = Constants.PENALTY
        assert isinstance(preconditioner, BasePreconditioner)
        self._preconditioner = preconditioner

    @property
    def lagrange_type(self):
        return self._lagrange_type

    @property
    def constraint_matrix(self):
        return self._constraint_matrix

    @property
    def lagrange_penalty(self):
        return self._lagrange_penalty

    @property
    def preconditioner(self):
        return self._preconditioner

    def _validate_lagrange_type(self, value: str):
        if not isinstance(value, str):
            raise TypeError(f"Lagrange type must be string, got {type(value)}")

        if value not in self._VALID_TYPES:
            raise ValueError(f"""
                Invalid tolerance type {value}.
                f"Must be one of {self._VALID_TYPES}
                """)
        return value

    def update(self, lagrange_penalty: Optional[float] = None):
        # TODO: implement update only if the new penalty differs
        # in more than 10% with respect the old value
        if np.isscalar(lagrange_penalty):
            self._lagrange_penalty = float(np.real(lagrange_penalty))
            if self.lagrange_type == "augmented":
                self._SMW_SS_FACT = None
                self._update_augmented()

    def _update_augmented(self):
        "Sherman-Morrison-Woodbury parameters"
        if self._SMW_SS_FACT is not None:
            return

        start = time()
        penalty = self.lagrange_penalty
        CC = self.constraint_matrix
        nr = CC.shape[0]

        SS = np.zeros((nr, nr))
        for i in range(nr):
            cT = CC[i, :].todense()
            invPcT = self.preconditioner.apply_spatial_preconditioner(cT)
            SS[:, i] = CC @ invPcT
        SS += np.eye(nr) / penalty

        self._SMW_SS_FACT = sclin.lu_factor(SS)
        logger.info(f"Update fast-diagonalization in {time() - start:.2e} seconds")

    def _apply_standard(self, array_in: np.ndarray) -> np.ndarray:
        return self.preconditioner.apply_spatial_preconditioner(array_in)

    def _apply_augmented(self, array_in: np.ndarray) -> np.ndarray:
        "Sherman-Morrison-Woodbury preconditioner"
        start = time()
        CC = self.constraint_matrix
        SS_FACT = self._SMW_SS_FACT
        assert SS_FACT is not None

        zz = self._apply_standard(array_in)
        ss = CC @ zz
        ww = sclin.lu_solve(SS_FACT, ss)
        tt = CC.T @ ww
        uu = zz - self._apply_standard(tt)
        logger.debug(f"Lagrange fast-diagonalization in {time() - start:.2e} seconds")
        return uu

    def apply_spatial_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        fun = {
            "augmented": self._apply_augmented,
            "standard": self._apply_standard,
        }[self.lagrange_type]
        return fun(array_in)
