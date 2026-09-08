from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.base.cls import BasePreconditioner
from yeti_iga.pymfiga.common.numerics.solvers.nystrom_preconditioner import (
    NystromSchurPreconditioner,
)
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from typing import Optional, Literal
from scipy import sparse as sp
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.FASTDIAG")


class LagrangeFastDiagonalization:

    _VALID_TYPES = ["standard", "augmented"]
    _SMW_MF_SS = _SMW_PRECOND = _linear_solver = None

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

    @property
    def linear_solver(self):
        if self._linear_solver is None:
            self._linear_solver = LinearSolver(
                tolerance=Constants.TINY,
                maxiters=100,  # Find a better value
                linear_type="cg",
                verbose=False,
            )
            # # TODO: Test if GCRODR reduces the number of iterations
            # self._linear_solver = GCRODR(
            #     tolerance=Constants.TINY,
            #     maxiters=100,
            #     verbose=False,
            # )
        return self._linear_solver

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
        if np.isscalar(lagrange_penalty):
            self._lagrange_penalty = float(np.real(lagrange_penalty))
            if self.lagrange_type == "augmented":
                self._SMW_MF_SS = None
                self._SMW_PRECOND = None
                self._update_augmented()

    def _update_augmented(self):
        "Sherman-Morrison-Woodbury parameters"
        if self._SMW_MF_SS is not None and self._SMW_PRECOND is not None:
            return

        start = time()
        CC = self.constraint_matrix
        penalty = self.lagrange_penalty
        inveigvals = self.preconditioner.global_space_eigenvalues_inverse
        U = lambda x: self.preconditioner.matvec_space_eigenvectors(
            x, is_transpose=False
        )
        UT = lambda x: self.preconditioner.matvec_space_eigenvectors(
            x, is_transpose=True
        )

        def mf_SS(x: np.ndarray) -> np.ndarray:
            r = CC.T @ x
            s = self.preconditioner.apply_spatial_preconditioner(r)
            t = CC @ s
            return t + x / penalty

        nystrom = NystromSchurPreconditioner(
            C=CC,
            U=U,
            UT=UT,
            Qinv=inveigvals,
            mu=1 / penalty,
            rank=int(np.ceil(np.sqrt(CC.shape[0]))),
        )
        nystrom.build()

        self._SMW_MF_SS = lambda x: mf_SS(x)
        self._SMW_PRECOND = lambda x: nystrom(x)
        logger.info(repr(nystrom))
        logger.info(f"Update fast-diagonalization in {time() - start:.2e} seconds")

    def _apply_standard(self, array_in: np.ndarray) -> np.ndarray:
        return self.preconditioner.apply_spatial_preconditioner(array_in)

    def _apply_augmented(self, array_in: np.ndarray) -> np.ndarray:
        "Sherman-Morrison-Woodbury preconditioner"
        start = time()
        CC = self.constraint_matrix
        SS = self._SMW_MF_SS
        PRECOND = self._SMW_PRECOND
        assert callable(SS) and callable(PRECOND)

        zz = self._apply_standard(array_in)
        ss = CC @ zz

        # TODO: find a good preconditioner !
        # For instance Jacobi preconditioner is not worth it
        # and even Nystrom lower-rank preconditioner is not reducing
        # the number of iterations as expected
        # NOTE: Eventually, we can include a GCRODR solver
        # to reuse Ritz eigenvectors and save some time

        # NOTE: the convergence of the problem depends on the accuracy
        # of the linear solver: with a relaxed tolerance it does not converge
        ww = self.linear_solver.solve(SS, ss, PRECOND)["sol"]
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
