from yeti_iga.pymfiga.common.base.enum import Constants
from .core import SolverArgs
from .linear_solver import LinearSolver
from typing import Union, Callable, List, Optional, Literal
from scipy.sparse import linalg as scsplin
from dataclasses import dataclass
from scipy import linalg as sclin
from scipy import sparse as sp
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.SOLVER")


@dataclass
class LagrangeArgs(SolverArgs):
    lagrange_penalty: float

    def __post_init__(self):
        """Validate arguments after initialization."""
        if self.lagrange_penalty <= 0:
            raise ValueError(f"Penalty must be > 0, got {self.lagrange_penalty}")


class LagrangeSolver:
    """
    Solve K u = f with constraints C u = g using Lagrange multipliers.
    Here,  K could be semi-definite and C should be full row rank
    Supports:
        - standard: Standard Lagrange method
        - augmented: Augmented Lagrangian method
    """

    _VALID_TYPES = ["augmented", "standard"]
    _cleandod = _kernel = _linear_solver = None

    def __init__(
        self,
        constraint_matrix: Union[np.ndarray, sp.csr_array],
        constraint_vector: np.ndarray,
        tolerance: float,
        maxiters: int,
        cleandod: Optional[List[int]] = None,
        lagrange_type: Literal["augmented", "standard"] = "augmented",
        verbose: bool = True,
    ):
        logger.debug("Initialize Lagrange solver")
        self._lagrange_type = self._validate_lagrange_type(lagrange_type)
        self._cleandod = cleandod
        self._config = LagrangeArgs(
            tolerance=tolerance,
            maxiters=maxiters,
            lagrange_penalty=Constants.PENALTY,
        )
        self._original_Cmat = constraint_matrix
        self._original_gvec = constraint_vector
        self._lag_mult = np.zeros_like(constraint_vector)
        self._verbose = verbose

        # Define C_mat, g_vec
        self._set_constraints(
            self._original_Cmat.copy(),
            self._original_gvec.copy(),
        )

    @property
    def lagrange_type(self) -> str:
        return self._lagrange_type

    @property
    def config(self) -> LagrangeArgs:
        return self._config

    @property
    def constraint_matrix(self) -> sp.csr_array:
        return self._C_mat

    @property
    def constraint_vector(self) -> np.ndarray:
        return self._g_vec

    @property
    def lagrange_multiplier(self) -> np.ndarray:
        return self._lag_mult.copy()

    @property
    def kernel(self):
        if self._kernel is None:
            self._kernel = LagrangeSolver._compute_kernel(self.constraint_matrix)
        return self._kernel

    @property
    def linear_solver(self):
        if not isinstance(self._linear_solver, LinearSolver):
            cleandod = self._cleandod if self.lagrange_type == "augmented" else None
            linear_solver = LinearSolver(
                tolerance=self.config.tolerance,
                maxiters=self.config.maxiters,
                cleandod=cleandod,
                linear_type="gmres",
                verbose=self._verbose,
            )
            self._linear_solver = linear_solver
        return self._linear_solver

    def update(
        self,
        tolerance: Optional[float] = None,
        maxiters: Optional[float] = None,
        lagrange_penalty: Optional[float] = None,
    ):
        """Update solver parameters and optionally the Dirichlet mask."""
        kwargs = self.config.__dict__.copy()
        if tolerance is not None:
            kwargs["tolerance"] = tolerance
        if maxiters is not None:
            kwargs["maxiters"] = maxiters
        if lagrange_penalty is not None:
            kwargs["lagrange_penalty"] = lagrange_penalty
        self._config = LagrangeArgs(**kwargs)

    def set_cleandod(self, cleandod: Optional[List[int]]):
        if cleandod is not None:
            assert isinstance(cleandod, list)
            self._cleandod = cleandod
            # Since we have update the cleandod,
            # linear solver is wrong for augmented lagrangian
            if self.lagrange_type == "augmented":
                self._linear_solver = None
            # kernel is wrong for standard lagrangian
            if self.lagrange_type == "standard":
                self._set_constraints(
                    self._original_Cmat.copy(), self._original_gvec.copy()
                )
                self._kernel = None

    def _apply_mask(self, x: np.ndarray):
        if self._cleandod is not None:
            x[self._cleandod] = 0.0

    def _validate_lagrange_type(self, value: str):
        if not isinstance(value, str):
            raise TypeError(f"Lagrange type must be string, got {type(value)}")

        if value not in self._VALID_TYPES:
            raise ValueError(
                f"Invalid tolerance type '{value}'. "
                f"Must be one of {self._VALID_TYPES}"
            )
        return value

    def _set_constraints(
        self,
        matrix: Union[np.ndarray, sp.csr_array],
        vector: np.ndarray,
    ):
        assert isinstance(matrix, (np.ndarray, sp.csr_array))
        assert isinstance(vector, np.ndarray) and vector.ndim == 1
        C_mat = matrix.copy()
        g_vec = vector.copy()
        if isinstance(self._cleandod, list) and len(self._cleandod) > 0:
            nr = len(self._cleandod)
            nc = matrix.shape[-1]
            mat = np.zeros((nr, nc))
            mat[:, self._cleandod] = np.eye(nr)
            vec = np.zeros(nr)
            C_mat = (
                np.vstack((matrix, mat))
                if isinstance(matrix, np.ndarray)
                else sp.vstack((matrix, mat))
            )
            g_vec = np.hstack((vector, vec))
        C_mat = sp.csr_array(C_mat)
        C_mat.eliminate_zeros()
        self._C_mat = sp.csr_array(C_mat)
        self._g_vec = g_vec
        self._lag_mult = np.zeros_like(g_vec)

    @staticmethod
    def _compute_kernel(C: Union[np.ndarray, sp.csr_array]):
        """
        Computes the nullspace (kernel) of the constraint matrix.

        Args:
            C (np.ndarray or sp.sparray): Constraint matrix.

        Returns:
            np.ndarray or sp.sparray: Nullspace basis.
        """
        # TODO: find kernel directly for sparse matrices.
        # There is a functionality in scipy but there are a lot of issues
        # when I was developing this class
        assert isinstance(C, np.ndarray) or sp.issparse(C)
        start = time()
        mat = C.toarray() if isinstance(C, sp.csr_array) else C
        ker = sclin.null_space(mat)
        logger.info(
            f"Null space of constraint matrix computed in {time() - start:.2e} seconds"
        )
        return LagrangeSolver._maybe_convert_sparse(ker)

    @staticmethod
    def _maybe_convert_sparse(
        Z: np.ndarray, density_threshold: float = 0.10, memory_saving: float = 0.8
    ) -> Union[np.ndarray, sp.csr_array]:
        """
        Converts a dense matrix to sparse if it is memory efficient.

        Args:
            Z (np.ndarray): Matrix to check.
            density_threshold (float): Density threshold for conversion.
            memory_saving (float): Minimum memory saving ratio.

        Returns:
            np.ndarray or sp.sparray: Converted matrix.
        """
        if sp.issparse(Z):
            return Z
        Z_csr = sp.csr_array(Z)
        Z_csr.eliminate_zeros()
        sparse_bytes = Z_csr.data.nbytes + Z_csr.indices.nbytes + Z_csr.indptr.nbytes
        dense_bytes = Z.nbytes
        density = Z_csr.nnz / Z.size
        if density < density_threshold or sparse_bytes < memory_saving * dense_bytes:
            logger.debug(
                f"Converting to sparse: density={100*density:.2f}%, dense={dense_bytes/1e6:.2f} MB, sparse={sparse_bytes/1e6:.2f} MB"
            )
            return Z_csr
        logger.debug(
            f"Keeping dense: density={100*density:.2f}%, dense={dense_bytes/1e6:.2f} MB, sparse={sparse_bytes/1e6:.2f} MB"
        )
        return Z

    def _solve(
        self,
        Afun: Callable,
        bvec: np.ndarray,
        Pfun: Optional[Callable] = None,
        **mf_args,
    ) -> np.ndarray:
        """
        Solves a linear system using the provided matrix-free operator.

        Args:
            Afun (callable): Linear operator function.
            bvec (np.ndarray): Right-hand side vector.
            Pfun (callable, optional): Preconditioner function.
            **mf_args: Additional arguments for the solver.

        Returns:
            np.ndarray: Solution vector.
        """

        linsolver = self.linear_solver
        linsolver.update(
            tolerance=self.config.tolerance,
            maxiters=self.config.maxiters,
        )
        return linsolver.solve(Afun, bvec, Pfun=Pfun, **mf_args)["sol"]

    @staticmethod
    def _lstsq(
        M: Union[np.ndarray, sp.csr_array],
        rhs: np.ndarray,
        is_transpose=False,
    ) -> np.ndarray:
        """
        Solves a least squares problem.

        Args:
            M (np.ndarray or sp.sparray): Matrix.
            rhs (np.ndarray): Right-hand side vector.
            is_transpose (bool): If True, uses M.T.

        Returns:
            np.ndarray: Solution vector.
        """
        assert isinstance(M, np.ndarray) or isinstance(M, sp.csr_array)
        if isinstance(M, np.ndarray):
            return np.linalg.lstsq(M.T if is_transpose else M, rhs, rcond=None)[0]
        return scsplin.lsqr(M.T if is_transpose else M, rhs)[0]

    def _compute_constraint_residual(self, x: np.ndarray) -> np.ndarray:
        "Computes g - C @ x"
        C = self.constraint_matrix
        g = self.constraint_vector
        return g - C @ x

    def _increment_standard(
        self,
        apply_T: Callable,
        apply_P: Optional[Callable],
        current_res: np.ndarray,
        current_sol: np.ndarray,
        **mf_args,
    ) -> np.ndarray:
        """
        The linear system to solve is:
        |T   C.T| |du  | |res_u  |
        |C    0 | |dlam|=|res_lam|
        Args:
            apply_T: it is the tangent matrix
            apply_P: it is a 'good' preconditioner for T
        """
        Z = self.kernel
        C = self.constraint_matrix

        res_lam = self._compute_constraint_residual(current_sol)
        delta_up = LagrangeSolver._lstsq(C, res_lam, is_transpose=False)
        self._apply_mask(delta_up)

        def red_matvec(x_red: np.ndarray, **args) -> np.ndarray:
            x = Z @ x_red
            x_mask = x.copy()
            self._apply_mask(x_mask)
            y = apply_T(x_mask, **args)
            self._apply_mask(y)
            return Z.T @ y

        def red_preconditioner(x_red: np.ndarray) -> np.ndarray:
            # NOTE: Z @ xred is equivalent to solve Z.T y = xred
            # since Z is orthogonal, ie, Z @ Z.T = Identity
            y = Z @ x_red
            self._apply_mask(y)
            w = apply_P(y) if callable(apply_P) else np.copy(y)
            self._apply_mask(w)
            # NOTE: Z.T @ w is equivalent to solve Z out = w
            # since Z is orthogonal, ie, Z.T @ Z = Identity
            return Z.T @ w

        rhs = current_res - apply_T(delta_up, **mf_args)
        self._apply_mask(rhs)
        rhs_red = Z.T @ rhs

        y = self._solve(red_matvec, rhs_red, red_preconditioner, **mf_args)
        delta_ug = Z @ y
        delta_u = delta_up + delta_ug
        self._apply_mask(delta_u)

        Tdug = apply_T(delta_ug, **mf_args)
        self._apply_mask(Tdug)
        self._lag_mult += LagrangeSolver._lstsq(C, rhs - Tdug, is_transpose=True)

        return delta_u

    def _increment_augmented(
        self,
        apply_T: Callable,
        apply_P: Optional[Callable],
        current_res: np.ndarray,
        current_sol: np.ndarray,
        **mf_args,
    ) -> np.ndarray:
        """
        The linear system to solve is:
        |T + rho*C.T@C    C.T| |du  | |res_u  |
        |      C           0 | |dlam|=|res_lam|
        Args:
            apply_T: it is the tangent matrix
            apply_P: it is a 'good' preconditioner for (T + rho*C.T@C)
        """
        C = self.constraint_matrix
        RHO = self.config.lagrange_penalty
        NR = C.shape[0]

        def matvec_11(x: np.ndarray, **args) -> np.ndarray:
            "Computes (T + rho*C.T@C)@x"
            Tv = apply_T(x, **args)
            Cv = C @ x
            CT_C_v = C.T @ Cv
            y = Tv + RHO * CT_C_v
            return y

        def alm_matvec(x: np.ndarray, **args) -> np.ndarray:
            x_mask = x.copy()
            self._apply_mask(x_mask)
            mvup = matvec_11(x_mask[:-NR], **args) + C.T @ x_mask[-NR:]
            mvdw = C @ x_mask[:-NR]
            y = np.hstack((mvup, mvdw))
            self._apply_mask(y)
            return y

        # TODO: propose a better preconditioner (?)
        def block_preconditioner(x: np.ndarray) -> np.ndarray:
            x_mask = x.copy()
            self._apply_mask(x_mask)
            y = np.zeros_like(x_mask)
            dlam = -RHO * x_mask[-NR:]
            rhs_u = x_mask[:-NR] - C.T @ dlam
            y[:-NR] = apply_P(rhs_u) if callable(apply_P) else rhs_u
            self._apply_mask(y)
            y[-NR:] = dlam
            return y

        # RHS already defined
        res_lam = self._compute_constraint_residual(current_sol)
        res = np.hstack((current_res.copy(), res_lam))
        self._apply_mask(res)

        delta = self._solve(alm_matvec, res, block_preconditioner, **mf_args)
        self._apply_mask(delta)

        # Update lambda (ALM)
        self._lag_mult += delta[-NR:]

        return delta[:-NR]

    def compute_residual(
        self,
        current_res: np.ndarray,
        current_sol: Optional[np.ndarray] = None,
    ):
        """
        Compute the residual incorporating the constraints.
        Args:
            current_res (np.ndarray): The original residual without constraints.
            current_sol (np.ndarray, optional): Current solution vector, required for augmented method.
            **kwargs: Additional arguments if needed.
        """

        # Current residual =  f - A @ U
        updated_res = current_res.copy()

        # For Lagrange residual f - A @ U - C^T @ lambda - rho C^T @ (C u - g)
        # For standard residual f - A @ U - C^T @ lambda
        # NOTE: The last statement is true because we will later multiply Z.T @ res
        # and by definition of the null space Z.T @ C.T = (C @ Z).T should be zero
        C = self.constraint_matrix
        lam = self.lagrange_multiplier
        updated_res -= C.T @ lam

        if self.lagrange_type == "augmented":
            if not isinstance(current_sol, np.ndarray):
                raise ValueError(
                    "Current solution should be array in augmented Lagrangian method"
                )
            u = current_sol.copy()
            self._apply_mask(u)
            rho = self.config.lagrange_penalty
            res_lam = self._compute_constraint_residual(u)
            updated_res += rho * (C.T @ res_lam)

        # Apply Dirichlet mask
        self._apply_mask(updated_res)
        return updated_res

    def compute_increment(
        self,
        apply_T: Callable,
        apply_P: Optional[Callable],
        current_res: np.ndarray,
        current_sol: np.ndarray,
        **mf_args,
    ):
        """
        Compute the increment based on the selected Lagrange method.
        Args:
            apply_T (callable): Function to apply the tangent operator T.
            apply_P (callable, optional): Function to apply the preconditioner P.
            current_sol (np.ndarray): Current solution vector.
            current_res (np.ndarray): Current residual vector.
            **mf_args: Additional arguments for the solver.
        """

        lag_type = self.lagrange_type
        func = {
            "standard": self._increment_standard,
            "augmented": self._increment_augmented,
        }[lag_type]
        residual = current_res.copy()
        solution = current_sol.copy()
        self._apply_mask(residual)
        self._apply_mask(solution)
        return func(apply_T, apply_P, residual, solution, **mf_args)
