from yeti_iga.pymfiga.common.base.enum import Constants
from .core import SolverArgs
from .convergence_manager import ConvergenceManager
from .inner_tolerance import InnerToleranceSetter
from .update_manager import UpdateManager
from typing import Callable, Dict, List, Optional
from dataclasses import dataclass
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.SOLVER")


@dataclass
class NonLinearArgs(SolverArgs):
    allow_acceleration: bool
    allow_line_search: bool
    linesearch_maxiters: int = 4
    linesearch_miniter: int = 0
    anderson_maxsize: int = 4
    anderson_miniter: int = 1

    def __post_init__(self):
        """Validate arguments after initialization."""
        if not isinstance(self.allow_acceleration, bool):
            raise ValueError(
                f"It should be boolean, got {type(self.allow_acceleration)}"
            )

        if not isinstance(self.allow_line_search, bool):
            raise ValueError(
                f"It should be boolean, got {type(self.allow_line_search)}"
            )

        # Other requirements
        assert self.linesearch_maxiters > 0
        assert self.linesearch_miniter >= 0
        assert self.anderson_maxsize > 0
        assert self.anderson_miniter >= 0
        assert self.anderson_miniter < self.anderson_maxsize


class NonLinearSolver:
    """
    A class for solving nonlinear problems of the form res(x) = 0 using iterative methods.
    The solver supports Anderson acceleration and backtracking line search to enhance convergence.
    """

    # Internal variables
    _update_manager = None
    _inner_tolerance_setter = None
    _convergence_manager = None

    def __init__(
        self,
        tolerance: float,
        maxiters: int,
        allow_acceleration: bool,
        allow_line_search: bool,
        verbose: bool = True,
    ):
        """
        Parameters:
            tolerance (float): Convergence tolerance for the nonlinear solver
            maxiters (int): Maximum number of iterations allowed
            allow_acceleration (float): Whether to use Anderson acceleration
            allow_line_search (float): Whether to use backtracking line search
            verbose (bool): Whether to print convergence information during iterations
        """
        logger.debug("Initialize nonlinear solver")
        # Internal variables
        self._config = NonLinearArgs(
            tolerance=tolerance,
            maxiters=maxiters,
            allow_acceleration=allow_acceleration,
            allow_line_search=allow_line_search,
        )
        self._verbose = verbose

    @property
    def config(self) -> NonLinearArgs:
        return self._config

    @property
    def update_manager(self):
        if not isinstance(self._update_manager, UpdateManager):
            self._update_manager = UpdateManager()
        return self._update_manager

    @property
    def inner_tolerance_setter(self) -> InnerToleranceSetter:
        if not isinstance(self._inner_tolerance_setter, InnerToleranceSetter):
            inner_tolerance = InnerToleranceSetter(
                inner_tolerance_type="exact picard",
                inner_tolerance_args={"default": Constants.TINY},
            )
            self._inner_tolerance_setter = inner_tolerance
        return self._inner_tolerance_setter

    @property
    def convergence_manager(self) -> ConvergenceManager:
        if self._convergence_manager is None:
            cm = ConvergenceManager(
                refe_relative=self.config.tolerance,
                refe_iteration=self.config.maxiters,
                refe_absolute=Constants.SAFEGUARD,
            )
            cm.add_criterion("curr_increment", "relative_error")
            self._convergence_manager = cm
        return self._convergence_manager

    def update(
        self,
        tolerance: Optional[float] = None,
        maxiters: Optional[float] = None,
        allow_acceleration: Optional[bool] = None,
        allow_line_search: Optional[bool] = None,
        inner_tolerance_setter: Optional[InnerToleranceSetter] = None,
        update_manager: Optional[UpdateManager] = None,
    ):
        """
        Update solver parameters.
        Parameters:
            tolerance (float, optional): New convergence tolerance
            maxiters (int, optional): New maximum number of iterations
            allow_acceleration (bool, optional): Whether to use Anderson acceleration
            allow_line_search (bool, optional): Whether to use backtracking line search
        """
        kwargs = self._config.__dict__.copy()
        if tolerance is not None:
            kwargs["tolerance"] = tolerance
        if maxiters is not None:
            kwargs["maxiters"] = maxiters
        if allow_acceleration is not None:
            kwargs["allow_acceleration"] = allow_acceleration
        if allow_line_search is not None:
            kwargs["allow_line_search"] = allow_line_search
        self._config = NonLinearArgs(**kwargs)
        self._convergence_manager = None
        if isinstance(inner_tolerance_setter, InnerToleranceSetter):
            self._inner_tolerance_setter = inner_tolerance_setter
        if isinstance(update_manager, UpdateManager):
            self._update_manager = update_manager

    def _anderson_acceleration(
        self, incr_hist: List[np.ndarray], sol_hist: List[np.ndarray]
    ) -> np.ndarray:
        """
        The problem to be solved is A(u) * u = b
        The fixed point is A(u_k) * u_{k+1} = b, with u_0 = 0
        The problem can be also written as Newton: u_{k+1} = u_k + incr_k
        where incr_k = (A(u_k))^{-1} @ res_k and res_k = b - A(u_k) * u_k

        For Anderson acceleration we define f_k = (A(u_k))^{-1} @ b,
        and g_k = f_k - u_k = (A(u_k))^{-1} @ b - u_k = (A(u_k))^{-1} @ res_k = incr_k

        Let define G = [g_{k-m+1} - g_{k-m}, ..., g_k - g_{k-1}], then we compute
        the solution of min || G @ alpha - g_k||^2 over alpha in R^m using QR decomposition
        Finally, let define X = [u_{k-m+1} - u_{k-m}, ..., u_k - u_{k-1}]
        the new increment is g_k - alpha @ (G + X)

        Parameters:
            incr_hist (list): history of increments g_i for i = k-m, ..., k
            sol_hist (list): history of solutions u_i for i = k-m, ..., k

        Returns accelerated solution.
        """

        logger.debug("Starting fixed-point acceleration with Anderson algorithm")

        # Number of previous directions available
        m = len(incr_hist) - 1
        if m < self.config.anderson_miniter:
            # Not enough history, fallback to fixed point
            return incr_hist[-1]

        # Stack g_k into matrix G (ndof × m))
        G = np.column_stack([incr_hist[i + 1] - incr_hist[i] for i in range(m)])
        X = np.column_stack([sol_hist[i + 1] - sol_hist[i] for i in range(m)])

        # Right-hand side is current increment g_k
        g_k = incr_hist[m]

        # Solve least squares problem: min ||G * gamma - g_k||
        alpha = np.linalg.lstsq(G, g_k, rcond=None)[0]

        if np.linalg.norm(alpha) > Constants.HIGH:
            # Value too high fallback to fixed point
            return incr_hist[-1]

        # Compute accelerated increment
        return g_k - (G + X) @ alpha

    def _backtracking_line_search(
        self,
        sol_current: np.ndarray,
        incr_current: np.ndarray,
        compute_residual: Callable,
        **res_args,
    ) -> float:
        """
        Backtracking line search using the Armijo condition.

        The merit function is:
            phi(x) = 0.5 * ||R(x)||^2

        For a Newton step, the directional derivative of phi at alpha=0
        can be approximated as:
            phi'(0) ≈ -2 * phi(0)

        Parameters:
            sol_current (np.ndarray): Current solution vector.
            incr_current (np.ndarray): Newton increment (search direction).
            compute_residual (Callable): Function that returns the residual vector.
            res_args (dict): Additional arguments passed to compute_residual.

        Returns:
            Step length alpha.
        """

        # 1. Line search parameters
        alpha = 1.0  # Initial full Newton step
        rho = 0.5  # Step reduction factor
        c1 = 1e-4  # Armijo sufficient decrease parameter
        alpha_min = Constants.TINY  # Minimum allowed step size

        logger.debug("Starting Backtracking Line Search (Armijo)")

        # 2. Evaluate residual at current solution
        res_0 = compute_residual(sol_current, **res_args)[0]
        phi_0 = 0.5 * np.linalg.norm(res_0) ** 2

        # Directional derivative approximation for Newton methods
        # This assumes a reasonably accurate Newton direction
        phi_prime_0 = -2.0 * phi_0

        for iteration in range(self.config.linesearch_maxiters):

            # 3. Trial step
            sol_trial = sol_current + alpha * incr_current

            # 4. Evaluate residual at trial solution
            res_trial = compute_residual(sol_trial, **res_args)[0]
            phi_trial = 0.5 * np.linalg.norm(res_trial) ** 2

            # 5. Armijo condition
            if phi_trial <= phi_0 + c1 * alpha * phi_prime_0:
                logger.debug(
                    "Line search accepted: alpha=%.3e, phi=%.3e, reductions=%d",
                    alpha,
                    phi_trial,
                    iteration,
                )
                return float(alpha)

            # 6. Reduce step length
            alpha *= rho

            logger.debug(
                "Line search reduction: alpha=%.3e, phi_trial=%.3e", alpha, phi_trial
            )

            # 7. Safeguard against excessively small steps
            if alpha < alpha_min:
                logger.debug(
                    "Line search reached minimum step size (alpha=%.3e). "
                    "Returning smallest step.",
                    alpha,
                )
                return float(alpha)

        logger.debug(
            "Line search reached maximum iterations. Returning alpha=%.3e", alpha
        )
        return float(alpha)

    def solve(
        self,
        solution: np.ndarray,
        compute_residual: Callable,
        compute_increment: Callable,
        residual_args: Dict = {},
        increment_args: Dict = {},
    ) -> dict:
        """
        Solve the nonlinear problem res(x) = 0 using an iterative method.
        Args:
            solution (np.ndarray): Initial guess for the solution, updated in-place.
            compute_residual (Callable): Function that computes the residual vector given a solution.
            compute_increment (Callable): Function that computes the increment (search direction) given a residual.
            residual_args (Dict): Additional arguments to pass to compute_residual.
            increment_args (Dict): Additional arguments to pass to compute_increment.
        """

        assert callable(compute_residual) and callable(compute_increment)
        # Clear current values of Convergence
        self.convergence_manager.clear()
        self.update_manager.reset_newton()

        # Variables for inner and outer tolerance
        norm_increment, norm_solution = 1.0, 1.0
        norm_residual_old = None
        inner_tolerance_old = None
        norm_residual_ref = 0.0  # Dummy initialization

        # Variable for Anderson acceleration
        incr_history: List = []
        sol_history: List = []

        # Save data
        nonlinear_residual_list = []
        solution_history_list = {}
        nonlinear_time_list = []
        nonlinear_rate_list = []
        linear_tolerance_list = []
        extra_args = {}

        start = time()
        if self._verbose:
            logger.info(repr(self))
            logger.info(f"ITERATION | ABS RESIDUAL")
            logger.info("=" * 30)

        for it in range(self.config.maxiters):
            residual, extra_args = compute_residual(solution, **residual_args)
            norm_residual = float(np.linalg.norm(residual))
            if self._verbose:
                logger.info(f" IT {it} : ABS {norm_residual:.2e}")
            else:
                logger.debug(f"Iteration {it} and residual {norm_residual:.2e}")
            nonlinear_residual_list.append(norm_residual)
            solution_history_list[f"noniter_{it}"] = np.copy(solution)
            nonlinear_time_list.append(time() - start)

            if it == 0:
                norm_residual_ref = norm_residual

            self.convergence_manager.update(
                curr_iteration=int(it),
                curr_absolute=float(norm_residual),
                curr_relative=float(
                    norm_residual / (norm_residual_ref + Constants.SAFEGUARD)
                ),
                curr_increment=float(
                    norm_increment / (norm_solution + Constants.SAFEGUARD)
                ),
            )
            if self.convergence_manager.has_converged():
                if self._verbose:
                    stat = self.convergence_manager.get_status()
                    it = stat["iteration"]["current"] or 0
                    abs = stat["absolute_error"]["current"] or 0.0
                    rel = stat["relative_error"]["current"] or 0.0
                    incr = stat["curr_increment"]["current"] or 0.0
                    message = f"""
                        Convergence summary in NONLINEAR SOLVER:
                        - iteration {it},
                        - abs residue {abs:.2e}
                        - rel residue {rel:.2e}
                        - rel increment {incr:.2e}
                        Execution time: {time() - start:.2e} seconds.
                    """
                    logger.info(message)
                break

            inner_tolerance = self.inner_tolerance_setter.compute(
                norm_residual,
                norm_residual_old,
                inner_tolerance_old,
            )
            norm_residual_old = norm_residual
            inner_tolerance_old = inner_tolerance
            linear_tolerance_list.append(inner_tolerance)

            increment_args.update(extra_args)
            increment_args.update({"inner_tolerance": inner_tolerance})
            increment_args.update({"current_solution": solution})
            increment: np.ndarray = compute_increment(residual, **increment_args)

            if self.config.allow_acceleration:
                if it >= self.config.anderson_maxsize:
                    incr_history.pop(0)
                    sol_history.pop(0)
                incr_history.append(increment.copy())
                sol_history.append(solution.copy())
                increment = self._anderson_acceleration(incr_history, sol_history)

            if self.config.allow_line_search and it > self.config.linesearch_miniter:
                alpha = self._backtracking_line_search(
                    solution,
                    increment,
                    compute_residual,
                    **residual_args,
                )
                increment *= alpha

            # Update active control points
            solution += increment

            # Compute norms for outer stopping criterion
            norm_increment = np.linalg.norm(increment)
            norm_solution = np.linalg.norm(solution)
            nonlinear_rate_list.append(float(norm_increment / norm_solution))

            # Update manager
            self.update_manager.increment_newton()
        else:
            it = self.convergence_manager.parameters["iteration"].current
            logger.warning(f"No convergence after {it} iterations")

        all_extra_args = {
            "nonlinear_residual_list": nonlinear_residual_list,
            "nonlinear_time_list": nonlinear_time_list,
            "nonlinear_rate_list": nonlinear_rate_list,
            "linear_tolerance_list": linear_tolerance_list,
            "solution_history_list": solution_history_list,
            "extra_args": extra_args,
        }

        return all_extra_args

    def __repr__(self) -> str:
        message = f""""
            \nNON LINEAR SOLVER:
            Max iterations: {self.config.maxiters}
            Tolerance: {self.config.tolerance}
            With Anderson acceleration: {self.config.allow_acceleration}
            With line search: {self.config.allow_line_search}
        """
        return message
