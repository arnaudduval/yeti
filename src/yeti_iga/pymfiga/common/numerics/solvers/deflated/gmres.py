"""
GCRO-DR: GCRO with Deflated Restarting.

Translated from the MATLAB implementation adapted from:
    Parks, Michael L., et al. "Recycling Krylov subspaces for sequences of
    linear systems." SIAM Journal on Scientific Computing 28.5 (2006): 1651-1674.
"""

from yeti_iga.pymfiga.common.base.enum import Constants
from ..core import SolverArgs
from ..convergence_manager import ConvergenceManager
from .helpers import gmres1, gmres2, getHarmVecs1, getHarmVecs2
from typing import Optional, Callable, Dict
from scipy.linalg import solve_triangular
from dataclasses import dataclass
import numpy as np
import logging

logger = logging.getLogger("SRC.SOLVER")


@dataclass
class GcrodrArgs(SolverArgs):
    krylov_size: int
    recycled_size: int
    max_cycles: int

    def __post_init__(self):
        if self.krylov_size <= 0:
            raise ValueError("Krylov size should be positive")

        if self.recycled_size < 0:
            # NOTE: 0 means that gcrodr does not save any information
            raise ValueError("Recycle size should be greater than or equal to 0")

        if self.max_cycles <= 0:
            raise ValueError("Max no. cycles should greater than 0")

        if self.krylov_size < self.recycled_size:
            raise ValueError("Not consistent with the method")


# ---------------------------------------------------------------------------
# Main solver class
# ---------------------------------------------------------------------------


class GCRODR:
    """
    GCRO with Deflated Restarting.

    Parameters
    ----------
    tolerance : relative residual tolerance
    maxiters  : total number of iterations (related to number of matvec)
    m         : maximum Krylov subspace dimension per cycle
    k         : number of approximate eigenvectors recycled across cycles / calls
    max_cycles: fixed number of cycles (to avoid a forever while loop)

    Usage
    -----
    solver = GCRODR(A, m=30, k=10)
    x, resvec, r, nmv, relres = solver.solve(b)
    x, resvec, r, nmv, relres = solver.solve(b2)   # reuses recycled subspace
    """

    _convergence_manager = None
    A = B = None

    def __init__(
        self,
        tolerance: float,
        maxiters: int,
        m: int = 40,
        k: int = 10,
        max_cycles: int = 10,
        verbose: bool = True,
    ):
        self.verbose = verbose
        self._U_persist: Dict[str, np.ndarray] = {}
        self._config = GcrodrArgs(
            tolerance=tolerance,
            maxiters=maxiters,
            krylov_size=m,
            recycled_size=k,
            max_cycles=max_cycles,
        )

    @property
    def convergence_manager(self) -> ConvergenceManager:
        if self._convergence_manager is None:
            cm = ConvergenceManager(
                refe_relative=self.config.tolerance,
                refe_iteration=self.config.maxiters,
                refe_absolute=Constants.SAFEGUARD,
            )
            self._convergence_manager = cm
        return self._convergence_manager

    @property
    def config(self):
        return self._config

    def update(
        self,
        tolerance: Optional[float] = None,
        maxiters: Optional[float] = None,
    ):
        """Update the solver parameters."""
        kwargs = self.config.__dict__.copy()
        if tolerance is not None:
            kwargs["tolerance"] = tolerance
        if maxiters is not None:
            kwargs["maxiters"] = maxiters
        self._config = GcrodrArgs(**kwargs)
        self._convergence_manager = None

    def reset(self, reuse_name: str = "default"):
        """Discard the recycled subspace for the given key."""
        self._U_persist.pop(reuse_name, None)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _apply_A(self, x):
        # A can not be None
        return self.A(x) if callable(self.A) else self.A @ x

    def _apply_M(self, v):
        """Apply preconditioner  M v."""
        if self.M is not None:
            v = self.M(v) if callable(self.M) else self.M @ v
        # Else identity
        return v

    def _send_message(self):
        if self.verbose:
            stat = self.convergence_manager.get_status()
            it = stat["iteration"]["current"]
            it = -1 if it is None else it
            abs = stat["absolute_error"]["current"] or 0.0
            rel = stat["relative_error"]["current"] or 0.0
            message = f"""
                Convergence summary in GCRODR:
                - iteration {it+1},
                - abs residue {abs:.2e}
                - rel residue {rel:.2e}
            """
            logger.info(message)

    # ------------------------------------------------------------------
    # Public solver
    # ------------------------------------------------------------------

    def solve(
        self,
        Afun: Callable,
        b: np.ndarray,
        Pfun: Optional[Callable] = None,
        x0: Optional[np.ndarray] = None,
        reuse_name: str = "default",
    ) -> Dict:
        """
        Solve  A x = b.

        Parameters
        ----------
        b           : right-hand side vector
        x0          : initial guess (zeros if None)
        reuse_name  : key for the recycled subspace

        Returns
        -------
        A dictionary containing:
        - 'sol': The solution vector x.
        - 'res': An array of absolute residuals at each iteration.
        """
        self.convergence_manager.clear()
        self.A, self.M = Afun, Pfun
        m, k = self.config.krylov_size, self.config.recycled_size
        tolerance, max_cycles = self.config.tolerance, self.config.max_cycles

        dtype = np.result_type(b, np.complex128)
        if x0 is None:
            x0 = np.zeros_like(b, dtype=dtype)

        x = np.zeros_like(x0, dtype=dtype)  # offset from x0
        nmv = 1

        r = self._apply_M(b - self._apply_A(x0))
        resvec = [float(np.linalg.norm(r))]
        if self.verbose:
            logger.info(f"||r|| = {resvec[0]:.6e}\t\tnmv = {nmv - 1}")

        self.convergence_manager.update(curr_absolute=float(resvec[-1]))
        if self.convergence_manager.has_converged():
            logger.info("External force almost zero. No iterations")
            return {"sol": np.real(x), "res": np.array([1.0])}

        bnorm = (
            np.linalg.norm(self._apply_M(b.copy()))
            if self.M is not None
            else np.linalg.norm(b)
        )

        # ---- Initialise / recycle U ------------------------------------
        U = C = None
        if reuse_name in self._U_persist:
            # ---- Branch A: warm-start from a previous call -------------
            U = self._U_persist[reuse_name].copy()

            # C = M⁻¹ A U  (recompute; handles A varying between calls)
            C = np.column_stack(
                [self._apply_M(self._apply_A(U[:, i])) for i in range(U.shape[1])]
            )
            # Orthonormalise C; adjust U so C = A U still holds: Q = A (U/R)
            C, R = np.linalg.qr(C, mode="reduced")
            U = solve_triangular(R.T, U.T, lower=True).T  # U = U / R

            Cr = C.conj().T.dot(r)
            x = x + U @ Cr
            r = r - C @ Cr
            resvec[0] = float(np.linalg.norm(r))

            self.convergence_manager.update(
                curr_absolute=float(resvec[-1]),
                curr_relative=float(resvec[-1] / bnorm),
                curr_iteration=1,
            )

        else:
            # ---- Branch B: prime the pump with one plain GMRES cycle ---
            x, r, V, H, p, rv = gmres1(
                self._apply_A, x, r, m, self._apply_M, float(tolerance * bnorm)
            )
            resvec.extend(rv.tolist())
            nmv += p

            self.convergence_manager.update(
                curr_absolute=float(resvec[-1]),
                curr_relative=float(resvec[-1] / bnorm),
                curr_iteration=nmv,
            )

            if self.verbose:
                logger.info(f"||r|| = {resvec[-1]:.6e}\t\tnmv = {nmv - 1}")

            if k < p:
                P = getHarmVecs1(p, k, H)
                U = V[:, :p] @ P
                C_coeff, R = np.linalg.qr(H[: p + 1, :p] @ P, mode="reduced")
                C = V[:, : p + 1] @ C_coeff  # lift back to full space
                U = solve_triangular(R.T, U.T, lower=True).T

            # Early convergence during the priming cycle
            if p < m:
                self._send_message()
                x = x0 + x
                nmv -= 1
                if k < p and U is not None:
                    self._U_persist[reuse_name] = U
                return {"sol": np.real(x), "res": np.asarray(resvec)}

        # ---- Main GCRO-DR loop -----------------------------------------

        assert U is not None and C is not None
        currcycle = 0
        for currcycle in range(max_cycles):

            if self.convergence_manager.has_converged():
                self._send_message()
                break

            V, H_inner, B, p, rv = gmres2(
                self._apply_A, r, m - k, self._apply_M, C, float(tolerance * bnorm)
            )
            resvec.extend(rv.tolist())
            nmv += p

            self.convergence_manager.update(
                curr_absolute=float(resvec[-1]),
                curr_relative=float(resvec[-1] / bnorm),
                curr_iteration=nmv,
            )

            if self.verbose:
                logger.info(f"||r|| = {resvec[-1]:.6e}\t\tnmv = {nmv - 1}")

            # Column-normalise U; record inverse norms in D
            col_norms = np.linalg.norm(U, axis=0)
            U = U / col_norms
            D = np.diag(1.0 / col_norms)

            # Augmented Hessenberg  H2  of shape (p+k+1, p+k)
            #
            #   [ D     B[:, :p] ]   ← k rows
            #   [ 0     H[:p+1]  ]   ← p+1 rows
            #
            H2 = np.block([[D, B[:, :p]], [np.zeros((p + 1, k)), H_inner[: p + 1, :p]]])

            # [C, V[:,0:p+1]] has p+k+1 columns matching H2's row count
            CV_full = np.hstack([C, V[:, : p + 1]])  # n × (p+k+1)
            # [U, V[:,0:p]] are the actual search-direction vectors for x update
            UV = np.hstack([U, V[:, :p]])  # n × (p+k)

            rhs = CV_full.conj().T.dot(r)  # (p+k+1,)
            y, *_ = np.linalg.lstsq(H2, rhs, rcond=None)

            # x update:  use UV  (actual search directions, NOT A-images)
            x = x + UV @ y
            # r update:  use CV_full  (consistent with the GCRODR relation A·UV = CV_full·H2)
            r = r - CV_full @ (H2 @ y)

            # Early convergence within a cycle
            if p < m - k:
                self._send_message()
                break

            # Harmonic Ritz extraction for the next recycled subspace
            P = getHarmVecs2(p + k, k, H2, V[:, : p + 1], U, C)

            U_new = UV @ P
            Q2, R2 = np.linalg.qr(H2 @ P, mode="reduced")
            C = CV_full @ Q2
            U = np.linalg.solve(R2.T, U_new.T).T

        else:
            logger.warning(f"No convergence after {currcycle} cycles")

        # ---- Converged -------------------------------------------------
        x = x0 + x
        nmv -= 1
        self._U_persist[reuse_name] = U

        return {"sol": np.real(x), "res": np.asarray(resvec)}
