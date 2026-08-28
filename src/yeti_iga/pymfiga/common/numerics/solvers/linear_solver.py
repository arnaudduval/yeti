from yeti_iga.pymfiga.common.base.math import Constants
from .core import SolverArgs
from .convergence_manager import ConvergenceManager
from typing import Callable, Tuple, Dict, List, Optional, Union, Literal
from scipy.sparse import linalg as scsplin
from scipy import sparse as sp
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.SOLVER")


class LinearSolver:
    """
    A class for solving linear systems of equations using various iterative methods,
    including Conjugate Gradient (CG), BiConjugate Gradient Stabilized (BiCGSTAB),
    and Generalized Minimal Residual Method (GMRES).
    The class also provides direct solvers using sparse LU decomposition.
    """

    _VALID_TYPES = ["cg", "bicgstab", "gmres"]
    _convergence_manager = None

    def __init__(
        self,
        tolerance: float,
        maxiters: int,
        cleandod: Optional[List[int]] = None,
        linear_type: Literal["cg", "bicgstab", "gmres"] = "gmres",
        verbose: bool = True,
    ):
        """
        Args:
            tolerance (float): The convergence tolerance for the iterative solver.
            maxiters (int): The maximum number of iterations for the iterative solver.
            cleandod (List[int], optional): A list of indices to be masked (set to zero)
                in the residual vector during the iterative process. Default is None.
            linear_type (str, optional): The type of iterative solver to use.
                Must be one of 'cg', 'bicgstab', or 'gmres'. Default is 'gmres'.
            verbose (bool, optional): If True, prints convergence information and execution time. Default is True.

        """
        logger.debug("Initialize linear solver")
        # Internal variables
        self._linear_type = self._set_linear_type(linear_type)
        self._cleandod = cleandod
        self._config = SolverArgs(
            tolerance=tolerance,
            maxiters=maxiters,
        )
        self._Pfun_default = lambda x: x
        self._verbose = verbose

    @property
    def linear_type(self):
        return self._linear_type

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        if isinstance(value, bool):
            self._verbose = value

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

    def set_cleandod(self, cleandod: Optional[List[int]]):
        if cleandod is not None and not isinstance(cleandod, list):
            raise ValueError("Cleandod should be list")
        self._cleandod = cleandod

    def _apply_mask(self, x):
        if self._cleandod is not None:
            x[self._cleandod] = 0.0

    def _set_linear_type(self, value: str):
        if not isinstance(value, str):
            raise TypeError(f"Linear type must be string, got {type(value)}")

        if value not in self._VALID_TYPES:
            raise ValueError(
                f"Invalid tolerance type '{value}'. "
                f"Must be one of {self._VALID_TYPES}"
            )
        return value

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
        self._config = SolverArgs(**kwargs)
        self._convergence_manager = None

    def solve(
        self,
        Afun: Callable,
        b: np.ndarray,
        Pfun: Optional[Callable] = None,
        **mf_args,
    ) -> Dict[str, np.ndarray]:
        """
        Solve the linear system Ax = b using the specified iterative method.

        Parameters:
        -----------
        Afun : callable
                A function that computes the matrix-vector product A*x.
        b : array_like
                The right-hand side vector.
        Pfun : callable, optional
                A function that computes the preconditioner matrix-vector product P*x. Default is the identity function.
        kwargs : dict, optional
                Additional arguments to pass to Afun.

        Returns:
        --------
        A dictionary containing:
        - 'sol': The solution vector x.
        - 'res': An array of relative residuals at each iteration.
        """
        self.convergence_manager.clear()
        # TODO: maybe use scipy iterative solver with callbacks
        solver_list = {
            "cg": self._CG,
            "bicgstab": self._BICGSTAB,
            "gmres": self._GMRES,
        }
        solver = solver_list[self.linear_type]
        start = time()
        output = solver(Afun, b, Pfun=Pfun, **mf_args)
        if self._verbose:
            stat = self.convergence_manager.get_status()
            it = stat["iteration"]["current"]
            it = -1 if it is None else it
            abs = stat["absolute_error"]["current"] or 0.0
            rel = stat["relative_error"]["current"] or 0.0
            message = f"""
                Convergence summary in LINEAR SOLVER:
                - iteration {it+1},
                - abs residue {abs:.2e}
                - rel residue {rel:.2e}
                Execution time: {time() - start:.2e} seconds.
            """
            logger.info(message)
        return output

    def _CG(
        self,
        Afun: Callable,
        b: np.ndarray,
        Pfun: Optional[Callable] = None,
        x0: Optional[np.ndarray] = None,
        **mf_args,
    ) -> Dict[str, np.ndarray]:
        """
        Conjugate Gradient (CG) solver for solving the linear system Ax = b.

        Args:

            Afun (Callable) :
                Function that computes the matrix-vector product A @ x.
            b (np.ndarray) :
                Right-hand side vector.
            Pfun (Optional, Callable) :
                Preconditioner function. Default is the identity function.
            x0 (Optional, np.ndarray) :
                Initial guess for the solution. Default is None, which means a zero vector will be used
            mf_args (Optional, dict) :
                Additional arguments to pass to Afun. Default is an empty dictionary.

        Returns:

            (dict) :
            - 'sol' (np.ndarray), Solution vector
            - 'res' (list), List of relative residuals at each iteration.

        """
        Pfun = self._Pfun_default if Pfun is None else Pfun
        assert callable(Afun) and callable(Pfun)

        if x0 is not None:
            x = np.copy(x0)
            self._apply_mask(x)
            r = b - Afun(x, **mf_args)
        else:
            x = np.zeros_like(b)
            r = b.copy()

        self._apply_mask(r)
        norm_0 = float(np.linalg.norm(r))
        self.convergence_manager.update(curr_absolute=norm_0)
        if self.convergence_manager.has_converged():
            logger.info("External force almost zero. No iterations")
            return {"sol": x, "res": np.array([1.0])}
        logger.debug(f"CG: Iteration {0} and residual {1.0}")

        ptilde = Pfun(r)
        self._apply_mask(ptilde)
        p = np.copy(ptilde)
        rsold = np.dot(r, ptilde)
        res = [1.0]

        maxiters = self.convergence_manager.parameters["iteration"].reference
        for it in range(int(maxiters)):
            Ap = Afun(p, **mf_args)
            self._apply_mask(Ap)
            alpha = rsold / np.dot(p, Ap)
            r -= alpha * Ap
            x += alpha * p

            norm_1 = float(np.linalg.norm(r))
            res.append(norm_1 / norm_0)
            logger.debug(f"CG: Iteration {it+1} and residual {res[-1]:.2e}")
            self.convergence_manager.update(
                curr_iteration=it, curr_absolute=norm_1, curr_relative=res[-1]
            )
            if self.convergence_manager.has_converged():
                break

            ptilde = Pfun(r)
            self._apply_mask(ptilde)
            rsnew = np.dot(r, ptilde)
            p = ptilde + rsnew / rsold * p
            rsold = np.copy(rsnew)
        else:
            it = self.convergence_manager.parameters["iteration"].current or 0
            logger.warning(f"No convergence after {it} iterations")
            
        if x0 is not None:
            x += x0
    

        return {"sol": x, "res": np.array(res)}  


    def _BICGSTAB(
        self,
        Afun: Callable,
        b: np.ndarray,
        Pfun: Optional[Callable] = None,
        x0: Optional[np.ndarray] = None,
        **mf_args,
    ) -> Dict[str, np.ndarray]:
        """
        BiConjugate Gradient Stabilized (BiCGSTAB) method for solving linear systems.

        Args:

            Afun (Callable) :
                Function that computes the matrix-vector product A @ x.
            b (np.ndarray) :
                Right-hand side vector.
            Pfun (Optional, Callable) :
                Preconditioner function. Default is the identity function.
            x0 (Optional, np.ndarray) :
                Initial guess for the solution. Default is None, which means a zero vector will be used
            mf_args (Optional, dict) :
                Additional arguments to pass to Afun. Default is an empty dictionary.

        Returns:

            (dict) :
            - 'sol' (np.ndarray), Solution vector
            - 'res' (list), List of relative residuals at each iteration.

        Notes:
        ------
            This method iteratively solves the linear system A*x = b using the BiCGSTAB algorithm.
            The algorithm is suitable for large, sparse, and non-symmetric linear systems.

        """
        Pfun = self._Pfun_default if Pfun is None else Pfun
        assert callable(Afun) and callable(Pfun)

        if x0 is not None:
            x = np.copy(x0)
            self._apply_mask(x)
            r = b - Afun(x, **mf_args)
        else:
            x = np.zeros_like(b)
            r = b.copy()

        self._apply_mask(r)
        norm_0 = float(np.linalg.norm(r))
        self.convergence_manager.update(curr_absolute=norm_0)
        if self.convergence_manager.has_converged():
            logger.info("External force almost zero. No iterations")
            return {"sol": x, "res": np.array([1.0])}
        logger.debug(f"BICGSTAB: Iteration {0} and residual {1.0}")

        rhat = r.copy()
        p = r.copy()
        rsold = np.dot(r, rhat)
        res = [1.0]

        maxiters = self.convergence_manager.parameters["iteration"].reference
        for it in range(int(maxiters)):
            ptilde = Pfun(p)
            self._apply_mask(ptilde)
            Aptilde = Afun(ptilde, **mf_args)
            self._apply_mask(Aptilde)
            alpha = rsold / np.dot(Aptilde, rhat)
            s = r - alpha * Aptilde
            x += alpha * ptilde

            norm_1 = float(np.linalg.norm(s))
            res.append(norm_1 / norm_0)
            self.convergence_manager.update(
                curr_iteration=it, curr_absolute=norm_1, curr_relative=res[-1]
            )
            if self.convergence_manager.has_converged():
                break

            stilde = Pfun(s)
            self._apply_mask(stilde)
            Astilde = Afun(stilde, **mf_args)
            self._apply_mask(Astilde)
            omega = np.dot(Astilde, s) / np.dot(Astilde, Astilde)
            r = s - omega * Astilde
            x += omega * stilde

            norm_1 = float(np.linalg.norm(r))
            res.append(norm_1 / norm_0)
            logger.debug(f"BICGSTAB: Iteration {it+1} and residual {res[-1]:.2e}")
            self.convergence_manager.update(
                curr_iteration=it, curr_absolute=norm_1, curr_relative=res[-1]
            )
            if self.convergence_manager.has_converged():
                break

            rsnew = np.dot(r, rhat)
            beta = (alpha / omega) * (rsnew / rsold)
            p = r + beta * (p - omega * Aptilde)
            rsold = np.copy(rsnew)
        else:
            it = self.convergence_manager.parameters["iteration"].current or 0
            logger.warning(f"No convergence after {it} iterations")
        if x0 is not None:
            x += x0

        return {"sol": x, "res": np.array(res)}

    def _GMRES(
        self,
        Afun: Callable,
        b: np.ndarray,
        Pfun: Optional[Callable] = None,
        x0: Optional[np.ndarray] = None,
        **mf_args,
    ) -> Dict[str, np.ndarray]:
        """
        Generalized Minimal Residual Method (GMRES) for solving a linear system of equations.

        Args:

            Afun (Callable) :
                Function that computes the matrix-vector product A @ x.
            b (np.ndarray) :
                Right-hand side vector.
            Pfun (Optional, Callable) :
                Preconditioner function. Default is the identity function.
            x0 (Optional, np.ndarray) :
                Initial guess for the solution. Default is None, which means a zero vector will be used
            mf_args (Optional, dict) :
                Additional arguments to pass to Afun. Default is an empty dictionary.

        Returns:

            (dict) :
            - 'sol' (np.ndarray), Solution vector
            - 'res' (list), List of relative residuals at each iteration.

        Notes:
        ------
            This implementation uses the Arnoldi process to build an orthonormal basis
            for the Krylov subspace and solves the least squares problem to minimize the residual.

        """
        Pfun = self._Pfun_default if Pfun is None else Pfun
        assert callable(Afun) and callable(Pfun)

        if x0 is not None:
            x = np.copy(x0)
            self._apply_mask(x)
            r = b - Afun(x, **mf_args)
        else:
            x = np.zeros_like(b)
            r = b.copy()

        self._apply_mask(r)
        norm_0 = float(np.linalg.norm(r))
        self.convergence_manager.update(curr_absolute=norm_0)
        if self.convergence_manager.has_converged():
            logger.info("External force almost zero. No iterations")
            return {"sol": x, "res": np.array([1.0])}
        logger.debug(f"GMRES: Iteration {0} and residual {1.0}")
        maxiters = int(self.convergence_manager.parameters["iteration"].reference)

        Hessenberg = np.zeros((maxiters + 1, maxiters))
        Vectors = np.zeros((maxiters + 1, len(b)))

        Vectors[0] = r / norm_0
        y = np.array([])
        res = [1.0]

        for it in range(maxiters):
            p = Pfun(Vectors[it])
            self._apply_mask(p)
            w = Afun(p, **mf_args)
            self._apply_mask(w)

            for j in range(it + 1):
                Hessenberg[j, it] = np.dot(w, Vectors[j])
                w -= Hessenberg[j, it] * Vectors[j]

            Hessenberg[it + 1, it] = np.linalg.norm(w)
            if Hessenberg[it + 1, it] != 0:
                Vectors[it + 1] = w / Hessenberg[it + 1, it]

            mat = Hessenberg[: it + 2, : it + 1]
            rhs = np.zeros(it + 2)
            rhs[0] = norm_0

            y = np.linalg.lstsq(mat, rhs, rcond=None)[0]
            norm_1 = float(np.linalg.norm(mat @ y - rhs))
            res.append(norm_1 / norm_0)
            logger.debug(f"GMRES: Iteration {it+1} and residual {res[-1]:.2e}")
            self.convergence_manager.update(
                curr_iteration=it, curr_absolute=float(norm_1), curr_relative=res[-1]
            )
            if self.convergence_manager.has_converged():
                break
        else:
            it = self.convergence_manager.parameters["iteration"].current or 0
            logger.warning(f"No convergence after {it} iterations")

        x += Pfun(Vectors[: len(y)].T @ y)
        if x0 is not None:
            x += x0

        return {"sol": x, "res": np.array(res)}

    @staticmethod
    def direct(
        A: Union[sp.csr_array, np.ndarray],
        b: np.ndarray,
        verbose: bool = True,
        **kwargs,
    ):
        """
        Direct solver for linear systems using sparse LU decomposition.

        Args:

            Afun (sp.csr_array, np.ndarray) :
                Coefficient matrix of the linear system.
            b (np.ndarray) :
                Right-hand side vector.
            verbose (Optional, bool) :
                If True, prints the time taken to solve the system. Default is True.
            kwargs (Optional, dict) :
                Additional arguments (not used in this function).

        Returns:

            (dict) :
            - 'sol' (np.ndarray), Solution vector
            - 'res' (list), An empty array (for consistency with iterative solvers).

        """
        start = time()
        if isinstance(A, np.ndarray):
            A = sp.csr_array(A)
            A.eliminate_zeros()
        sol = scsplin.spsolve(A, b)
        if verbose:
            logger.info(f"Direct solver took {time() - start:.2e} seconds.")
        return {"sol": sol, "res": np.array([])}

    @staticmethod
    def iterative(
        A: Union[sp.csr_array, np.ndarray],
        b: np.ndarray,
        x0: Optional[np.ndarray] = None,
        verbose: bool = True,
        linear_type: Literal["gmres", "cg", "bicgstab"] = "gmres",
        **kwargs,
    ):
        """
        Iterative direct solver for linear systems using sparse LU decomposition as a preconditioner.

        Args:

            Afun (sp.csr_array, np.ndarray) :
                Coefficient matrix of the linear system.
            b (np.ndarray) :
                Right-hand side vector.
            x0 (Optional, np.ndarray) :
                Initial guess for the solution. Default is None, which means a zero vector will be used
            verbose (Optional, bool) :
                    If True, prints the time taken to solve the system. Default is True.
            iterative_type (str) : A string specifying the type of iterative solver to use ('gmres', 'cg', or 'bicgstab').
            kwargs (Optional, dict) :
                Additional arguments (not used in this function).

        Returns:

            (dict) :
            - 'sol' (np.ndarray), Solution vector
            - 'res' (list), An empty array (for consistency with iterative solvers).
        """
        start = time()
        if isinstance(A, np.ndarray):
            A = sp.csr_array(A)
            A.eliminate_zeros()
        solver_list = {
            "cg": scsplin.cg,
            "gmres": scsplin.gmres,
            "bicgstab": scsplin.bicgstab,
        }
        assert solver_list.get(linear_type) is not None, NotImplementedError(
            f"Linear solver {linear_type} not implemented"
        )
        ilu = scsplin.spilu(A)
        matvec = lambda x: ilu.solve(x)
        precond = scsplin.LinearOperator(shape=A.shape, matvec=matvec, dtype=A.dtype)
        sol = solver_list[linear_type](A, b, x0=x0, M=precond)[0]
        if verbose:
            logger.info(
                f"Iterative direct solver ({linear_type}) took {time() - start:.2e} seconds."
            )
        return {"sol": sol, "res": np.array([])}

    @staticmethod
    def eigs(
        N: int,
        Afun: Callable,
        Bfun: Optional[Callable] = None,
        Pfun: Optional[Callable] = None,
        k: int = 5,
        which: Literal["LM", "SM"] = "LM",
        **kwargs,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the eigenvalues and eigenvectors of a linear operator.

        Args:

            N (int):
                The size of the square matrix.
            Afun (callable):
                A function that performs the matrix-vector multiplication for the matrix A.
            Bfun (callable, optional):
                A function that performs the matrix-vector multiplication for the matrix B (default is None).
            Pfun (callable, optional):
                A function that performs the matrix-vector multiplication for the preconditioner matrix P (default is None).

        Returns:
            (np.ndarray, np.ndarray):
            - eigvals_sorted, The sorted eigenvalues of the matrix.
            - eigvecs_sorted, The eigenvectors corresponding to the sorted eigenvalues.

        Notes:
        ------
            This function uses the `scipy.sparse.linalg.eigs` method to compute the eigenvalues and eigenvectors.
            The eigenvalues and eigenvectors are sorted in ascending order based on the eigenvalues.

        """
        if k > N - 2:
            k = N - 2 
        ALinOp = scsplin.LinearOperator(shape=(N, N), matvec=Afun, dtype=float)
        BLinOp = (
            None
            if Bfun is None
            else scsplin.LinearOperator(shape=(N, N), matvec=Bfun, dtype=float)
        )
        PLinOp = (
            None
            if Pfun is None
            else scsplin.LinearOperator(shape=(N, N), matvec=Pfun, dtype=float)
        )
        kwargs.update({"k": k, "which": which})
        output = scsplin.eigs(A=ALinOp,M=BLinOp, Minv=PLinOp, **kwargs)

        eigvals = np.real(output[0])
        eigvecs = np.real(output[1])
        sorted_indices = np.argsort(eigvals)
        eigvals_sorted = eigvals[sorted_indices]
        eigvecs_sorted = eigvecs[:, sorted_indices]
        return eigvals_sorted, eigvecs_sorted
    

    @staticmethod
    def power_iteration(
        N: int,
        Afun: Callable,
        Bfun: Optional[Callable] = None,
        Pfun: Optional[Callable] = None,
        maxiter: int = 5000,
        tol: float = 1e-8,
        seed: int = 0,
        **kwargs,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the largest eigenvalue and its associated eigenvector
        using the Power Iteration method.

        For the generalized eigenvalue problem

            Kx = λMx,

        the iteration is performed on M^{-1}K.

        If Pfun is provided, it is used to apply M^{-1}.
        If Pfun is None and Bfun is provided, the system

            Mz = y

        is solved with SciPy CG without preconditioner.
        """

        rng = np.random.default_rng(seed)

        x = rng.random(N)
        x /= np.linalg.norm(x)

        lam_old = 0.0

        # --------------------------------------------------------
        # Opérateur de masse pour le CG sans préconditionneur
        # --------------------------------------------------------

        if Bfun is not None:

            M_operator = scsplin.LinearOperator(
                shape=(N, N),
                matvec=Bfun,
                dtype=float,
            )

        else:

            M_operator = None

        # --------------------------------------------------------
        # POWER ITERATION
        # --------------------------------------------------------

        for it in range(maxiter):


            y = Afun(x)

            if Pfun is not None:

                z = Pfun(y)

            elif Bfun is not None:

                # Aucun préconditionneur :
                # résolution de Mz = y avec CG SciPy
                z, info = scsplin.cg(
                    M_operator,
                    y,
                    rtol=1e-12,
                    atol=0.0,
                    maxiter=10000,
                )

                if info != 0:
                    raise RuntimeError(
                        f"CG non convergé pour Mz = y, info={info}"
                    )

            else:

                z = y

            norm_z = np.linalg.norm(z)

            if norm_z == 0.0:
                raise ValueError("Vecteur nul.")

            x = z / norm_z

            # ----------------------------------------------------
            # Quotient de Rayleigh
            # ----------------------------------------------------

            Kx = Afun(x)

            if Bfun is not None:

                Mx = Bfun(x)

                denominator = np.dot(x, Mx)

                if denominator == 0.0:
                    raise ZeroDivisionError(
                        "Le dénominateur x^T M x est nul."
                    )

                lam = (
                    np.dot(x, Kx)
                    / denominator
                )

            else:

                lam = np.dot(x, Kx)


            if abs(lam - lam_old) < (
                tol * max(1.0, abs(lam))
            ):
                break

            lam_old = lam

        return (
            np.array([lam]),
            x.reshape(-1, 1),
        )


    









