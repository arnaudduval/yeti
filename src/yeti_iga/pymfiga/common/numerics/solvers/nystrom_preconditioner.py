from yeti_iga.pymfiga.common.base.enum import Constants
from scipy import sparse as sp
from numpy.linalg import norm, svd
from scipy.linalg import cholesky, solve_triangular, qr
from typing import Callable, Optional
import numpy as np


def gaussian_test_matrix(n: int, sample_size: int, orthonormal=False) -> np.ndarray:
    Omega = (1 / np.sqrt(n)) * np.random.randn(n, sample_size)
    if orthonormal:
        Omega = qr(Omega, mode="economic")[0]
    return np.array(Omega)


def rangefinder(A: Callable, Omega: np.ndarray, orthogonalize=False) -> np.ndarray:
    Y = np.column_stack([A(Omega[:, i]) for i in range(Omega.shape[1])])
    if orthogonalize:
        Y = qr(Y, mode="economic")[0]
    return np.array(Y)


def nystrom_sketch(A: Callable, n: int, rank: int, sample_size: int):

    for param, name in zip([n, rank, sample_size], ["n", "rank", "sample_size"]):
        if not isinstance(param, int):
            raise TypeError(f"{name} must be an integer")
        if param <= 0:
            raise ValueError(f"{name} must be a positive integer")

    if rank > sample_size:
        raise ValueError("k must be less than r")

    # Set random vectors
    Omega = gaussian_test_matrix(n, sample_size, orthonormal=True)

    # Set the stability shift
    nu = Constants.TINY * min(1.0, norm(Omega))
    # TODO: get better estimate of nu

    # Regularize
    A_reg = lambda v: A(v) + nu * v

    Y = rangefinder(A_reg, Omega)

    Z = Omega.T @ Y

    R = cholesky(Z, lower=False)

    # Solve Y / R safely
    B = solve_triangular(R.T, Y.T, lower=True).T

    U, S, _ = svd(B, full_matrices=False)

    Lambda = np.maximum(0, S**2 - nu)

    return U[:, :rank], Lambda[:rank]


class NystromSchurPreconditioner:
    """
    Preconditioner for the Schur-complement-like SPD system

        M = mu * I_m  +  C @ A^{-1} @ C.T          (m x m)

    where

        A^{-1} = U @ diag(Qinv) @ U.T                     (n x n)

    and U, U.T are accessible only as matvec callables.

    Two modes
    ---------
    Dense (exact)
        Built by applying M to m unit vectors.  Used automatically when
        m < small_threshold or rank >= m.  Gives exact Cholesky factorisation.

    Nyström (approximate)
        Sketches  S = C @ A^{-1} @ C.T  (the low-rank part, without mu * I)
        using `rank` random vectors.  The full preconditioner is then recovered
        via the Woodbury identity with shift.

        Sketching S rather than M is intentional: S has faster spectral decay
        than M (no identity floor), so fewer sketch vectors are needed to
        resolve the spectrum.

    Parameters
    ----------
    C : sparse matrix, shape (m, n)
        Fixed constraint / coupling matrix.
    U : callable
        v -> U @ v  (n,) -> (n,).  Matvec of the eigenvector matrix.
    UT : callable
        v -> U.T @ v  (n,) -> (n,).  Transpose matvec.
    Qinv : array-like (n,)
        Diagonal of Q^{-1} (reciprocal eigenvalues of A)
    mu : float
        Spectral shift parameter.
    rank : int
        Number of sketch vectors for Nyström mode.
        Use rank >= m to force dense mode regardless of small_threshold.
    small_threshold : int
        Switch to dense mode when m < small_threshold.
    random_state : int or None
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        C: sp.csr_array,
        U: Callable,
        UT: Callable,
        Qinv: np.ndarray,
        mu: float,
        rank: int = 24,
        small_threshold: int = 576,
        random_state: Optional[int] = None,
    ):
        self.C = C
        self.U = U
        self.UT = UT
        self.Qinv = Qinv
        self.sqrtQinv = np.sqrt(Qinv)
        self.mu = mu
        self.rank = rank
        self.small_threshold = small_threshold
        self.rng = np.random.default_rng(random_state)

        self._built = False
        self._dense = False

        # Dense mode
        self._L = np.array([])

        # Nyström mode
        self._U_S = np.array([])
        self._Lambda_S = np.array([])

    # ------------------------------------------------------------------
    # Internal matvecs
    # ------------------------------------------------------------------

    def _apply_Ainv(self, v: np.ndarray) -> np.ndarray:
        """Apply A^{-1} = U Q^{-1} U^T to v."""
        return self.U(self.Qinv * self.UT(v))

    def _apply_S(self, v: np.ndarray):
        """Apply S = C A^{-1} C^T to v  (the low-rank part of M)."""
        return self.C @ self._apply_Ainv(self.C.T @ v)

    # ------------------------------------------------------------------
    # Dense build  (exact Cholesky of M)
    # ------------------------------------------------------------------

    def _build_dense(self):
        CC = self.C
        m = CC.shape[0]
        mu = self.mu

        SS = np.zeros((m, m))
        for i in range(m):
            cT = CC[i, :].todense()
            invPcT = self._apply_Ainv(cT)
            SS[:, i] = CC @ invPcT

        # Symmetrise against floating-point noise
        SS = 0.5 * (SS + SS.T)
        SS += mu * np.eye(m)
        self._L: np.ndarray = cholesky(SS, lower=True)
        self._dense = True
        self._built = True

    # ------------------------------------------------------------------
    # Nyström build
    # ------------------------------------------------------------------

    def _build_nystrom(self, oversampling: int = 10):
        m = self.C.shape[0]
        r = min(self.rank, m)
        k = min(self.rank + oversampling, m)

        UU, LambdaS = nystrom_sketch(self._apply_S, m, r, k)

        self._U_S = UU  # approx eigenvectors of S
        self._Lambda_S = LambdaS  # approx eigenvalues of S
        self._dense = False
        self._built = True

    # ------------------------------------------------------------------
    # Public build
    # ------------------------------------------------------------------

    def build(self):
        """Build the preconditioner (dense or Nyström based on m and rank)."""
        m = self.C.shape[0]
        if m < self.small_threshold or self.rank >= m:
            self._build_dense()
        else:
            self._build_nystrom()

    # ------------------------------------------------------------------
    # Apply  P^{-1}
    # ------------------------------------------------------------------

    def __call__(self, x: np.ndarray) -> np.ndarray:
        """
        Apply P⁻¹ to vector x.

        Dense mode: two triangular solves — exact inverse of M.

        Nyström mode: implements idea from Frangella, Tropp & Udell (2023),
        adapted to the Schur complement M = μI + S where S ≈ U Λ Uᵀ.

            P⁻¹ x = (λ_k + μ) · U (Λ + μI)⁻¹ Uᵀ x  +  (I - U Uᵀ) x

        where λ_k = LambdaS[-1] is the smallest retained eigenvalue of S.
        This is the paper's formula with Ŝ_nys = U Λ Uᵀ playing the role
        of the low-rank approximation to S.

        Quality degrades when μ >> λ_max(S) (identity dominates, sketch irrelevant)
        or when rank is too small to capture the spectrum of S.
        """
        if not self._built:
            raise RuntimeError("Call build() before applying the preconditioner.")

        if self._dense:
            y = solve_triangular(self._L, x, lower=True)
            return solve_triangular(self._L.T, y, lower=False)

        # Nyström / Woodbury
        U = self._U_S
        Lambda = self._Lambda_S
        lam = Lambda[-1]
        mu = self.mu

        cache = U.T @ x
        cache *= (lam + mu) / (Lambda + mu) - 1
        y = U @ cache
        return x + y

    @property
    def rank_used(self):
        if not self._built or self._dense:
            return None
        return len(self._Lambda_S)

    def __repr__(self):
        if not self._built:
            return "NystromSchurPreconditioner(not built)"
        m = self.C.shape[0]
        if self._dense:
            return f"NystromSchurPreconditioner(dense, m={m}, mu={self.mu:.2e})"
        ratio = self._Lambda_S[-1] / (self._Lambda_S[0] + Constants.SAFEGUARD)
        return (
            f"NystromSchurPreconditioner("
            f"m={m}, mu={self.mu:.2e}, "
            f"rank={self.rank_used}, "
            f"spectral_ratio={ratio:.2e})"
        )
