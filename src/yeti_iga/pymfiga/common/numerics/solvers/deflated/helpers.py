from typing import Callable
import numpy as np

# ---------------------------------------------------------------------------
# gmres1  –  plain preconditioned GMRES (priming cycle)
# ---------------------------------------------------------------------------


def gmres1(
    A_apply: Callable,
    x: np.ndarray,
    r: np.ndarray,
    m: int,
    M_apply: Callable,
    tol: float,
):
    """
    Run up to m steps of preconditioned GMRES starting from (x, r).
    Generates  A V[:, :m] = V[:, :m+1] H.

    Returns  x, r, V, H, k, resvec
    """
    dtype = np.result_type(r, np.complex128)
    n = len(r)
    V = np.zeros((n, m + 1), dtype=dtype)
    H = np.zeros((m + 1, m), dtype=dtype)
    resvec = np.zeros(m)

    beta = np.linalg.norm(r)
    V[:, 0] = r / beta
    y = np.zeros(1)
    res = np.array([beta])

    for k in range(1, m + 1):
        w = M_apply(A_apply(V[:, k - 1]))

        for j in range(k):  # classical Gram-Schmidt
            H[j, k - 1] = V[:, j].conj().dot(w)
            w -= H[j, k - 1] * V[:, j]

        H[k, k - 1] = np.linalg.norm(w)
        V[:, k] = w / H[k, k - 1]

        rhs = np.zeros(k + 1)
        rhs[0] = beta
        y, *_ = np.linalg.lstsq(H[: k + 1, :k], rhs, rcond=None)
        res = rhs - H[: k + 1, :k] @ y
        resvec[k - 1] = np.linalg.norm(res)

        if resvec[k - 1] < tol:
            x = x + V[:, :k] @ y
            r = V[:, : k + 1] @ res
            return x, r, V[:, : k + 1], H[: k + 1], k, resvec[:k]

    x = x + V[:, :m] @ y
    r = V[:, : m + 1] @ res
    return x, r, V, H, m, resvec


# ---------------------------------------------------------------------------
# gmres2  –  deflated Arnoldi for the main GCRO-DR cycle
# ---------------------------------------------------------------------------


def gmres2(
    A_apply: Callable,
    r: np.ndarray,
    m: int,
    M_apply: Callable,
    C: np.ndarray,
    tol: float,
):
    """
    Run up to m deflated Arnoldi steps.
    Generates  (I - C Cᵀ) M⁻¹ A V[:, :m] = V[:, :m+1] H.

    B[:, j] = Cᵀ (M⁻¹ A v_j)  is the cross-term coupling C to V.

    Returns  V, H, B, k, resvec
    """
    dtype = np.result_type(r, np.complex128)
    n = len(r)
    k_c = C.shape[1]
    V = np.zeros((n, m + 1), dtype=dtype)
    H = np.zeros((m + 1, m), dtype=dtype)
    B = np.zeros((k_c, m), dtype=dtype)
    resvec = np.zeros(m)

    beta = np.linalg.norm(r)
    V[:, 0] = r / beta
    y = np.zeros(1)
    res = np.array([beta])

    for k in range(1, m + 1):
        w = M_apply(A_apply(V[:, k - 1]))

        B[:, k - 1] = C.conj().T @ w  # record C-component before deflation
        w = w - C @ B[:, k - 1]  # deflate

        for j in range(k):  # classical Gram-Schmidt
            H[j, k - 1] = V[:, j].conj().dot(w)
            w -= H[j, k - 1] * V[:, j]

        H[k, k - 1] = np.linalg.norm(w)
        V[:, k] = w / H[k, k - 1]

        rhs = np.zeros(k + 1)
        rhs[0] = beta
        y, *_ = np.linalg.lstsq(H[: k + 1, :k], rhs, rcond=None)
        res = rhs - H[: k + 1, :k] @ y
        resvec[k - 1] = np.linalg.norm(res)

        if resvec[k - 1] < tol:
            return V[:, : k + 1], H[: k + 1], B[:, :k], k, resvec[:k]

    return V, H, B, m, resvec


# ---------------------------------------------------------------------------
# getHarmVecs1  –  harmonic Ritz extraction after a plain GMRES run
# ---------------------------------------------------------------------------


def getHarmVecs1(m: int, k: int, H: np.ndarray):
    """
    Extract k harmonic Ritz vectors from the (m+1)xm Hessenberg H.

    Harmonic Ritz matrix:
        G = Hm  +  h_{m+1,m}²  (Hm^{-T} eₘ) eₘᵀ

    Returns the k eigenvectors of G with *smallest* |eigenvalue|.

    Parameters
    ----------
    m : number of GMRES steps  (H is (m+1)xm)
    k : number of vectors to return
    H : (m+1, m) upper Hessenberg

    Returns  harmVecs : (m, k)
    """
    Hm = H[:m, :m]  # m × m  (square upper Hessenberg)

    em = np.zeros(m)
    em[-1] = 1.0
    Hm_invT_em = np.linalg.solve(Hm.conj().T, em)  # Hm^{-T} eₘ

    harmRitzMat = Hm + H[m, m - 1] ** 2 * np.outer(Hm_invT_em, em)

    eigvals, eigvecs = np.linalg.eig(harmRitzMat)
    order = np.argsort(np.abs(eigvals))  # ascending |λ|
    harmVecs = eigvecs[:, order[:k]]
    return harmVecs


# ---------------------------------------------------------------------------
# getHarmVecs2  –  harmonic Ritz extraction inside the main GCRO-DR loop
# ---------------------------------------------------------------------------


def getHarmVecs2(
    m: int, k: int, H2: np.ndarray, V: np.ndarray, U: np.ndarray, C: np.ndarray
):
    """
    Extract k harmonic Ritz vectors from the augmented (m+1)xm Hessenberg H2.

    Parameters
    ----------
    m  : p + k_recycle  (total joint subspace dimension)
    k  : number of vectors to extract
    H2 : (m+1, m) augmented Hessenberg
    V  : (n, p+1) Krylov basis including the (p+1)-th column
    U  : (n, k_recycle) recycled subspace (column-normalised)
    C  : (n, k_recycle) A-orthonormal basis, C = A U

    Returns  harmVecs : (m, k)
    """
    k_c = U.shape[1]
    p = m - k_c

    B_eig = H2.conj().T @ H2  # m × m

    # Phi is (m+1) × m, encoding the GCRO-DR Petrov-Galerkin condition
    #
    #          [ Cᵀ U      0   ]   k_c rows
    # Phi  =   [ Vᵀ U      I_p ]   (p+1) rows
    #
    Phi = np.zeros((m + 1, m), dtype=np.complex128)
    Phi[:k_c, :k_c] = C.conj().T @ U
    Phi[k_c:, :k_c] = V.conj().T @ U  # V is n × (p+1)
    Phi[k_c : k_c + p, k_c:m] = np.eye(p)

    A_eig = H2.conj().T @ Phi

    harmVals, harmVecs = np.linalg.eig(np.linalg.solve(B_eig, A_eig))
    iperm = np.argsort(np.abs(harmVals))[::-1]  # ascending |λ|
    return harmVecs[:, iperm[:k]]
