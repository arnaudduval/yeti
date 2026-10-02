"""
Matrix-free iterative solve for mechanics problems: K @ u = F, K never
assembled.

Wraps any `apply(v) -> K @ v` operator over the FULL dof vector (e.g. a
`WQMatrixFreeStiffness` instance's `.apply`) into a
`scipy.sparse.linalg.LinearOperator` restricted to the FREE dofs, and solves
via an iterative Krylov method -- GMRES by default (see the caveat below for
why, and why not CG despite K being SPD).

The BC handling is the standard matrix-free "zero-padding" trick: a
free-dof vector is extended to a full dof vector by placing zeros at the
fixed dofs, `apply()` is called on that, and the result is restricted back
to the free dofs. For a HOMOGENEOUS Dirichlet BC (fixed value 0, e.g. a
clamped/blocked support) this directly gives K_ff @ v_free -- exactly the
reduced operator CG needs, without ever slicing rows/columns out of an
assembled matrix (there is no assembled matrix). A non-zero prescribed
value at the fixed dofs only costs one extra `apply()` call, to move the
K_fc @ u_c term to the right-hand side once, before iterating.

**Known caveat -- `WQMatrixFreeStiffness.apply()` is only approximately
symmetric on curved geometry.** Weighted quadrature's "gather" matrices
(`B0`/`B1`) and "scatter" matrices (`W00..W11`) are fit independently (see
`WeightedQuadrature1D`'s docstring) and do not form an exact adjoint pair;
for an AFFINE patch this doesn't matter (the pulled-back material tensor is
constant and the construction stays exactly symmetric -- confirmed to
machine precision in `tests/future/test_wq_stiffness.py`), but for a CURVED
patch the tensor varies per quadrature point and sits between an
independently-fit gather/scatter pair, breaking exact symmetry by a margin
that shrinks under mesh refinement (same convergence as the K-approximation
error itself -- observed empirically from ~5% relative asymmetry on a
single-element, degree-4 quarter-ring down to ~1e-5 after 2 uniform
refinements). Plain CG can fail to converge outright on the worst (coarse,
high-degree, strongly curved) cases -- confirmed to hit `maxiter` with zero
progress on a 1-element, degree-4 quarter-ring, where GMRES converges
cleanly in ~200 iterations.

This is not a `future`-specific defect: `pymfiga`'s own reference
implementation has the identical structural asymmetry (confirmed
numerically -- its matrix-free *mass* operator, same
gather/multiply/scatter shape, is comparably asymmetric on the same curved
case), but its `LinearSolver` (`common/numerics/solvers/linear_solver.py`)
defaults to `linear_type="gmres"`, not `"cg"` (`PhysicsArgs.linear_solver_type
= "gmres"`, `common/physics/core.py`), and always applies a spatial
preconditioner in `MechanicalModel.solve_linearized_system`. pymfiga simply
never exercises plain, unpreconditioned CG on this operator -- so this
module defaults to `"gmres"` too, matching that validated choice, rather
than defaulting to CG's speed and leaving the caller to discover the
failure mode. Pass `method="cg"` explicitly when the problem is known to be
affine/well-refined enough for it to be safe and faster.

Iteration count still grows with mesh refinement even once convergence
itself isn't in question (unpreconditioned GMRES: ~400 iterations at ~100
dofs, ~9000+ at ~800 dofs on a quarter-ring). Pass `precondition_fn` (e.g.
`FastDiagonalizationPreconditioner.apply`, see `fast_diagonalization.py` --
ported from pymfiga's own `SingleFastDiagonalization`, the preconditioner
pymfiga's static elasticity solve always applies) to flatten that growth:
confirmed to cut iteration counts 4x-80x, growing with mesh size, keeping
the *preconditioned* count roughly constant (~100-120) across refinement
levels where the unpreconditioned count explodes into the thousands.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import LinearOperator, bicgstab, cg, gmres

_METHODS = {"cg": cg, "gmres": gmres, "bicgstab": bicgstab}


def free_dof_operator(apply_fn, ndof, fixed_dofs):
    """
    Build a `LinearOperator` acting on the FREE-dof subspace only, from any
    `apply(v) -> K @ v` matrix-free operator over the FULL dof vector.

    Parameters
    ----------
    apply_fn : callable(np.ndarray) -> np.ndarray
        Matrix-free K @ v, e.g. a `WQMatrixFreeStiffness` instance's `.apply`.
    ndof : int
        Total (full) number of dofs `apply_fn` expects/returns.
    fixed_dofs : array_like of int
        Global dof indices held at a prescribed (Dirichlet) value.

    Returns
    -------
    op : scipy.sparse.linalg.LinearOperator
        Shape `(n_free, n_free)`, symmetric positive-definite whenever the
        full K is.
    free_dofs : np.ndarray of int
        The complementary (unconstrained) dof indices, sorted -- `op`'s
        rows/columns are in this order.
    """
    fixed_dofs = np.asarray(fixed_dofs, dtype=np.int64)
    free_dofs = np.setdiff1d(np.arange(ndof), fixed_dofs, assume_unique=False)
    n_free = free_dofs.size

    def matvec(v_free):
        v_full = np.zeros(ndof)
        v_full[free_dofs] = v_free
        return apply_fn(v_full)[free_dofs]

    op = LinearOperator((n_free, n_free), matvec=matvec, dtype=np.float64)
    return op, free_dofs


def solve(apply_fn, ndof, F, fixed_dofs, prescribed=None, method="gmres",
          precondition_fn=None, rtol=1e-8, atol=0.0, maxiter=None, callback=None,
          **solver_kwargs):
    """
    Solve `K @ u = F` for a mechanics problem with Dirichlet BCs, matrix-free
    (K is never assembled), via an iterative Krylov method.

    Parameters
    ----------
    apply_fn : callable(np.ndarray) -> np.ndarray
        Matrix-free K @ v over the FULL dof vector (length `ndof`).
    ndof : int
        Total number of dofs.
    F : np.ndarray, shape (ndof,)
        Global load vector.
    fixed_dofs : array_like of int
        Global dof indices with a prescribed (Dirichlet) value.
    prescribed : array_like of float, optional
        Prescribed value at each of `fixed_dofs`, same order/length (default:
        all zero -- the common clamped/blocked-support case).
    method : {"cg", "gmres", "bicgstab"}, default "gmres"
        Which `scipy.sparse.linalg` solver to use. Defaults to "gmres",
        matching pymfiga's own default -- see the module docstring's caveat
        about `WQMatrixFreeStiffness.apply()` only being approximately
        symmetric on curved geometry, which can make "cg" fail outright on
        coarse/high-degree/strongly-curved patches. Pass "cg" explicitly for
        the speed win on problems known not to hit that (affine geometry, or
        well-refined curved meshes).
    precondition_fn : callable(np.ndarray) -> np.ndarray, optional
        Approximate K^{-1} @ v over the FULL dof vector, same convention as
        `apply_fn` -- e.g. `FastDiagonalizationPreconditioner(...).apply`
        (`fast_diagonalization.py`). Cuts iteration counts substantially,
        and keeps them roughly constant under mesh refinement instead of
        growing -- see the module docstring. `None` (default): unpreconditioned.
    rtol, atol, maxiter, callback :
        Passed straight through to the chosen `scipy.sparse.linalg` solver.
    **solver_kwargs :
        Any other solver-specific keyword, passed straight through -- e.g.
        `restart=` for `method="gmres"`. scipy's own GMRES default
        (`restart=20`) discards Krylov information every 20 iterations,
        which can meaningfully inflate the iteration count on a
        preconditioned system that would otherwise converge in fewer steps
        than that (confirmed on the quarter-ring elasticity benchmark:
        `restart=20` needed ~55 iterations where `restart=50`+ needed only
        ~49, unchanged for any larger restart tried) -- pass a larger
        `restart` (or `len(free_dofs)` for effectively unrestarted GMRES) if
        iteration count matters more than per-iteration memory.

    Returns
    -------
    u : np.ndarray, shape (ndof,)
        Full displacement field: free dofs solved, fixed dofs set to
        `prescribed`.
    info : int
        The solver's convergence flag (0 = converged to the requested
        tolerance within `maxiter`).
    """
    if method not in _METHODS:
        raise ValueError(f"method must be one of {sorted(_METHODS)}, got {method!r}")
    solver = _METHODS[method]

    if method == "gmres" and callback is not None:
        # Silence scipy's DeprecationWarning about a future default change --
        # "legacy" is today's actual default (callback receives the current
        # iterate xk), and every callback in this codebase either ignores its
        # argument (iteration counters) or explicitly requests
        # callback_type="pr_norm" itself via **solver_kwargs, which this
        # setdefault leaves untouched.
        solver_kwargs.setdefault("callback_type", "legacy")

    fixed_dofs = np.asarray(fixed_dofs, dtype=np.int64)
    prescribed = np.zeros(fixed_dofs.shape[0]) if prescribed is None \
        else np.asarray(prescribed, dtype=np.float64)

    op, free_dofs = free_dof_operator(apply_fn, ndof, fixed_dofs)
    M = None
    if precondition_fn is not None:
        M, _ = free_dof_operator(precondition_fn, ndof, fixed_dofs)

    # K @ [0 at free ; prescribed at fixed], restricted to free dofs, is
    # exactly K_fc @ u_c -- move it to the RHS once before iterating (a
    # no-op extra apply() when prescribed is all zero).
    u_c_full = np.zeros(ndof)
    u_c_full[fixed_dofs] = prescribed
    rhs = F[free_dofs] - apply_fn(u_c_full)[free_dofs]

    u_free, info = solver(op, rhs, rtol=rtol, atol=atol, maxiter=maxiter, M=M, callback=callback,
                         **solver_kwargs)

    u = np.zeros(ndof)
    u[fixed_dofs] = prescribed
    u[free_dofs] = u_free
    return u, info
