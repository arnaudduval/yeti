"""
Fast-diagonalization (FD) preconditioner for WQMatrixFreeStiffness, ported
from pymfiga's iga/fastdiagonalization/single/cspace.py
(SpaceFD/SingleFastDiagonalization) -- the preconditioner pymfiga's own
MechanicalModel.solve_linearized_system() always applies for a static
elasticity solve (as opposed to "scaled_mass", which is a *different*,
mass-specific preconditioner used only by pymfiga's explicit-dynamics
models -- not applicable here, since there is no matrix-free mass action in
`future` yet).

The idea (classical "fast diagonalization method", e.g. Lynch/Rice/Thomas
generalized to IGA via Sangalli & Tani): approximate K, for ONE displacement
component at a time (block-diagonal in the components -- the coupling terms
between u_x/u_y are dropped for preconditioning purposes only, not for the
actual operator), by a Kronecker SUM `K_I ~= K_u (x) M_v + M_u (x) K_v` built
from the *parametric* 1D weighted-quadrature mass/stiffness pencils. A
Kronecker sum of this form is diagonalized exactly by the *generalized*
eigenvectors of the 1D pencils: solving K_d @ v = lambda * M_d @ v per
direction gives eigenvector matrices V_u, V_v such that

    P^{-1} @ r = V_u @ ((V_u^T @ R @ V_v) / Lambda) @ V_v^T

applies the exact inverse of that Kronecker-sum approximation in O(n^1.5)
work (dense 1D eigendecompositions, size = ctrlpts per direction, done once)
instead of ever forming or factorizing an n x n matrix.

`method=` controls how much geometry/material information is folded into
that Kronecker sum before diagonalizing it -- matching the four variants
implemented in `iga_wq_mf/fsrc/mf_wq_solver.f90` and
`tensor_algebra.f90` (the actual Fortran code behind Cornejo Fuentes'
thesis, vendored in this repo at `src/yeti_iga/iga_wq_mf`), from weakest to
strongest:

- `"classic"` (the default, no `law` needed): no geometry or material
  information at all -- `Lambda[iv, iu] = eigval_v[iv] + eigval_u[iu]`. The
  thesis's Eq. (2.13)/(2.14).
- `"mean"` (needs `law=`): `Lambda[iv, iu] = c[I,1]*eigval_v[iv] +
  c[I,0]*eigval_u[iu]`, with `c[I, l]` a single geometry+material-aware
  scalar per direction and displacement component -- the thesis's own
  printed Eq. (2.15) formula (`compute_stiffness_correctors()`).
- `"tensor"` (needs `law=`): Montardini's *per-quadrature-point* alternating
  min-max tensor decomposition ("Preconditioners for Isogeometric Analysis"
  -- the Fortran code's "TD" method, `compute_tensor_decomposition_correctors()`
  and `tensor_decomposition_3d`'s direct 2D port), replacing the single mean
  coefficient with a full coefficient *vector* per direction (one value per
  WQ quadrature point) chosen to best separably approximate the true,
  spatially-varying coefficient field -- a strictly better approximation
  than a single mean, and the method the Fortran code's naming suggests
  Table 2.3's ~39-iteration numbers actually come from (not the simpler
  "JM"/mean variant the printed thesis text alone describes).

`diag_physical=` (any method) layers the Fortran code's "S" suffix
(TDS/JMS) on top: a symmetric point-wise rescaling by
`sqrt(diag(P)/diag(K))` before and after applying the base FD inverse,
intended (per `scaling_FastDiag` in `tensor_algebra.f90`) to correct the
preconditioner's diagonal to match the true operator's more closely.
**Status: implemented and unit-tested (the scale factors it computes are
sane -- close to 1, no blow-ups -- and `compute_operator_diagonal()`'s
Gauss-based `diag_physical` agrees with an independent WQ dense-probe to
~1%), but empirically it makes GMRES need MORE iterations, not fewer, on
every case tried here (e.g. quarter-ring, degree 4, `mean`: 57 -> 105 with
scaling added; `tensor`: 51 -> 60) -- the opposite of the Fortran code's own
"TDS"/"JMS" framing.** The root cause hasn't been pinned down (plausibly a
subtlety in how the scaling should interact with the free/fixed-DOF
restriction, or with the eigenbasis this class's `apply()` already
projects into, that the Fortran's own CSR-based construction sidesteps
differently) -- treat `diag_physical=` as experimental/not recommended
until this is resolved, not as a working feature.

Dirichlet boundary conditions are handled the same way pymfiga does: any
direction/side that is *entirely* fixed for a given component (i.e. every
control point on that edge has that component in `fixed_dofs`) has its
corresponding row/column dropped from the 1D pencil before the
eigendecomposition -- this only supports edge-uniform Dirichlet BCs (the
same scope pymfiga's own BoundaryCondition API has), not arbitrary fixed-dof
sets.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as sclin

from .bspline import WeightedQuadrature1D

# 1D composite trapezoidal rule on [0, 1] with 3 points {0, 0.5, 1} (h=0.5):
# integral(f) ~= 0.25*f(0) + 0.5*f(0.5) + 0.25*f(1). Tensor-producted across
# parametric directions below, this is exactly the thesis's own "3^d grid,
# weights 8/4/2/1 (interior/face/edge/corner)" rule (Eq. (2.15)'s own
# footnote) written per-point instead of per-category -- same rule, same
# result, just without needing to classify each of the 3^d points by how
# many of its coordinates land on a domain boundary. Used only by
# compute_stiffness_correctors() (the "mean" method).
_TRAPEZOID_POINTS = (0.0, 0.5, 1.0)
_TRAPEZOID_WEIGHTS = (0.25, 0.5, 0.25)


def _extract_lame(law):
    """(lambda, mu) from any isotropic ConstitutiveLaw, by numerically
    probing stiffness_density() with a unit vector -- same trick as
    WQMatrixFreeStiffness.cpp's own extract_lame(), so this always matches
    whatever tensor convention the stiffness operator itself actually uses."""
    e0 = np.array([1.0, 0.0])
    x0 = np.array([0.0, 0.0])
    K = np.asarray(law.stiffness_density(e0, e0, x0))
    mu = K[1, 1]
    lam = K[0, 0] - 2.0 * mu
    return lam, mu


def _evaluate_jacobian(patch, xi, eta):
    """Physical Jacobian (J11, J12, J21, J22) = (dx/du, dx/dv, dy/du, dy/dv)
    at one parametric point (xi, eta) -- NURBS-aware (quotient rule), same
    convention as WQMatrixFreeStiffness's constructor. Pure Python, built
    from already-exposed primitives (1D basis derivatives + control point
    coordinates/weights) since evaluating the Jacobian at one arbitrary
    point isn't otherwise exposed."""
    bsp_u, bsp_v = patch.tensor.components
    deg_u, deg_v = bsp_u.degree, bsp_v.degree
    span_u = bsp_u.find_span(xi)
    span_v = bsp_v.find_span(eta)
    Ndu = bsp_u.basis_funs_derivatives(span_u, xi, 1)
    Ndv = bsp_v.basis_funs_derivatives(span_v, eta, 1)
    Nu, dNu = Ndu[0], Ndu[1]
    Nv, dNv = Ndv[0], Ndv[1]

    rational = patch.cp_manager.is_rational
    weights_view = patch.cp_manager.weights_view() if rational else None
    nu = patch.local_shape[0]

    W = dWdu = dWdv = 0.0
    Xs = Ys = dXdu = dXdv = dYdu = dYdv = 0.0
    for b in range(deg_v + 1):
        iv = span_v - deg_v + b
        for a in range(deg_u + 1):
            iu = span_u - deg_u + a
            idx = iv * nu + iu
            P = patch.control_point(idx)
            w = weights_view[patch.global_indices[idx]] if rational else 1.0
            Nval, dNu_val, dNv_val = Nu[a] * Nv[b], dNu[a] * Nv[b], Nu[a] * dNv[b]
            W += w * Nval
            dWdu += w * dNu_val
            dWdv += w * dNv_val
            Xs += w * Nval * P[0]
            Ys += w * Nval * P[1]
            dXdu += w * dNu_val * P[0]
            dXdv += w * dNv_val * P[0]
            dYdu += w * dNu_val * P[1]
            dYdv += w * dNv_val * P[1]

    if rational:
        W2 = W * W
        J11 = (dXdu * W - Xs * dWdu) / W2
        J12 = (dXdv * W - Xs * dWdv) / W2
        J21 = (dYdu * W - Ys * dWdu) / W2
        J22 = (dYdv * W - Ys * dWdv) / W2
    else:
        J11, J12, J21, J22 = dXdu, dXdv, dYdu, dYdv

    return J11, J12, J21, J22


def _diag_tensor_at_point(lam, mu, patch, xi, eta):
    """D[m, l] = D_ll^(m,m)(xi, eta) -- the diagonal block (displacement
    component m) of the pulled-back elasticity tensor WQMatrixFreeStiffness
    itself builds (see its .cpp: stiffness_property_[I][J][a][b], same
    lambda/mu convention), restricted to I=J=m, a=b=l."""
    J11, J12, J21, J22 = _evaluate_jacobian(patch, xi, eta)
    detJ = J11 * J22 - J12 * J21
    # invJ[a][l] = d(param_a)/d(phys_l), matching WQMatrixFreeStiffness's own convention.
    invJ = np.array([[J22, -J12], [-J21, J11]]) / detJ
    absDetJ = abs(detJ)
    D = np.zeros((2, 2))
    for m in range(2):
        for l in range(2):
            D[m, l] = absDetJ * (
                (lam + mu) * invJ[l, m] ** 2 + mu * (invJ[l, 0] ** 2 + invJ[l, 1] ** 2)
            )
    return D


def compute_stiffness_correctors(patch, law):
    """
    Geometry- and material-aware scalar coefficients `c[m, l]` for the
    fast-diagonalization preconditioner -- Cornejo Fuentes' thesis, Section
    2.3.2 / Eq. (2.15) extended to elasticity ("mean" method below).
    `c[m, l]` approximates `integral_{domain} D_ll^(m,m)(xi) dOmega_hat` (see
    `_diag_tensor_at_point`), via the thesis's own cheap 3^d-point (here 2D:
    3x3 = 9 points) composite-trapezoidal rule at parametric points
    {0, 0.5, 1} per direction.

    A single mean coefficient per (component, direction) -- see
    `compute_tensor_decomposition_correctors()` for the more faithful,
    per-quadrature-point variant ("tensor"/"TD") the Fortran code behind the
    thesis's own reported numbers actually uses.

    Returns
    -------
    c : ndarray, shape (2, 2). `c[m, l]`: displacement component `m`,
        parametric direction `l` (0=u, 1=v).
    """
    lam, mu = _extract_lame(law)
    c = np.zeros((2, 2))
    for xi, wu in zip(_TRAPEZOID_POINTS, _TRAPEZOID_WEIGHTS):
        for eta, wv in zip(_TRAPEZOID_POINTS, _TRAPEZOID_WEIGHTS):
            c += (wu * wv) * _diag_tensor_at_point(lam, mu, patch, xi, eta)
    return c


def _tensor_decomposition_2d(CC0, CC1, n_iter=2):
    """
    Montardini's alternating min-max tensor decomposition -- direct,
    line-for-line port of `tensor_decomposition_2d` in
    `iga_wq_mf/fsrc/tensor_algebra.f90` (the actual Fortran code behind
    Cornejo Fuentes' thesis). `CC0[i1, i2]`, `CC1[i1, i2]`: the coefficient
    tensor's own diagonal entries `D_00`, `D_11`, sampled at the
    (WQ quadrature point) tensor grid. Returns per-quadrature-point
    coefficient vectors `(K_u, K_v, M_u, M_v)`, one value per quadrature
    point in each direction, used to rescale the WQ gather matrices' ROWS
    (the quadrature-point axis) before forming the univariate mass/stiffness
    pencils -- a strictly finer approximation than a single mean coefficient
    (`compute_stiffness_correctors`), since it lets the Kronecker-sum
    preconditioner track the true coefficient field's variation along each
    direction instead of averaging it away.
    """
    nqu, nqv = CC0.shape
    Mu, Mv = np.ones(nqu), np.ones(nqv)
    Ku, Kv = np.ones(nqu), np.ones(nqv)
    CC = (CC0, CC1)

    for _ in range(n_iter):
        # --- update K from M ---
        for k in range(2):
            V = np.empty((nqu, nqv))
            for i2 in range(nqv):
                for i1 in range(nqu):
                    UU = (Mu[i1], Mv[i2])
                    V[i1, i2] = CC[k][i1, i2] * UU[k] / (UU[0] * UU[1])
            if k == 0:
                for i1 in range(nqu):
                    Ku[i1] = np.sqrt(V[i1, :].min() * V[i1, :].max())
            else:
                for i2 in range(nqv):
                    Kv[i2] = np.sqrt(V[:, i2].min() * V[:, i2].max())

        # --- update M from K ---
        for k in range(2):
            l = 1 - k
            W = np.empty((nqu, nqv))
            for i2 in range(nqv):
                for i1 in range(nqu):
                    UU = (Mu[i1], Mv[i2])
                    WW = (Ku[i1], Kv[i2])
                    W[i1, i2] = CC[k][i1, i2] * UU[k] * UU[l] / (UU[0] * UU[1] * WW[k])
            if k == 0:
                for i1 in range(nqu):
                    Mu[i1] = np.sqrt(W[i1, :].min() * W[i1, :].max())
            else:
                for i2 in range(nqv):
                    Mv[i2] = np.sqrt(W[:, i2].min() * W[:, i2].max())

    return Ku, Kv, Mu, Mv


def compute_tensor_decomposition_correctors(patch, law, quad_u, quad_v, n_iter=2):
    """
    Per-quadrature-point coefficient vectors for each displacement
    component -- thesis Table 2.3's "TD" method (Montardini's tensor
    decomposition, `_tensor_decomposition_2d`), evaluated on the same
    diagonal elasticity tensor as `compute_stiffness_correctors` but sampled
    at the WQ quadrature points instead of a fixed 9-point grid, and kept as
    a full per-point vector instead of collapsed to a single mean.

    Parameters
    ----------
    quad_u, quad_v : array_like
        WQ quadrature point parametric positions in each direction (e.g.
        `WeightedQuadrature1D.build(...).quadpts`).

    Returns
    -------
    Kcoef, Mcoef : dict[int, tuple[ndarray, ndarray]]
        `Kcoef[m] = (Ku, Kv)`, `Mcoef[m] = (Mu, Mv)`, keyed by displacement
        component `m`; `Ku`/`Mu` have length `len(quad_u)`, `Kv`/`Mv` length
        `len(quad_v)`.
    """
    lam, mu = _extract_lame(law)
    nqu, nqv = len(quad_u), len(quad_v)
    D = np.zeros((2, 2, nqu, nqv))
    for i, xi in enumerate(quad_u):
        for j, eta in enumerate(quad_v):
            D[:, :, i, j] = _diag_tensor_at_point(lam, mu, patch, xi, eta)

    Kcoef, Mcoef = {}, {}
    for m in range(2):
        Ku, Kv, Mu, Mv = _tensor_decomposition_2d(D[m, 0], D[m, 1], n_iter=n_iter)
        Kcoef[m] = (Ku, Kv)
        Mcoef[m] = (Mu, Mv)
    return Kcoef, Mcoef


def compute_operator_diagonal(patch, law, basis_u=None, basis_v=None):
    """
    Diagonal of the assembled stiffness matrix (Gauss quadrature, via
    `PatchIntegrator`), for use as `diag_physical` in
    `FastDiagonalizationPreconditioner` -- the thesis's "S" (scaling)
    correction (`scaling_FastDiag` in `tensor_algebra.f90`) needs the TRUE
    operator's diagonal, which this reuses from `future`'s own trusted
    Gauss assembly rather than reimplementing a WQ-native sum-factorized
    diagonal extraction. Same interleaved DOF layout as `apply()`.
    """
    from .bspline import IGABasis1D, PatchIntegrator
    if basis_u is None:
        basis_u = IGABasis1D.build(patch.tensor.components[0], patch.tensor.components[0].degree + 1)
    if basis_v is None:
        basis_v = IGABasis1D.build(patch.tensor.components[1], patch.tensor.components[1].degree + 1)
    K = PatchIntegrator(patch, basis_u, basis_v, law).integrate_stiffness()
    return np.asarray(K.diagonal())


class FastDiagonalizationPreconditioner:
    """
    Fast-diagonalization preconditioner for a `WQMatrixFreeStiffness`
    problem. Build once per (patch, dof_manager, fixed_dofs); `apply(r)`
    approximates `K^{-1} @ r`, suitable as `M=` for
    `scipy.sparse.linalg.cg`/`gmres`/`bicgstab`, or as `Pfun` composed by
    hand (see `matrix_free_solver.solve`).

    `method`: `"classic"` (default, no geometry/material info), `"mean"`
    (needs `law=`, Eq. (2.15)'s single scalar per direction), or `"tensor"`
    (needs `law=`, Montardini's per-quadrature-point tensor decomposition --
    see the module docstring for the full comparison). `diag_physical=`
    (any method) layers the additional diagonal-scaling correction
    ("TDS"/"JMS" in the Fortran code) on top -- see
    `compute_operator_diagonal()`.
    """

    _METHODS = ("classic", "mean", "tensor")

    def __init__(self, patch, dof_manager, n_dofs_per_cp, fixed_dofs, quadtype="2",
                 law=None, method=None, diag_physical=None):
        if method is None:
            # Backwards-compatible default: law=None -> "classic" (unchanged
            # behavior); law= given -> "mean" (this class's only correction
            # before "tensor" existed). Pass method="tensor" explicitly for
            # the finer, per-quadrature-point Montardini decomposition.
            method = "classic" if law is None else "mean"
        if method not in self._METHODS:
            raise ValueError(f"method must be one of {self._METHODS}, got {method!r}")
        if method != "classic" and law is None:
            raise ValueError(f"method={method!r} requires law=")

        self.n_dofs_per_cp = n_dofs_per_cp
        self.nu = patch.local_shape[0]
        self.nv = patch.local_shape[1]
        self.nctrl = self.nu * self.nv

        wq_u = WeightedQuadrature1D.build(patch.tensor.components[0], quadtype)
        wq_v = WeightedQuadrature1D.build(patch.tensor.components[1], quadtype)
        Mu0 = (wq_u.W00 @ wq_u.B0).toarray()
        Ku0 = (wq_u.W10 @ wq_u.B1).toarray()
        Mv0 = (wq_v.W00 @ wq_v.B0).toarray()
        Kv0 = (wq_v.W10 @ wq_v.B1).toarray()

        mean_correctors = compute_stiffness_correctors(patch, law) if method == "mean" else None
        if method == "tensor":
            quad_u = np.asarray(wq_u.quadpts)
            quad_v = np.asarray(wq_v.quadpts)
            Kcoef, Mcoef = compute_tensor_decomposition_correctors(patch, law, quad_u, quad_v)

        fixed_dofs = set(int(d) for d in fixed_dofs)

        def edge_fully_fixed(direction, side, component):
            cps = patch.boundary_control_points(direction, side)
            return all(
                dof_manager.get_global_dof_indices(cp)[component] in fixed_dofs
                for cp in cps
            )

        diag_physical_full = None if diag_physical is None else np.asarray(diag_physical)

        self._free_u = []
        self._free_v = []
        self._Vu = []
        self._Vv = []
        self._lambda = []
        self._diag_scale = []

        for I in range(n_dofs_per_cp):
            free_u = np.arange(self.nu)
            if edge_fully_fixed(0, 0, I):
                free_u = free_u[1:]
            if edge_fully_fixed(0, 1, I):
                free_u = free_u[:-1]
            free_v = np.arange(self.nv)
            if edge_fully_fixed(1, 0, I):
                free_v = free_v[1:]
            if edge_fully_fixed(1, 1, I):
                free_v = free_v[:-1]

            if method == "tensor":
                # Per-quadrature-point row scaling of the WQ gather matrices --
                # a genuinely different (component-I-dependent) 1D pencil,
                # unlike "classic"/"mean" which only rescale eigenvalues of a
                # SHARED (component-independent) pencil after the fact.
                Ku_q, Kv_q = Kcoef[I]
                Mu_q, Mv_q = Mcoef[I]
                Mu = (wq_u.W00 @ wq_u.B0.multiply(Mu_q[:, None])).toarray()
                Ku = (wq_u.W10 @ wq_u.B1.multiply(Ku_q[:, None])).toarray()
                Mv = (wq_v.W00 @ wq_v.B0.multiply(Mv_q[:, None])).toarray()
                Kv = (wq_v.W10 @ wq_v.B1.multiply(Kv_q[:, None])).toarray()
            else:
                Mu, Ku, Mv, Kv = Mu0, Ku0, Mv0, Kv0

            eigval_u, Vu = sclin.eigh(Ku[np.ix_(free_u, free_u)], Mu[np.ix_(free_u, free_u)])
            eigval_v, Vv = sclin.eigh(Kv[np.ix_(free_v, free_v)], Mv[np.ix_(free_v, free_v)])

            tiny = 1e-12 * max(eigval_u.max(), eigval_v.max(), 1.0)
            eigval_u = np.where(eigval_u <= tiny, tiny, eigval_u)
            eigval_v = np.where(eigval_v <= tiny, tiny, eigval_v)

            if method == "mean":
                c0, c1 = mean_correctors[I, 0], mean_correctors[I, 1]
            else:
                c0, c1 = 1.0, 1.0
            # Kronecker SUM: Lambda[iv, iu] = c1*eigval_v[iv] + c0*eigval_u[iu]
            # ("tensor" method folds its own scaling into Ku/Kv/Mu/Mv already,
            # so c0=c1=1 there too), shape (len(free_v), len(free_u)) to
            # match the (nv, nu) reshape used in apply().
            Lambda = c1 * eigval_v[:, None] + c0 * eigval_u[None, :]

            self._free_u.append(free_u)
            self._free_v.append(free_v)
            self._Vu.append(Vu)
            self._Vv.append(Vv)
            self._lambda.append(Lambda)

            if diag_physical_full is None:
                self._diag_scale.append(None)
            else:
                Mu_f = Mu[np.ix_(free_u, free_u)]
                Ku_f = Ku[np.ix_(free_u, free_u)]
                Mv_f = Mv[np.ix_(free_v, free_v)]
                Kv_f = Kv[np.ix_(free_v, free_v)]
                diag_param_2d = (
                    c1 * np.outer(Kv_f.diagonal(), Mu_f.diagonal())
                    + c0 * np.outer(Mv_f.diagonal(), Ku_f.diagonal())
                )
                diag_phys_2d = diag_physical_full[I::n_dofs_per_cp].reshape(
                    self.nv, self.nu)[np.ix_(free_v, free_u)]
                self._diag_scale.append(np.sqrt(diag_param_2d / diag_phys_2d))

    def apply(self, r):
        """
        Approximate K^{-1} @ r. `r`/output: length `nctrl * n_dofs_per_cp`,
        interleaved per control point (component fastest) -- the same
        convention as `WQMatrixFreeStiffness.apply()` and
        `matrix_free_solver.solve()`.
        """
        r = np.asarray(r)
        out = np.zeros_like(r)
        for I in range(self.n_dofs_per_cp):
            r_I = r[I::self.n_dofs_per_cp].reshape(self.nv, self.nu)  # (iv, iu), C-order
            free_u, free_v = self._free_u[I], self._free_v[I]
            R = r_I[np.ix_(free_v, free_u)]

            scale = self._diag_scale[I]
            if scale is not None:
                R = scale * R

            Vu, Vv = self._Vu[I], self._Vv[I]
            Z = Vv.T @ R @ Vu
            Z /= self._lambda[I]
            R_out = Vv @ Z @ Vu.T

            if scale is not None:
                R_out = scale * R_out

            out_I = np.zeros((self.nv, self.nu))
            out_I[np.ix_(free_v, free_u)] = R_out
            out[I::self.n_dofs_per_cp] = out_I.ravel()
        return out
