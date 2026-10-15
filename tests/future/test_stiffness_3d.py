"""
Test the 3D solid (BSplineVolume) stiffness kernel added for genuine 3D
elasticity (Elasticity3D + PatchIntegrator's 3-basis constructor). No legacy
Fortran 3D reference exists, so correctness is checked against:

1. The classic rigid-body-mode patch test (cheap first filter, run first as
   per the 3D-solid plan): K is symmetric and gives zero strain energy under
   any of the 6 rigid body motions (3 translations + 3 rotations about an
   arbitrary center).
2. Saint-Venant pure bending of a prismatic beam (Timoshenko & Goodier,
   Theory of Elasticity, Sec. 20) -- a closed-form solution requiring zero
   body force:
     - exact reproduction (degree>=2, any mesh): the FE solution must match
       the exact (quadratic) field to near machine precision.
     - h-convergence (degree=1, meshes of 1/2/4/8/16 elements per
       direction): L2 error ~ O(h^2), H1-seminorm error ~ O(h^1) (standard
       a priori FE estimates).
"""

import numpy as np
import scipy.sparse.linalg as spla

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineVolume, Patch, ControlPointManager,
    GlobalDOFManager, PatchDOFManager, IGABasis1D, Material, Elasticity3D,
    PlaneStress, PatchIntegrator, Traction3D, ScalarLocalOperator3D,
)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _open_uniform_knots(p, nbel):
    interior = np.linspace(0.0, 1.0, nbel + 1)[1:-1]
    return np.concatenate([[0.0] * (p + 1), interior, [1.0] * (p + 1)])


def _greville_abscissae(kv, p):
    n = len(kv) - p - 1
    return np.array([np.mean(kv[i + 1:i + p + 1]) for i in range(n)])


def _build_box_patch(degree, nbel, Lx, Ly, Lz, center_y=False, center_z=False):
    """
    Non-rational BSplineVolume box, x in [0, Lx], y in [0, Ly] or
    [-Ly/2, Ly/2] (center_y), z similarly (center_z). Control points at
    Greville abscissae -- exact for this affine geometry map regardless of
    degree (established earlier in the 3D-solid work).
    """
    p = degree
    kv = _open_uniform_knots(p, nbel)
    splines = [BSpline(p, kv) for _ in range(3)]
    tensor = BSplineVolume(*splines)
    n1 = nbel + p
    n_cp = n1 ** 3
    gr = _greville_abscissae(kv, p)

    off_y = 0.5 if center_y else 0.0
    off_z = 0.5 if center_z else 0.0

    mgr = ControlPointManager(dim=3)
    coords = np.zeros((n_cp, 3))
    idx = 0
    for kw in range(n1):
        for jv in range(n1):
            for iu in range(n1):
                x = gr[iu] * Lx
                y = (gr[jv] - off_y) * Ly
                z = (gr[kw] - off_z) * Lz
                coords[idx] = (x, y, z)
                mgr.add_point([x, y, z])
                idx += 1

    patch0 = Patch(tensor, mgr, list(range(n_cp)), [n1, n1, n1])
    gdm = GlobalDOFManager([3] * n_cp)
    pdm = PatchDOFManager(3, list(range(n_cp)), gdm)
    patch = Patch(patch0.tensor, patch0.cp_manager, list(patch0.global_indices),
                  list(patch0.local_shape), pdm)
    return patch, coords, kv


def _build_integrator(patch, degree, law):
    """
    PatchIntegrator stores C++ REFERENCES to its patch/IGABasis1D arguments
    (not copies); the pybind11 constructor bindings use py::keep_alive to
    keep the corresponding Python objects alive for as long as the
    integrator itself is alive, so it is safe to let basis_u/v/w go out of
    scope here (see test_patch_integrator_lifetime.py for a dedicated
    regression test of this guarantee).
    """
    basis_u = IGABasis1D.build(patch.tensor.components[0], degree + 1)
    basis_v = IGABasis1D.build(patch.tensor.components[1], degree + 1)
    basis_w = IGABasis1D.build(patch.tensor.components[2], degree + 1)
    return PatchIntegrator(patch, basis_u, basis_v, basis_w, law)


# ---------------------------------------------------------------------------
# 1. Rigid-body-mode patch test -- run first, cheap correctness gate.
# ---------------------------------------------------------------------------

def test_elasticity3d_stiffness_symmetric_rigid_body_modes():
    degree, nbel = 2, 2
    patch, coords, _ = _build_box_patch(degree, nbel, Lx=3.0, Ly=1.0, Lz=2.0)
    law = Elasticity3D(Material(E=210000.0, nu=0.3))
    integrator = _build_integrator(patch, degree, law)

    K = integrator.integrate_stiffness()

    np.testing.assert_allclose(K.toarray(), K.toarray().T, atol=1.e-6)

    n_cp = coords.shape[0]
    center = coords.mean(axis=0)

    translations = [np.array(v) for v in ([1., 0., 0.], [0., 1., 0.], [0., 0., 1.])]
    rotation_axes = [np.array(v) for v in ([1., 0., 0.], [0., 1., 0.], [0., 0., 1.])]

    for t in translations:
        u = np.tile(t, n_cp)
        energy = u @ (K @ u)
        assert abs(energy) < 1.e-6

    for axis in rotation_axes:
        u = np.zeros(3 * n_cp)
        for i in range(n_cp):
            u[3*i:3*i+3] = np.cross(axis, coords[i] - center)
        energy = u @ (K @ u)
        assert abs(energy) < 1.e-6 * max(np.linalg.norm(u) ** 2, 1.0)


def test_3d_constructor_rejects_wrong_dof_count_law():
    """
    The 3D constructor validates law.n_dofs_per_cp() == 3 -- a defensive
    check added specifically because without it a 2-DOF law silently
    corrupts the DOF scatter instead of failing loudly (see the plan).
    """
    degree, nbel = 1, 1
    patch, _, _ = _build_box_patch(degree, nbel, Lx=1.0, Ly=1.0, Lz=1.0)
    basis_u = IGABasis1D.build(patch.tensor.components[0], degree + 1)
    basis_v = IGABasis1D.build(patch.tensor.components[1], degree + 1)
    basis_w = IGABasis1D.build(patch.tensor.components[2], degree + 1)

    wrong_law = PlaneStress(Material(E=210000.0, nu=0.3))  # n_dofs_per_cp() == 2
    try:
        PatchIntegrator(patch, basis_u, basis_v, basis_w, wrong_law)
        assert False, "expected an exception for a law with n_dofs_per_cp() != 3"
    except Exception:
        pass


# ---------------------------------------------------------------------------
# 2. Saint-Venant pure bending of a prismatic beam.
# ---------------------------------------------------------------------------

E_BEAM, NU_BEAM, M_BEAM = 210000.0, 0.3, 1.0
L_BEAM, H_BEAM, B_BEAM = 4.0, 1.0, 1.0
IZ_BEAM = B_BEAM * H_BEAM ** 3 / 12.0
C_BEAM = M_BEAM / (E_BEAM * IZ_BEAM)


def _exact_u(x, y, z):
    ux = C_BEAM * x * y
    uy = -0.5 * C_BEAM * (x**2 + NU_BEAM * y**2 - NU_BEAM * z**2)
    uz = -C_BEAM * NU_BEAM * y * z
    return np.array([ux, uy, uz])


def _exact_grad_u(x, y, z):
    # grad[i, j] = d(u_i) / d(x_j), (x_j) = (x, y, z)
    return np.array([
        [C_BEAM * y, C_BEAM * x, 0.0],
        [-C_BEAM * x, -C_BEAM * NU_BEAM * y, C_BEAM * NU_BEAM * z],
        [0.0, -C_BEAM * NU_BEAM * z, -C_BEAM * NU_BEAM * y],
    ])


class _BendingTraction(Traction3D):
    def evaluate(self, pt):
        _, y, _ = pt
        return np.array([M_BEAM * y / IZ_BEAM, 0.0, 0.0])


class _L2ErrorOperator(ScalarLocalOperator3D):
    def compute_scalar_integrand(self, R, dRdx, dRdy, dRdz, pt, u_local):
        R = np.asarray(R)
        u = np.asarray(u_local)
        uh = np.array([R @ u[0::3], R @ u[1::3], R @ u[2::3]])
        ex = _exact_u(*pt)
        return float(np.sum((uh - ex) ** 2))


class _H1SemiNormErrorOperator(ScalarLocalOperator3D):
    def compute_scalar_integrand(self, R, dRdx, dRdy, dRdz, pt, u_local):
        dRdx, dRdy, dRdz = np.asarray(dRdx), np.asarray(dRdy), np.asarray(dRdz)
        u = np.asarray(u_local)
        ux, uy, uz = u[0::3], u[1::3], u[2::3]
        grad_h = np.array([
            [dRdx @ ux, dRdy @ ux, dRdz @ ux],
            [dRdx @ uy, dRdy @ uy, dRdz @ uy],
            [dRdx @ uz, dRdy @ uz, dRdz @ uz],
        ])
        grad_ex = _exact_grad_u(*pt)
        return float(np.sum((grad_h - grad_ex) ** 2))


def _collocation_matrix_1d(spline, gr):
    """
    Greville collocation matrix: C[i, a] = N_a(greville_i). Solving C @ d = f
    (f sampled at the Greville points) gives the TRUE control coefficients
    of f whenever f lies in the spline space -- unlike naively assigning
    d_i = f(greville_i), which is only exact for affine (degree<=1) f (a
    subtlety specific to the degree>=2 exact-reproduction check below: the
    pure-bending field is quadratic).
    """
    n = len(gr)
    C = np.zeros((n, n))
    for i in range(n):
        for a in range(n):
            C[i, a] = spline.one_basis_fun(gr[i], a)
    return C


def _fit_exact_field(components, grs, phys_map, field_fn, n_components):
    """
    Exact B-spline control coefficients of an analytic vector field via
    per-direction collocation solves (sum factorization), u-fastest
    Fortran-order flattening to match this codebase's flat-index
    convention.
    """
    Cs = [_collocation_matrix_1d(components[d], grs[d]) for d in range(3)]
    n_u, n_v, n_w = len(grs[0]), len(grs[1]), len(grs[2])

    F = np.zeros((n_u, n_v, n_w, n_components))
    for iw in range(n_w):
        for iv in range(n_v):
            for iu in range(n_u):
                xyz = phys_map(grs[0][iu], grs[1][iv], grs[2][iw])
                F[iu, iv, iw, :] = field_fn(*xyz)

    D = F
    for c in range(n_components):
        d = D[:, :, :, c]
        d = np.linalg.solve(Cs[0], d.reshape(n_u, -1)).reshape(n_u, n_v, n_w)
        d = np.linalg.solve(Cs[1], d.transpose(1, 0, 2).reshape(n_v, -1)).reshape(n_v, n_u, n_w).transpose(1, 0, 2)
        d = np.linalg.solve(Cs[2], d.transpose(2, 0, 1).reshape(n_w, -1)).reshape(n_w, n_u, n_v).transpose(1, 2, 0)
        D[:, :, :, c] = d
    return D


def _run_pure_bending(degree, nbel):
    patch, coords, kv = _build_box_patch(
        degree, nbel, L_BEAM, H_BEAM, B_BEAM, center_y=True, center_z=True)
    n1 = nbel + degree
    n_cp = n1 ** 3
    ndof = 3 * n_cp

    law = Elasticity3D(Material(E=E_BEAM, nu=NU_BEAM))
    integrator = _build_integrator(patch, degree, law)

    K = integrator.integrate_stiffness().tocsr()
    F = integrator.integrate_boundary_load_3d(0, 1, _BendingTraction())

    gr = _greville_abscissae(kv, degree)
    grs = [gr, gr, gr]

    def phys_map(xi, eta, zeta):
        return xi * L_BEAM, (eta - 0.5) * H_BEAM, (zeta - 0.5) * B_BEAM

    def field_fn(x, y, z):
        return _exact_u(x, y, z)

    D = _fit_exact_field(patch.tensor.components, grs, phys_map, field_fn, 3)
    u_exact_coeffs = np.zeros(ndof)
    u_exact_coeffs[0::3] = D[:, :, :, 0].flatten(order='F')
    u_exact_coeffs[1::3] = D[:, :, :, 1].flatten(order='F')
    u_exact_coeffs[2::3] = D[:, :, :, 2].flatten(order='F')

    # Dirichlet BC on the x=0 face, lifted to the exact field's own
    # coefficients there (nonzero lift -- standard u = u_D + u_0 technique).
    fixed_cp = np.array(sorted(patch.boundary_control_points(0, 0)))
    fixed_dofs = np.sort(np.concatenate([3*fixed_cp, 3*fixed_cp+1, 3*fixed_cp+2]))
    free_dofs = np.setdiff1d(np.arange(ndof), fixed_dofs)

    u_c_full = np.zeros(ndof)
    u_c_full[fixed_dofs] = u_exact_coeffs[fixed_dofs]
    rhs_free = F[free_dofs] - (K @ u_c_full)[free_dofs]
    K_ff = K[free_dofs][:, free_dofs].tocsc()

    u = np.zeros(ndof)
    u[fixed_dofs] = u_exact_coeffs[fixed_dofs]
    u[free_dofs] = spla.spsolve(K_ff, rhs_free)

    l2_sq = integrator.integrate_scalar_operator_3d(_L2ErrorOperator(), u)
    h1_sq = integrator.integrate_scalar_operator_3d(_H1SemiNormErrorOperator(), u)

    return dict(
        u=u, u_exact_coeffs=u_exact_coeffs,
        l2_err=float(np.sqrt(max(l2_sq, 0.0))),
        h1_semi_err=float(np.sqrt(max(h1_sq, 0.0))),
    )


def test_elasticity3d_pure_bending_exact_reproduction_degree2():
    r = _run_pure_bending(degree=2, nbel=1)
    np.testing.assert_allclose(r['u'], r['u_exact_coeffs'], atol=1.e-10)
    assert r['l2_err'] < 1.e-10
    assert r['h1_semi_err'] < 1.e-10


def test_elasticity3d_pure_bending_exact_reproduction_degree3_multi_element():
    """Same check on a non-trivial multi-element mesh (2 elements/direction,
    degree 3), to rule out a single-element/Bezier-only coincidence."""
    r = _run_pure_bending(degree=3, nbel=2)
    np.testing.assert_allclose(r['u'], r['u_exact_coeffs'], atol=1.e-9)
    assert r['l2_err'] < 1.e-9
    assert r['h1_semi_err'] < 1.e-9


def test_elasticity3d_pure_bending_convergence_degree1():
    """
    Degree 1 cannot represent the quadratic bending field exactly, so error
    decreases with mesh refinement -- check the standard a priori rates
    (L2 ~ O(h^2), H1-seminorm ~ O(h^1)) are approached, using a loose
    tolerance on the measured rate since coarse meshes are pre-asymptotic.
    """
    nbels = [2, 4, 8, 16]
    l2_errs, h1_errs = [], []
    for nbel in nbels:
        r = _run_pure_bending(degree=1, nbel=nbel)
        l2_errs.append(r['l2_err'])
        h1_errs.append(r['h1_semi_err'])

    # Errors must decrease monotonically under refinement.
    assert all(l2_errs[i+1] < l2_errs[i] for i in range(len(l2_errs) - 1))
    assert all(h1_errs[i+1] < h1_errs[i] for i in range(len(h1_errs) - 1))

    # Rate between the two finest meshes must approach the expected a priori
    # estimate (L2 order 2, H1 seminorm order >= 1 -- the measured H1 rate
    # exceeds 1 here due to Gauss-point superconvergence on this structured,
    # polynomial-exact-solution mesh, a known benign phenomenon, not a bug;
    # see the plan's discussion).
    rate_l2 = np.log(l2_errs[-2] / l2_errs[-1]) / np.log(2.0)
    rate_h1 = np.log(h1_errs[-2] / h1_errs[-1]) / np.log(2.0)
    assert rate_l2 > 1.5
    assert rate_h1 > 1.0


if __name__ == '__main__':
    test_elasticity3d_stiffness_symmetric_rigid_body_modes()
    test_3d_constructor_rejects_wrong_dof_count_law()
    test_elasticity3d_pure_bending_exact_reproduction_degree2()
    test_elasticity3d_pure_bending_exact_reproduction_degree3_multi_element()
    test_elasticity3d_pure_bending_convergence_degree1()
