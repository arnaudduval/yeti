"""
matrix_free_solver.solve(): iterative solve of K @ u = F with Dirichlet BCs,
K never assembled -- validated against a direct solve of the reference
Gauss-assembled (and sliced) system.
"""

import numpy as np
import pytest
import scipy.sparse.linalg as spla

from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager,
    Patch, GlobalDOFManager, PatchDOFManager, IGABasis1D,
    PatchIntegrator, Material, PlaneStress, WQMatrixFreeStiffness,
)
from yeti_iga.future.matrix_free_solver import solve as mf_solve


def _greville(degree, kv):
    n = len(kv) - degree - 1
    return np.array([np.mean(kv[i + 1:i + degree + 1]) for i in range(n)])


def _build_flat_patch(degree, nbel, Lx=4.0, Ly=1.0):
    """Exactly affine (Greville) flat rectangle -- see test_wq_stiffness.py."""
    kv_u = np.array([0.] * (degree + 1) + list(np.linspace(0, 1, nbel + 1)[1:-1]) + [1.] * (degree + 1))
    kv_v = np.array([0.] * (degree + 1) + list(np.linspace(0, 1, nbel + 1)[1:-1]) + [1.] * (degree + 1))
    su, sv = BSpline(degree, kv_u), BSpline(degree, kv_v)
    nu, nv = len(kv_u) - degree - 1, len(kv_v) - degree - 1
    xs, ys = _greville(degree, kv_u) * Lx, _greville(degree, kv_v) * Ly

    mgr = ControlPointManager(dim=2)
    mapping = []
    for jv in range(nv):
        for iu in range(nu):
            mgr.add_point([xs[iu], ys[jv]])
            mapping.append(jv * nu + iu)

    dof_manager = GlobalDOFManager([2] * mgr.n_points)
    pdm = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [nu, nv], pdm)
    return patch, pdm


def _fixed_dofs_on_edge(patch, pdm, direction=0, side=0):
    cps = patch.boundary_control_points(direction, side)
    return np.array(sorted(d for cp in cps for d in pdm.get_global_dof_indices(cp)))


@pytest.mark.parametrize("prescribed_nonzero", [False, True])
def test_solve_matches_direct_solve_of_reference_system(prescribed_nonzero):
    degree, nbel = 3, 6
    patch, pdm = _build_flat_patch(degree, nbel)
    basis_u = IGABasis1D.build(patch.tensor.components[0], degree + 1)
    basis_v = IGABasis1D.build(patch.tensor.components[1], degree + 1)
    law = PlaneStress(Material(E=210000., nu=0.3))

    K = PatchIntegrator(patch, basis_u, basis_v, law).integrate_stiffness()
    ndof = K.shape[0]

    fixed_dofs = _fixed_dofs_on_edge(patch, pdm, direction=0, side=0)
    free_dofs = np.setdiff1d(np.arange(ndof), fixed_dofs)

    rng = np.random.default_rng(0)
    F = rng.standard_normal(ndof)
    prescribed = rng.standard_normal(len(fixed_dofs)) * 0.01 if prescribed_nonzero else None
    prescribed_vals = np.zeros(len(fixed_dofs)) if prescribed is None else prescribed

    # Reference: direct solve of K_ff @ u_f = F_f - K_fc @ u_c on the
    # explicitly assembled (and sliced) Gauss stiffness matrix.
    K_csr = K.tocsr()
    K_ff = K_csr[free_dofs][:, free_dofs]
    K_fc = K_csr[free_dofs][:, fixed_dofs]
    u_ref = np.zeros(ndof)
    u_ref[fixed_dofs] = prescribed_vals
    u_ref[free_dofs] = spla.spsolve(K_ff.tocsc(), F[free_dofs] - K_fc @ prescribed_vals)

    wq_stiff = WQMatrixFreeStiffness(patch, law)
    u_mf, info = mf_solve(wq_stiff.apply, ndof, F, fixed_dofs,
                          prescribed=prescribed, rtol=1e-10, maxiter=2000)

    assert info == 0
    assert np.allclose(u_mf[fixed_dofs], prescribed_vals)
    # GMRES (the default -- see matrix_free_solver's module docstring for
    # why not CG) settles a bit short of CG's machine-precision match on
    # this well-conditioned affine case; still several orders of magnitude
    # tighter than anything this test cares about distinguishing.
    assert np.linalg.norm(u_mf - u_ref) / np.linalg.norm(u_ref) < 1e-6
