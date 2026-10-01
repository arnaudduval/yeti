"""
WQMatrixFreeStiffness: matrix-free elasticity stiffness action K @ v (Phase
C), validated against PatchIntegrator.integrate_stiffness() @ v -- the
already-validated Gauss-based assembly -- on a flat patch and a curved
(NURBS) patch, per the plan's acceptance criteria.
"""

import numpy as np
import pytest

from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager,
    Patch, GlobalDOFManager, PatchDOFManager, IGABasis1D,
    PatchIntegrator, Material, PlaneStress, PlaneStrain,
    WQMatrixFreeStiffness, SubdivisionRefiner, PRefiner,
)


def _greville(degree, kv):
    n = len(kv) - degree - 1
    return np.array([np.mean(kv[i + 1:i + degree + 1]) for i in range(n)])


def _build_flat_patch(degree, nbel_u, nbel_v, Lx=2.0, Ly=1.0):
    """
    An EXACTLY affine (constant-Jacobian) flat rectangle: control points at
    the Greville abscissae scaled by the physical dimensions -- B-spline
    bases reproduce a linear map exactly at Greville points, so this is not
    just "flat-looking" but truly affine (zero curvature), unlike a naive
    equally-spaced control net (which is NOT exactly affine for degree >= 2
    and introduces spurious geometric error unrelated to WQ itself).
    """
    kv_u = np.array([0.] * (degree + 1) + list(np.linspace(0, 1, nbel_u + 1)[1:-1]) + [1.] * (degree + 1))
    kv_v = np.array([0.] * (degree + 1) + list(np.linspace(0, 1, nbel_v + 1)[1:-1]) + [1.] * (degree + 1))
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
    return patch


def _build_quarter_ring(refine_u=0, refine_v=0, elevate_v=0):
    """NURBS quarter ring (inner radius 1, outer radius 2), same fixture as
    test_nurbs.py's _build_quarter_ring, with an added optional degree
    elevation in v (radial direction starts at degree 1)."""
    r1, r2 = 1.0, 2.0
    w_c = 1.0 / np.sqrt(2.0)
    mgr = ControlPointManager(dim=2)
    for (x, y, w) in [
        (r1, 0., 1.0), (r1, r1, w_c), (0., r1, 1.0),
        (r2, 0., 1.0), (r2, r2, w_c), (0., r2, 1.0),
    ]:
        mgr.add_point([x, y], w)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping = list(range(6))
    dof_mgr = GlobalDOFManager([2] * 6)
    pdm = PatchDOFManager(2, mapping, dof_mgr)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 2], pdm)

    if elevate_v:
        PRefiner(direction=1, n_elevations=elevate_v).refine(patch)
    for _ in range(refine_u):
        SubdivisionRefiner(direction=0, n_levels=1).refine(patch)
    for _ in range(refine_v):
        SubdivisionRefiner(direction=1, n_levels=1).refine(patch)

    return patch


def _max_rel_err(patch, law, n_trials=3, seed=0):
    basis_u = IGABasis1D.build(patch.tensor.components[0], patch.tensor.components[0].degree + 1)
    basis_v = IGABasis1D.build(patch.tensor.components[1], patch.tensor.components[1].degree + 1)
    K = PatchIntegrator(patch, basis_u, basis_v, law).integrate_stiffness()
    wq_stiff = WQMatrixFreeStiffness(patch, law)

    rng = np.random.default_rng(seed)
    ndof = K.shape[0]
    errs = []
    for _ in range(n_trials):
        v = rng.standard_normal(ndof)
        Kv_gauss = K @ v
        Kv_mf = wq_stiff.apply(v)
        errs.append(np.linalg.norm(Kv_mf - Kv_gauss) / np.linalg.norm(Kv_gauss))
    return max(errs)


@pytest.mark.parametrize("degree,nbel", [(2, 4), (3, 3)])
@pytest.mark.parametrize("law_cls", [PlaneStress, PlaneStrain])
def test_flat_patch_matches_gauss_to_machine_precision(degree, nbel, law_cls):
    """
    An exactly affine geometry keeps every integrand inside the polynomial
    space WQ is constructed to reproduce exactly -- so, unlike the curved
    case below, K @ v should match the Gauss assembly to machine precision,
    not just approximately.
    """
    patch = _build_flat_patch(degree, nbel, nbel)
    law = law_cls(Material(E=210000., nu=0.3))
    assert _max_rel_err(patch, law) < 1e-10


def test_curved_nurbs_patch_converges_under_refinement():
    """
    NURBS quarter-ring: the physical integrand is rational (not polynomial)
    after the Jacobian pullback, so WQ is only approximate here -- but the
    approximation must be a genuine, convergent one (error shrinking under
    mesh refinement), not "close by accident" at one mesh size.
    """
    law = PlaneStress(Material(E=210000., nu=0.3))
    errors = [
        _max_rel_err(_build_quarter_ring(refine_u=r, refine_v=r, elevate_v=1), law)
        for r in range(4)
    ]
    # Strictly decreasing, and the finest mesh is well under the coarsest.
    assert all(errors[i + 1] < errors[i] for i in range(len(errors) - 1))
    assert errors[-1] < 0.05 * errors[0]


def test_rational_patch_applies_quotient_rule_to_the_displacement_field():
    """
    Regression guard for a real bug found by comparing against pymfiga_jcf
    (an external, independent WQ implementation): on a NURBS patch (weights
    != 1, e.g. this quarter-ring's w=1/sqrt(2) control points), apply() must
    apply the SAME quotient-rule correction to the displacement field that
    the constructor already applies to the geometry map -- R_a = w_a*N_a/W(xi)
    is rational, so both the trial field and the implicit test function need
    it. Before that correction existed, this exact fixture (1-element
    quarter-ring) was off from Gauss by ~23% (see WQMatrixFreeStiffness.hpp's
    "NURBS RATIONALITY" note) -- comfortably tighter than the ordinary WQ
    approximation error a rational patch is expected to have.
    """
    law = PlaneStress(Material(E=210000., nu=0.3))
    patch = _build_quarter_ring(refine_u=0, refine_v=0, elevate_v=1)
    assert _max_rel_err(patch, law) < 0.10
