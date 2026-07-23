"""
Tests for PatchEvaluator.evaluate_solution().

Properties verified:
  1. Missing DOF manager raises an exception at construction.
  2. Output shape is (n_pts, n_dofs_per_cp).
  3. Zero solution vector returns all-zero field.
  4. B-spline partition of unity: a constant field (u_x=c, u_y=0) is reproduced
     exactly at every evaluation point.
  5. NURBS partition of unity: same property holds for a rational patch.
  6. Linear field completeness (B-spline): u_x = x, u_y = 2y is reproduced
     exactly on a rectangular patch because linear fields are in the B-spline
     space for any degree >= 1.
  7. Uniform NURBS weights: w_a = const => R_a = N_a (weights cancel) =>
     PatchEvaluator must return the same values as on the B-spline patch.
  8. NURBS quarter ring: numerical agreement with the reference manual loop
     (the same computation expanded step by step in pure NumPy).
"""

import numpy as np
import pytest

from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager,
    Patch, GlobalDOFManager, PatchDOFManager,
    PatchEvaluator, SubdivisionRefiner,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _rect_patch(weights=None):
    """
    3x3-CP rectangular patch on [0,3]x[0,1], degree 2 in both directions.
    weights: list of 9 floats (NURBS) or None (pure B-spline).
    Returns (patch, pdm).
    """
    pts = [(x, y) for y in (0.0, 0.5, 1.0) for x in (0.0, 1.5, 3.0)]
    mgr = ControlPointManager(dim=2)
    for i, (x, y) in enumerate(pts):
        w = weights[i] if weights is not None else 1.0
        mgr.add_point([x, y], w)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    mapping = list(range(9))
    gm  = GlobalDOFManager([2] * 9)
    pdm = PatchDOFManager(2, mapping, gm)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 3], pdm)
    return patch, pdm


def _rect_patch_no_dof():
    """Same rectangle but WITHOUT a PatchDOFManager (geometry-only patch)."""
    mgr = ControlPointManager(dim=2)
    for y in (0.0, 0.5, 1.0):
        for x in (0.0, 1.5, 3.0):
            mgr.add_point([x, y])
    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    return Patch(BSplineSurface(su, sv), mgr, list(range(9)), [3, 3])


def _quarter_ring_patch():
    """
    NURBS quarter-ring: inner radius 1, outer radius 2, corner weights 1/sqrt(2).
    Subdivided 2 levels in each direction.
    Returns (patch, pdm) where pdm is the post-refinement dof manager.
    """
    r1, r2, wc = 1.0, 2.0, 1.0 / np.sqrt(2.0)
    mgr = ControlPointManager(dim=2)
    for (x, y, w) in [(r1, 0., 1.), (r1, r1, wc), (0., r1, 1.),
                      (r2, 0., 1.), (r2, r2, wc), (0., r2, 1.)]:
        mgr.add_point([x, y], w)
    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping = list(range(6))
    gm  = GlobalDOFManager([2] * 6)
    pdm = PatchDOFManager(2, mapping, gm)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 2], pdm)
    SubdivisionRefiner(direction=0, n_levels=2).refine(patch)
    SubdivisionRefiner(direction=1, n_levels=2).refine(patch)
    # After refinement, patch.dof_manager has been updated internally to cover
    # all the new CPs; the original pdm is stale.
    return patch, patch.dof_manager


def _eval_grid(patch, n=6):
    """Return a (n*n, 2) parameter array on a regular interior grid."""
    u_vals = np.linspace(0.05, 0.95, n)
    UU, VV = np.meshgrid(u_vals, u_vals)
    return np.column_stack([UU.ravel(), VV.ravel()])


def _manual_evaluate(patch, pdm, u_sol, params):
    """
    Reference implementation: manual NURBS interpolation loop in pure NumPy.
    Equivalent to the pedagogical loop in notebook 10, used as ground truth.
    """
    su, sv = patch.tensor.components
    p_u, p_v = su.degree, sv.degree
    n_u = patch.local_shape[0]
    gi  = patch.global_indices
    rational = patch.cp_manager.is_rational
    all_w = patch.cp_manager.weights_view() if rational else None

    result = np.zeros((len(params), 2))
    for k, (u_val, v_val) in enumerate(params):
        i_u = su.find_span(u_val)
        i_v = sv.find_span(v_val)
        N = np.outer(sv.basis_funs(i_v, v_val),
                     su.basis_funs(i_u, u_val)).ravel()
        af = [(i_u - p_u + du) + (i_v - p_v + dv) * n_u
              for dv in range(p_v + 1) for du in range(p_u + 1)]
        if rational:
            w_act = all_w[[gi[f] for f in af]]
            W = w_act @ N
            R = w_act * N / W
        else:
            R = N
        result[k, 0] = R @ [u_sol[pdm.get_global_dof_indices(f)[0]] for f in af]
        result[k, 1] = R @ [u_sol[pdm.get_global_dof_indices(f)[1]] for f in af]
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

def test_requires_dof_manager():
    """PatchEvaluator must raise when the patch has no DOF manager."""
    patch = _rect_patch_no_dof()
    with pytest.raises(Exception):
        PatchEvaluator(patch)


def test_output_shape():
    """evaluate_solution returns (n_pts, 2) for a 2-DOF-per-CP patch."""
    patch, pdm = _rect_patch()
    params = _eval_grid(patch, n=5)
    u_sol = np.zeros(2 * patch.n_cp)
    u_h = PatchEvaluator(patch).evaluate_solution(params, u_sol)
    assert u_h.shape == (len(params), 2)


def test_zero_solution():
    """A zero displacement vector must produce a zero field everywhere."""
    patch, pdm = _rect_patch()
    params = _eval_grid(patch, n=7)
    u_h = PatchEvaluator(patch).evaluate_solution(params, np.zeros(2 * patch.n_cp))
    np.testing.assert_allclose(u_h, 0.0, atol=1e-15)


def test_constant_field_bspline():
    """
    B-spline partition of unity: setting u_x = 1, u_y = 0 at every CP must
    reproduce (1, 0) at every evaluation point.
    """
    patch, pdm = _rect_patch()
    u_sol = np.tile([1.0, 0.0], patch.n_cp)  # u_x=1, u_y=0 for all CPs
    params = _eval_grid(patch, n=8)
    u_h = PatchEvaluator(patch).evaluate_solution(params, u_sol)
    np.testing.assert_allclose(u_h[:, 0], 1.0, atol=1e-14,
        err_msg="B-spline constant field u_x=1 not reproduced")
    np.testing.assert_allclose(u_h[:, 1], 0.0, atol=1e-14,
        err_msg="B-spline constant field u_y=0 not reproduced")


def test_constant_field_nurbs():
    """
    NURBS partition of unity: sum_a R_a = 1 (rational basis sums to 1), so a
    uniform displacement (u_x=3, u_y=5) must be reproduced exactly.
    """
    patch, pdm = _quarter_ring_patch()
    c_x, c_y = 3.0, 5.0
    u_sol = np.tile([c_x, c_y], patch.n_cp)
    params = _eval_grid(patch, n=7)
    u_h = PatchEvaluator(patch).evaluate_solution(params, u_sol)
    np.testing.assert_allclose(u_h[:, 0], c_x, atol=1e-13,
        err_msg="NURBS partition of unity violated for u_x")
    np.testing.assert_allclose(u_h[:, 1], c_y, atol=1e-13,
        err_msg="NURBS partition of unity violated for u_y")


def test_linear_field_bspline():
    """
    Linear completeness: u_x = x, u_y = 2y is in the B-spline space for any
    degree >= 1, so PatchEvaluator must reproduce it exactly.

    DOF values at each CP are set to (x_cp, 2*y_cp); at any evaluation point
    the interpolated field must match the physical coordinates (x, 2y).
    """
    patch, pdm = _rect_patch()
    n_cp = patch.n_cp

    # Set u_x_a = x_a, u_y_a = 2 * y_a at each control point
    u_sol = np.zeros(2 * n_cp)
    for local_idx in range(n_cp):
        cp = patch.control_point(local_idx)   # physical coords of this CP
        dofs = pdm.get_global_dof_indices(local_idx)
        u_sol[dofs[0]] = cp[0]        # u_x = x
        u_sol[dofs[1]] = 2.0 * cp[1]  # u_y = 2y

    params = _eval_grid(patch, n=8)
    su, sv = patch.tensor.components
    spans  = np.array([[su.find_span(u), sv.find_span(v)] for u, v in params],
                      dtype=np.int32)
    phys_pts = patch.evaluate_patch_nd_omp(spans, params)   # shape (n_pts, 2)

    u_h = PatchEvaluator(patch).evaluate_solution(params, u_sol)

    np.testing.assert_allclose(u_h[:, 0], phys_pts[:, 0], atol=1e-12,
        err_msg="Linear field u_x=x not reproduced by B-spline PatchEvaluator")
    np.testing.assert_allclose(u_h[:, 1], 2.0 * phys_pts[:, 1], atol=1e-12,
        err_msg="Linear field u_y=2y not reproduced by B-spline PatchEvaluator")


def test_uniform_nurbs_weights_equal_bspline():
    """
    When all NURBS weights are equal (w_a = 2.0), the rational basis satisfies
    R_a = w N_a / (w sum N_b) = N_a (weights cancel). PatchEvaluator must then
    return the same field values as the equivalent B-spline patch.
    """
    patch_bs, pdm_bs = _rect_patch(weights=None)
    patch_nr, pdm_nr = _rect_patch(weights=[2.0] * 9)

    assert not patch_bs.cp_manager.is_rational
    assert     patch_nr.cp_manager.is_rational

    rng = np.random.default_rng(42)
    u_sol = rng.standard_normal(2 * 9)

    params = _eval_grid(patch_bs, n=7)
    u_h_bs = PatchEvaluator(patch_bs).evaluate_solution(params, u_sol)
    u_h_nr = PatchEvaluator(patch_nr).evaluate_solution(params, u_sol)

    np.testing.assert_allclose(u_h_nr, u_h_bs, atol=1e-13,
        err_msg="Uniform NURBS weights should give identical result to B-spline")


def test_nurbs_quarter_ring_vs_manual_loop():
    """
    On the NURBS quarter ring, PatchEvaluator.evaluate_solution must agree with
    the reference manual NumPy loop to machine precision.
    """
    patch, pdm = _quarter_ring_patch()
    n_cp = patch.n_cp

    rng = np.random.default_rng(0)
    u_sol = rng.standard_normal(2 * n_cp)

    params = _eval_grid(patch, n=8)

    u_h_new = PatchEvaluator(patch).evaluate_solution(params, u_sol)
    u_h_ref = _manual_evaluate(patch, pdm, u_sol, params)

    np.testing.assert_allclose(u_h_new, u_h_ref, atol=1e-14,
        err_msg="PatchEvaluator disagrees with manual NURBS loop on quarter ring")


def test_nurbs_quarter_ring_refined_vs_manual_loop():
    """
    Same agreement check after additional k-refinement (larger CP count,
    more complex DOF layout).
    """
    from yeti_iga.future.bspline import PRefiner

    r1, r2, wc = 1.0, 2.0, 1.0 / np.sqrt(2.0)
    mgr = ControlPointManager(dim=2)
    for (x, y, w) in [(r1, 0., 1.), (r1, r1, wc), (0., r1, 1.),
                      (r2, 0., 1.), (r2, r2, wc), (0., r2, 1.)]:
        mgr.add_point([x, y], w)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping = list(range(6))
    gm  = GlobalDOFManager([2] * 6)
    pdm = PatchDOFManager(2, mapping, gm)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 2], pdm)

    # k-refinement: elevate to degree 3, then subdivide
    PRefiner(direction=0, n_elevations=1).refine(patch)
    PRefiner(direction=1, n_elevations=2).refine(patch)
    SubdivisionRefiner(direction=0, n_levels=2).refine(patch)
    SubdivisionRefiner(direction=1, n_levels=2).refine(patch)

    # Rebuild pdm for the refined patch
    n_cp = patch.n_cp
    gm2  = GlobalDOFManager([2] * n_cp)
    pdm2 = PatchDOFManager(2, list(range(n_cp)), gm2)
    from yeti_iga.future.bspline import BSplineSurface as _BSurf
    su_r, sv_r = patch.tensor.components
    patch_k = Patch(BSplineSurface(su_r, sv_r), patch.cp_manager,
                    list(patch.global_indices), list(patch.local_shape), pdm2)

    rng = np.random.default_rng(7)
    u_sol = rng.standard_normal(2 * n_cp)
    params = _eval_grid(patch_k, n=6)

    u_h_new = PatchEvaluator(patch_k).evaluate_solution(params, u_sol)
    u_h_ref = _manual_evaluate(patch_k, pdm2, u_sol, params)

    np.testing.assert_allclose(u_h_new, u_h_ref, atol=1e-14,
        err_msg="PatchEvaluator disagrees with manual loop after k-refinement")


if __name__ == '__main__':
    test_requires_dof_manager()
    test_output_shape()
    test_zero_solution()
    test_constant_field_bspline()
    test_constant_field_nurbs()
    test_linear_field_bspline()
    test_uniform_nurbs_weights_equal_bspline()
    test_nurbs_quarter_ring_vs_manual_loop()
    test_nurbs_quarter_ring_refined_vs_manual_loop()
    print("All PatchEvaluator tests passed.")
