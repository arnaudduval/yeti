"""
Tests for NURBS (Non-Uniform Rational B-Spline) support in PatchIntegrator.

Key properties verified:
  1. B-spline fast path unchanged -- ControlPointManager.is_rational stays False
     and the entire existing test suite (test_stiffness.py, test_mass.py etc.)
     remains unaffected.
  2. Uniform weights (all w == 1.0) -> NURBS must reproduce B-spline exactly.
  3. Quarter-circle geometry (canonical NURBS reference): degree-2, 3 CPs with
     weights (1, 1/sqrt(2), 1), physical points at (1,0), (1/sqrt(2),1/sqrt(2)),
     (0,1) for u in [0, 0.5, 1].
  4. Constant non-unit weights: forces rational mode but all weights cancel in
     the formula -> same result as B-spline (sanity check of the rationalization).
"""

import os
import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineSurface, ControlPointManager,
    Patch, GlobalDOFManager, PatchDOFManager, IGABasis1D,
    PatchIntegrator, Material, PlaneStress,
    HRefiner, PRefiner, SubdivisionRefiner,
)


def test_bspline_is_not_rational():
    """Creating a manager without weights must keep is_rational False."""
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.0, 0.0])
    assert not mgr.is_rational


def test_add_point_weighted_activates_rational():
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])           # w=1.0 -- B-spline so far
    mgr.add_point([1.0, 0.0], 0.7071)  # w!=1.0 -- activates rational mode
    mgr.add_point([1.0, 1.0])          # w defaults to 1.0, but mode stays on
    assert mgr.is_rational
    assert len(mgr.weights_view()) == 3
    np.testing.assert_allclose(mgr.weights_view(), [1.0, 0.7071, 1.0])


def test_set_weight_activates_and_reverts():
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1.0, 0.0])
    assert not mgr.is_rational

    mgr.set_weight(0, 2.0)
    assert mgr.is_rational

    mgr.set_weight(0, 1.0)  # back to 1.0 -- should revert to B-spline fast path
    assert not mgr.is_rational


def _build_rect_patch():
    """Simple 3x1 rectangle, degree 2, used for uniform-weight comparisons."""
    mgr = ControlPointManager(dim=2)
    for (x, y) in [(0.0, 0.0), (1.5, 0.0), (3.0, 0.0),
                   (0.0, 0.5), (1.5, 0.5), (3.0, 0.5),
                   (0.0, 1.0), (1.5, 1.0), (3.0, 1.0)]:
        mgr.add_point([x, y])

    dof_manager = GlobalDOFManager([2] * mgr.n_points)
    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    mapping = list(range(9))
    pdm = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 3], pdm)
    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)
    return patch, basis_u, basis_v, mgr


def _build_rect_patch_nurbs(w_uniform):
    """Same rectangle but with all weights set to w_uniform != 1."""
    mgr = ControlPointManager(dim=2)
    for (x, y) in [(0.0, 0.0), (1.5, 0.0), (3.0, 0.0),
                   (0.0, 0.5), (1.5, 0.5), (3.0, 0.5),
                   (0.0, 1.0), (1.5, 1.0), (3.0, 1.0)]:
        mgr.add_point([x, y], w_uniform)

    dof_manager = GlobalDOFManager([2] * mgr.n_points)
    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    mapping = list(range(9))
    pdm = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 3], pdm)
    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)
    return patch, basis_u, basis_v, mgr


def test_uniform_weights_reproduce_bspline_stiffness():
    """
    NURBS with all weights = 2.0 must give the same stiffness matrix as the
    equivalent B-spline: the rational basis R_a = 2*N_a/(sum 2*N_b) = N_a
    when all weights are equal (they cancel in the quotient).
    """
    patch_bs, bu_bs, bv_bs, _ = _build_rect_patch()
    patch_nr, bu_nr, bv_nr, _ = _build_rect_patch_nurbs(2.0)

    mat = PlaneStress(Material(E=210000, nu=0.3))

    K_bs = PatchIntegrator(patch_bs, bu_bs, bv_bs, mat).integrate_stiffness()
    K_nr = PatchIntegrator(patch_nr, bu_nr, bv_nr, mat).integrate_stiffness()

    np.testing.assert_allclose(K_nr.toarray(), K_bs.toarray(), rtol=1e-12, atol=1e-10)


def test_uniform_weights_reproduce_bspline_mass():
    """Same as above but for the mass matrix."""
    patch_bs, bu_bs, bv_bs, _ = _build_rect_patch()
    patch_nr, bu_nr, bv_nr, _ = _build_rect_patch_nurbs(3.0)

    mat = PlaneStress(Material(E=210000, nu=0.3, rho=7800.0))

    M_bs = PatchIntegrator(patch_bs, bu_bs, bv_bs, mat).integrate_mass()
    M_nr = PatchIntegrator(patch_nr, bu_nr, bv_nr, mat).integrate_mass()

    np.testing.assert_allclose(M_nr.toarray(), M_bs.toarray(), rtol=1e-12, atol=1e-10)


def test_quarter_circle_geometry():
    """
    Canonical NURBS quarter-circle (degree 2, 3 CPs, open knot vector).

    Control points in physical space: P0=(1,0), P1=(1,1), P2=(0,1).
    Weights: w0=1, w1=1/sqrt(2), w2=1.

    Exact evaluation check: at u=0 -> (1,0), u=0.5 -> (1/sqrt(2), 1/sqrt(2)),
    u=1 -> (0,1) on the unit circle.
    """
    kv = np.array([0., 0., 0., 1., 1., 1.])
    su = BSpline(2, kv)

    w1 = 1.0 / np.sqrt(2.0)

    mgr = ControlPointManager(dim=2)
    mgr.add_point([1.0, 0.0], 1.0)
    mgr.add_point([1.0, 1.0], w1)
    mgr.add_point([0.0, 1.0], 1.0)
    assert mgr.is_rational

    dof_mgr = GlobalDOFManager([2] * 3)
    mapping = [0, 1, 2]
    local_shape = [3]
    # Build a 1D "patch" as a 2D patch with trivial v direction (single point).
    sv = BSpline(0, np.array([0., 1.]))
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 1],
                  PatchDOFManager(2, mapping, dof_mgr))

    # Evaluate at u=0, 0.5, 1 (v=0 throughout, single span)
    for u, expected in [(0.0,   [1.0,          0.0         ]),
                        (0.5,   [1/np.sqrt(2), 1/np.sqrt(2)]),
                        (1.0,   [0.0,          1.0         ])]:
        span_u = su.find_span(u)
        span_v = sv.find_span(0.0)
        pt = patch.evaluate_patch_nd(
            np.array([[span_u, span_v]], dtype=np.int32),
            np.array([[u, 0.0]]))[0]
        np.testing.assert_allclose(pt, expected, atol=1e-12,
            err_msg=f"Quarter-circle evaluation failed at u={u}")
        # Also verify the point lies on the unit circle
        np.testing.assert_allclose(np.linalg.norm(pt), 1.0, atol=1e-12,
            err_msg=f"Point at u={u} is not on the unit circle")


def _build_quarter_circle_patch():
    """
    Canonical NURBS quarter-circle (degree 2, 3 CPs).
    P0=(1,0) w=1, P1=(1,1) w=1/sqrt(2), P2=(0,1) w=1.
    Returns (patch, mgr) ready to be refined.
    """
    kv = np.array([0., 0., 0., 1., 1., 1.])
    su = BSpline(2, kv)
    sv = BSpline(0, np.array([0., 1.]))

    w1 = 1.0 / np.sqrt(2.0)
    mgr = ControlPointManager(dim=2)
    mgr.add_point([1.0, 0.0], 1.0)
    mgr.add_point([1.0, 1.0], w1)
    mgr.add_point([0.0, 1.0], 1.0)

    dof_mgr = GlobalDOFManager([2] * 3)
    mapping = [0, 1, 2]
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 1],
                  PatchDOFManager(2, mapping, dof_mgr))
    return patch, mgr


def _eval_circle_patch(patch):
    """Evaluate the quarter-circle patch at u in {0, 0.5, 1}."""
    splines = patch.tensor.components
    results = []
    for u in [0.0, 0.5, 1.0]:
        span_u = splines[0].find_span(u)
        span_v = splines[1].find_span(0.0)
        pt = patch.evaluate_patch_nd(
            np.array([[span_u, span_v]], dtype=np.int32),
            np.array([[u, 0.0]]))[0]
        results.append(pt)
    return results


def test_nurbs_knot_insertion_preserves_quarter_circle():
    """
    After inserting knot 0.5 into the quarter-circle NURBS patch, the
    geometry must be unchanged: evaluated points still lie on the unit circle
    and the physical coordinates at u=0, 0.5, 1 are unchanged.
    """
    patch, mgr = _build_quarter_circle_patch()
    pts_before = _eval_circle_patch(patch)

    HRefiner(0, 0.5).refine(patch)

    pts_after = _eval_circle_patch(patch)

    # Geometry unchanged
    for pb, pa in zip(pts_before, pts_after):
        np.testing.assert_allclose(pa, pb, atol=1e-12,
            err_msg="Knot insertion moved a point on the quarter-circle")

    # All points lie on the unit circle
    for pt in pts_after:
        np.testing.assert_allclose(np.linalg.norm(pt), 1.0, atol=1e-12,
            err_msg=f"Point {pt} is not on the unit circle after knot insertion")

    # Weights are correctly stored (4 CPs after one knot insertion)
    assert mgr.is_rational
    assert len(mgr.weights_view()) == len(patch.global_indices)


def test_nurbs_degree_elevation_preserves_quarter_circle():
    """
    After elevating the quarter-circle NURBS patch to degree 3, the
    geometry must be unchanged.
    """
    patch, mgr = _build_quarter_circle_patch()
    pts_before = _eval_circle_patch(patch)

    PRefiner(0, 1).refine(patch)

    pts_after = _eval_circle_patch(patch)

    for pb, pa in zip(pts_before, pts_after):
        np.testing.assert_allclose(pa, pb, atol=1e-12,
            err_msg="Degree elevation moved a point on the quarter-circle")

    for pt in pts_after:
        np.testing.assert_allclose(np.linalg.norm(pt), 1.0, atol=1e-12,
            err_msg=f"Point {pt} is not on the unit circle after degree elevation")

    # Weights are correctly stored (4 CPs after elevation from degree 2 to 3)
    assert mgr.is_rational
    assert len(mgr.weights_view()) == len(patch.global_indices)


def _build_quarter_ring(refine_u=0, refine_v=0):
    """
    NURBS quarter ring: inner radius 1, outer radius 2.
    pu=2 (angular), pv=1 (radial).  Optional subdivision refinement.
    """
    r1, r2 = 1.0, 2.0
    w_c = 1.0 / np.sqrt(2.0)

    mgr = ControlPointManager(dim=2)
    for (x, y, w) in [
        (r1, 0.,  1.0), (r1, r1,  w_c), (0., r1,  1.0),
        (r2, 0.,  1.0), (r2, r2,  w_c), (0., r2,  1.0),
    ]:
        mgr.add_point([x, y], w)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    mapping = list(range(6))
    dof_mgr = GlobalDOFManager([2] * 6)
    pdm = PatchDOFManager(2, mapping, dof_mgr)
    patch = Patch(BSplineSurface(su, sv), mgr, mapping, [3, 2], pdm)

    for _ in range(refine_u):
        SubdivisionRefiner(direction=0, n_levels=1).refine(patch)
    for _ in range(refine_v):
        SubdivisionRefiner(direction=1, n_levels=1).refine(patch)

    return patch, mgr, r1, r2


def test_nurbs_stiffness_symmetry():
    """NURBS quarter-ring stiffness matrix must be symmetric."""
    patch, mgr, _, _ = _build_quarter_ring()
    basis_u = IGABasis1D.build(patch.tensor.components[0], 4)
    basis_v = IGABasis1D.build(patch.tensor.components[1], 3)
    mat = PlaneStress(Material(E=210000, nu=0.3))
    K = PatchIntegrator(patch, basis_u, basis_v, mat).integrate_stiffness()
    K_arr = K.toarray()
    np.testing.assert_allclose(K_arr, K_arr.T, atol=1e-10,
        err_msg="NURBS stiffness matrix is not symmetric")


def test_nurbs_mass_symmetry():
    """NURBS quarter-ring mass matrix must be symmetric."""
    patch, mgr, _, _ = _build_quarter_ring()
    basis_u = IGABasis1D.build(patch.tensor.components[0], 4)
    basis_v = IGABasis1D.build(patch.tensor.components[1], 3)
    mat = PlaneStress(Material(E=210000, nu=0.3, rho=1.0))
    M = PatchIntegrator(patch, basis_u, basis_v, mat).integrate_mass()
    M_arr = M.toarray()
    np.testing.assert_allclose(M_arr, M_arr.T, atol=1e-10,
        err_msg="NURBS mass matrix is not symmetric")


def test_nurbs_mass_area():
    """
    NURBS bases form a partition of unity (sum_a R_a = 1), so:
      sum(M) / 2 / rho = integral dΩ = Area of the quarter ring.
    Quarter ring area = pi * (r2^2 - r1^2) / 4 = 3*pi/4.
    Uses a refined mesh for accurate Gauss quadrature.
    """
    patch, mgr, r1, r2 = _build_quarter_ring(refine_u=3, refine_v=2)
    basis_u = IGABasis1D.build(patch.tensor.components[0], 4)
    basis_v = IGABasis1D.build(patch.tensor.components[1], 3)
    rho = 7800.0
    mat = PlaneStress(Material(E=210000, nu=0.3, rho=rho))
    M = PatchIntegrator(patch, basis_u, basis_v, mat).integrate_mass()

    area_computed = np.sum(M.toarray()) / 2.0 / rho
    area_exact = np.pi * (r2**2 - r1**2) / 4.0   # 3*pi/4 ≈ 2.3562

    np.testing.assert_allclose(area_computed, area_exact, rtol=1e-4,
        err_msg=f"NURBS mass area check: got {area_computed:.8f}, expected {area_exact:.8f}")


def test_nurbs_stiffness_differs_from_bspline():
    """
    K_NURBS must differ from K_B-spline on the same control-point network.
    If the NURBS code path were not taken, the two would be equal.
    """
    patch_nurbs, mgr_nurbs, _, _ = _build_quarter_ring()

    # Build identical B-spline patch (same CPs, but no weights)
    mgr_bs = ControlPointManager(dim=2)
    for i in range(mgr_nurbs.n_points):
        mgr_bs.add_point(list(patch_nurbs.control_point(i)))   # w defaults to 1.0
    assert not mgr_bs.is_rational

    su = patch_nurbs.tensor.components[0]
    sv = patch_nurbs.tensor.components[1]
    mapping = list(range(mgr_bs.n_points))
    dof_mgr = GlobalDOFManager([2] * mgr_bs.n_points)
    pdm = PatchDOFManager(2, mapping, dof_mgr)
    patch_bs = Patch(BSplineSurface(su, sv), mgr_bs, mapping,
                     list(patch_nurbs.local_shape), pdm)

    basis_u = IGABasis1D.build(su, 4)
    basis_v = IGABasis1D.build(sv, 3)
    mat = PlaneStress(Material(E=210000, nu=0.3))

    K_nurbs = PatchIntegrator(patch_nurbs, basis_u, basis_v, mat).integrate_stiffness()
    K_bs    = PatchIntegrator(patch_bs,    basis_u, basis_v, mat).integrate_stiffness()

    diff = np.max(np.abs(K_nurbs.toarray() - K_bs.toarray()))
    assert diff > 1e-4, (
        f"NURBS and B-spline stiffness matrices are too similar (max diff={diff:.2e}); "
        "NURBS rationalization may not be active")


if __name__ == '__main__':
    test_bspline_is_not_rational()
    test_add_point_weighted_activates_rational()
    test_set_weight_activates_and_reverts()
    test_uniform_weights_reproduce_bspline_stiffness()
    test_uniform_weights_reproduce_bspline_mass()
    test_quarter_circle_geometry()
    test_nurbs_knot_insertion_preserves_quarter_circle()
    test_nurbs_degree_elevation_preserves_quarter_circle()
    test_nurbs_stiffness_symmetry()
    test_nurbs_mass_symmetry()
    test_nurbs_mass_area()
    test_nurbs_stiffness_differs_from_bspline()
    print("All NURBS tests passed.")
