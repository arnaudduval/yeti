"""
Test boundary control point selection (Patch.boundary_control_points) and
distributed boundary load assembly (PatchIntegrator.integrate_boundary_load /
assemble_boundary_load).
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (BSpline, BSplineSurface, ControlPointManager,
    Patch, GlobalDOFManager, PatchDOFManager, PatchAssembly, IGABasis1D,
    PatchIntegrator, MaterialProperties, ConstantTraction, BoundaryLoadSpec)


def _build_2_elts_c0_patch():
    """
    Same fixture/layout as test_stiffness.py::test_integration_2_elements_C0:
    a 6x1 rectangle, degree 2, two elements (C0 at the interior knot 0.5).
    """
    mgr = ControlPointManager(dim=2)
    for y in (0.0, 0.5, 1.0):
        for x in (0.0, 1.5, 3.0, 4.5, 6.0):
            mgr.add_point([x, y])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 0.5, 0.5, 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = list(range(15))
    local_shape = [5, 3]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)
    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    basis_u = IGABasis1D.build(su, 3)
    basis_v = IGABasis1D.build(sv, 3)
    return patch, basis_u, basis_v


def test_boundary_control_points_whole_edges():
    patch, _, _ = _build_2_elts_c0_patch()

    # local_shape = [5, 3], u-fastest flat index = iv*5 + iu.
    assert patch.boundary_control_points(0, 0) == [0, 5, 10]    # x=0 edge
    assert patch.boundary_control_points(0, 1) == [4, 9, 14]    # x=6 edge
    assert patch.boundary_control_points(1, 0) == [0, 1, 2, 3, 4]   # y=0 edge
    assert patch.boundary_control_points(1, 1) == [10, 11, 12, 13, 14]  # y=1 edge


def test_boundary_control_points_span_range():
    patch, _, _ = _build_2_elts_c0_patch()

    # Bottom edge (direction=1, side=0), degree-2 in u, two elements (knots
    # 0..0.5 and 0.5..1). The interior knot 0.5 has multiplicity 2 (C0), so
    # the raw span indices of the two VALID elements are 2 and 4 (span index
    # 3 is the zero-length [0.5, 0.5] span skipped by the knot vector's
    # double knot -- not a real element). Restricting to the first element
    # only must give local CP positions [0, 1, 2] (its local support); the
    # second element only must give [2, 3, 4].
    first_span_cps = patch.boundary_control_points(1, 0, span_min=2, span_max=2)
    second_span_cps = patch.boundary_control_points(1, 0, span_min=4, span_max=4)
    assert first_span_cps == [0, 1, 2]
    assert second_span_cps == [2, 3, 4]

    both_spans_cps = patch.boundary_control_points(1, 0, span_min=2, span_max=4)
    assert both_spans_cps == [0, 1, 2, 3, 4]


def test_integrate_boundary_load_resultant_force():
    """
    A constant traction (0, -p) on a straight edge of length L must produce
    a load vector whose y-components sum to -p*L (partition of unity: the
    active basis functions on the edge sum to 1 everywhere), x-components
    all zero, and zero entries away from the loaded edge.
    """
    patch, basis_u, basis_v = _build_2_elts_c0_patch()
    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))

    p = 1000.0
    traction = ConstantTraction(np.array([0.0, -p]))

    # Right edge (x=6, direction=0, side=1), length L=1.0 (y in [0, 1]).
    F = integrator.integrate_boundary_load(0, 1, traction)

    np.testing.assert_allclose(F[0::2].sum(), 0.0, atol=1.e-10)
    np.testing.assert_allclose(F[1::2].sum(), -p * 1.0, rtol=1.e-10)

    # Only control points on the right edge (global ids 4, 9, 14) carry a
    # nonzero contribution.
    loaded_cps = {4, 9, 14}
    for cp in range(patch.n_cp):
        dofs = patch.dof_manager.get_global_dof_indices(cp)
        if cp not in loaded_cps:
            assert F[dofs[0]] == 0.0
            assert F[dofs[1]] == 0.0


def test_integrate_boundary_load_span_range_matches_partial_resultant():
    """
    Restricting integrate_boundary_load() to one of the two spans of a
    straight edge of length L (made of 2 equal elements) must produce half
    the total resultant force of the whole edge.
    """
    patch, basis_u, basis_v = _build_2_elts_c0_patch()
    integrator = PatchIntegrator(patch, basis_u, basis_v, MaterialProperties(210000, 0.3))

    p = 1000.0
    traction = ConstantTraction(np.array([0.0, -p]))

    # Bottom edge (y=0, direction=1, side=0), length L=6.0, two equal
    # elements (x in [0,3] and [3,6]) -- raw span indices 2 and 4 (span 3 is
    # the zero-length span skipped by the double knot at 0.5), same knot
    # vector as test_boundary_control_points_span_range().
    F_first_half = integrator.integrate_boundary_load(1, 0, traction, span_min=2, span_max=2)
    F_second_half = integrator.integrate_boundary_load(1, 0, traction, span_min=4, span_max=4)
    F_whole = integrator.integrate_boundary_load(1, 0, traction)

    np.testing.assert_allclose(F_first_half[1::2].sum(), -p * 3.0, rtol=1.e-10)
    np.testing.assert_allclose(F_second_half[1::2].sum(), -p * 3.0, rtol=1.e-10)
    np.testing.assert_allclose(F_first_half + F_second_half, F_whole, rtol=1.e-10, atol=1.e-10)


def test_assemble_boundary_load_on_domain_split():
    """
    Splitting the same 6x1 rectangle into two patches sharing the middle
    edge (same construction as test_multipatch_stiffness.py) and loading
    only the right patch's far edge must give the same resultant force as
    loading the single-patch reference's far edge -- and zero contribution
    from the unloaded left patch.
    """
    mgr = ControlPointManager(dim=2)
    for y in (0.0, 0.5, 1.0):
        for x in (0.0, 1.5, 3.0, 4.5, 6.0):
            mgr.add_point([x, y])

    dofs_per_control_point = [2] * mgr.n_points
    global_dof_manager = GlobalDOFManager(dofs_per_control_point)

    su_left = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv_left = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    mapping_left = [0, 1, 2, 5, 6, 7, 10, 11, 12]
    dof_manager_left = PatchDOFManager(2, mapping_left, global_dof_manager)
    patch_left = Patch(BSplineSurface(su_left, sv_left), mgr, mapping_left, [3, 3],
                        dof_manager_left)

    su_right = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv_right = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    mapping_right = [2, 3, 4, 7, 8, 9, 12, 13, 14]
    dof_manager_right = PatchDOFManager(2, mapping_right, global_dof_manager)
    patch_right = Patch(BSplineSurface(su_right, sv_right), mgr, mapping_right, [3, 3],
                         dof_manager_right)

    assembly = PatchAssembly()
    assembly.add_patch(patch_left)
    assembly.add_patch(patch_right)
    assembly.detect_shared_control_points()

    p = 1000.0
    traction = ConstantTraction(np.array([0.0, -p]))
    spec = BoundaryLoadSpec(patch_index=1, direction=0, side=1, traction=traction)
    F = PatchIntegrator.assemble_boundary_load(assembly, [spec])

    n_dofs = global_dof_manager.get_dof_indices(mgr.n_points - 1)[-1] + 1
    assert F.shape == (n_dofs,)

    np.testing.assert_allclose(F[0::2].sum(), 0.0, atol=1.e-10)
    np.testing.assert_allclose(F[1::2].sum(), -p * 1.0, rtol=1.e-10)

    # Control points exclusive to the (unloaded) left patch carry no load.
    for cp in (0, 1, 5, 6, 10, 11):
        dofs = global_dof_manager.get_dof_indices(cp)
        assert F[dofs[0]] == 0.0
        assert F[dofs[1]] == 0.0


def test_assemble_boundary_load_rejects_out_of_range_patch_index():
    patch, _, _ = _build_2_elts_c0_patch()
    assembly = PatchAssembly()
    assembly.add_patch(patch)

    traction = ConstantTraction(np.array([0.0, -1.0]))
    spec = BoundaryLoadSpec(patch_index=5, direction=0, side=1, traction=traction)
    try:
        PatchIntegrator.assemble_boundary_load(assembly, [spec])
        assert False, "expected an exception for an out-of-range patch_index"
    except Exception:
        pass


if __name__ == '__main__':
    test_boundary_control_points_whole_edges()
    test_boundary_control_points_span_range()
    test_integrate_boundary_load_resultant_force()
    test_integrate_boundary_load_span_range_matches_partial_resultant()
    test_assemble_boundary_load_on_domain_split()
    test_assemble_boundary_load_rejects_out_of_range_patch_index()
