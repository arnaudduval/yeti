"""
Test Patch's helper methods on a genuine 3D (BSplineVolume) patch:
control_points_for_span() and boundary_control_points() -- the generic N-D
rewrite described in the 3D-solid plan (Patch.active_local_flat_indices(),
shared by both, plus PatchIntegrator::buildSpanLocalIndices()).

weights_for_span() is not exposed to Python (internal to the C++ NURBS
assembly path) -- not tested directly here; 3D NURBS is out of scope for
this pass anyway (see test_stiffness_3d.py's rejection test).
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import BSpline, BSplineVolume, ControlPointManager, Patch


def _build_1elt_lin_cube_patch():
    """
    Single-element (degree=1, one knot span per direction) unit cube, 8
    control points at the corners, u-fastest local numbering
    (local_shape=[2,2,2] -> flat = iu + 2*iv + 4*iw).
    """
    mgr = ControlPointManager(dim=3)
    for iw in range(2):
        for iv in range(2):
            for iu in range(2):
                mgr.add_point([float(iu), float(iv), float(iw)])

    kv = np.array([0., 0., 1., 1.])
    su, sv, sw = BSpline(1, kv), BSpline(1, kv), BSpline(1, kv)
    tensor = BSplineVolume(su, sv, sw)
    return Patch(tensor, mgr, list(range(8)), [2, 2, 2])


def _build_2elt_lin_box_patch():
    """
    Two elements along u only (degree=1), one element along v, w: a 3x2x2
    box of control points, local_shape=[3,2,2] -> flat = iu + 3*iv + 6*iw.
    """
    mgr = ControlPointManager(dim=3)
    for iw in range(2):
        for iv in range(2):
            for iu in range(3):
                mgr.add_point([float(iu), float(iv), float(iw)])

    kv_u = np.array([0., 0., 0.5, 1., 1.])
    kv_vw = np.array([0., 0., 1., 1.])
    su, sv, sw = BSpline(1, kv_u), BSpline(1, kv_vw), BSpline(1, kv_vw)
    tensor = BSplineVolume(su, sv, sw)
    return Patch(tensor, mgr, list(range(12)), [3, 2, 2])


def test_control_points_for_span_3d_single_element():
    patch = _build_1elt_lin_cube_patch()
    spans = list(patch.spans())
    assert len(spans) == 1  # single element in each of the 3 directions

    pts = patch.control_points_for_span(spans[0])

    assert pts.shape == (8, 3)
    expected = [(iu, iv, iw) for iw in range(2) for iv in range(2) for iu in range(2)]
    np.testing.assert_allclose(pts, np.array(expected, dtype=float))


def test_control_points_for_span_3d_two_elements():
    patch = _build_2elt_lin_box_patch()
    spans = list(patch.spans())
    assert len(spans) == 2  # 2 elements along u, 1 along v and w

    # u-fastest span order: spans[0] is u in [0, 0.5], spans[1] is u in [0.5, 1].
    pts_first = patch.control_points_for_span(spans[0])
    expected_first = [(iu, iv, iw) for iw in range(2) for iv in range(2) for iu in (0, 1)]
    np.testing.assert_allclose(pts_first, np.array(expected_first, dtype=float))

    pts_second = patch.control_points_for_span(spans[1])
    expected_second = [(iu, iv, iw) for iw in range(2) for iv in range(2) for iu in (1, 2)]
    np.testing.assert_allclose(pts_second, np.array(expected_second, dtype=float))


def test_boundary_control_points_3d_faces():
    patch = _build_1elt_lin_cube_patch()

    # local_shape = [2, 2, 2], flat = iu + 2*iv + 4*iw.
    assert sorted(patch.boundary_control_points(0, 0)) == [0, 2, 4, 6]  # x=0 face
    assert sorted(patch.boundary_control_points(0, 1)) == [1, 3, 5, 7]  # x=1 face
    assert sorted(patch.boundary_control_points(1, 0)) == [0, 1, 4, 5]  # y=0 face
    assert sorted(patch.boundary_control_points(1, 1)) == [2, 3, 6, 7]  # y=1 face
    assert sorted(patch.boundary_control_points(2, 0)) == [0, 1, 2, 3]  # z=0 face
    assert sorted(patch.boundary_control_points(2, 1)) == [4, 5, 6, 7]  # z=1 face


def test_boundary_control_points_3d_face_is_whole_face_not_single_edge():
    """A 3D boundary is a full 2D face (4 CPs for a linear cube), not a
    single edge as it would degenerate to if the 2D-only `varying = 1 -
    direction` formula had been left unchanged."""
    patch = _build_1elt_lin_cube_patch()
    assert len(patch.boundary_control_points(0, 0)) == 4


def test_boundary_control_points_3d_span_subrange_not_implemented():
    """3D sub-range selection (span_min/span_max >= 0) is explicitly out of
    scope for this pass -- must raise, not silently return a wrong subset."""
    patch = _build_2elt_lin_box_patch()
    try:
        patch.boundary_control_points(1, 0, span_min=0, span_max=0)
        assert False, "expected an exception for a 3D span sub-range request"
    except Exception:
        pass


if __name__ == '__main__':
    test_control_points_for_span_3d_single_element()
    test_control_points_for_span_3d_two_elements()
    test_boundary_control_points_3d_faces()
    test_boundary_control_points_3d_face_is_whole_face_not_single_edge()
    test_boundary_control_points_3d_span_subrange_not_implemented()
