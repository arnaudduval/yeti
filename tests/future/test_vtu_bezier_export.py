"""
Regression tests for write_bezier_patch_vtu()'s 3D (hexahedron) export path.

Two distinct, compounding VTK/writer bugs were found and fixed here (see the
comment block above _perm_3d() in vtu.py for the full diagnosis):

1. _perm_3d()'s face-interior point group was ordered wrong (an unverified
   guess at VTK's VTK_BEZIER_HEXAHEDRON convention) -- fixed by reading the
   real order directly off vtkBezierHexahedron.PointIndexFromIJK.
2. The hand-rolled XML declared `<VTKFile ... version="2.0">`, which VTK's
   reader silently mis-parses VTK_BEZIER_HEXAHEDRON connectivity under --
   fixed by declaring version "2.2" (matching what vtkXMLUnstructuredGridWriter
   itself emits). This was also the actual root cause of a previously
   suspected "anisotropic HigherOrderDegrees misread" crash -- confirmed by
   writing an anisotropic-degree (2,2,1) patch directly (no degree-elevation
   workaround) and finding both the correct order is read back AND
   warp_by_vector()/clip() no longer segfault.

The tests below check the *actual* connectivity read back by VTK against
ground truth from vtkBezierHexahedron.PointIndexFromIJK, for several
isotropic and anisotropic degree combinations -- the only way earlier tests
in this file could have caught either bug (the previous version only
checked the pre-export NURBS clone's own geometry, never anything VTK itself
parsed).

Every patch above is a SINGLE Bezier element (open knot vectors with no
internal knots), which never exercises `_local_flat_to_global()` -- the
function mapping a local per-element point to the shared global connectivity
ID across multiple elements (it degenerates to an identity map for one
element). The multi-element tests below close that gap: they check
watertightness (zero non-manifold "crack" edges, the exact shared-point
count a crack-free shared connectivity implies) rather than hand-deriving
per-point expected coordinates, which is easy to get subtly wrong for
non-trivial element layouts.
"""

import numpy as np
import pytest

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import (
    BSpline, BSplineVolume, Patch, ControlPointManager, SubdivisionRefiner,
)
from yeti_iga.future.postprocessing.vtu import write_bezier_patch_vtu

W0 = 1.0 / np.sqrt(2.0)


def _build_shell_octant_patch():
    """Anisotropic-degree (2,2,1) NURBS geometry (octant of a hollow
    sphere), same construction as test_nurbs_solid_3d.py."""
    def profile(R):
        return [(0.0, R, 1.0), (R, R, W0), (R, 0.0, 1.0)]
    rot = [(1.0, 0.0, 1.0), (1.0, 1.0, W0), (0.0, 1.0, 1.0)]

    mgr = ControlPointManager(dim=3)
    n_cp = 3 * 3 * 2
    for radius in (1.0, 2.0):
        prof = profile(radius)
        for iv in range(3):
            rx, ry, wv = rot[iv]
            for iu in range(3):
                x_p, z_p, wu = prof[iu]
                mgr.add_point([x_p * rx, x_p * ry, z_p], w=wu * wv)

    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sw = BSpline(1, np.array([0., 0., 1., 1.]))
    tensor = BSplineVolume(su, sv, sw)
    return Patch(tensor, mgr, list(range(n_cp)), [3, 3, 2])


def _build_box_patch(pu, pv, pw):
    """Trivial non-rational box on a single Bezier element: degree (pu,pv,pw),
    control points sitting exactly on an integer (iu,iv,iw) grid -- the
    simplest geometry for which "is the exported connectivity/point order
    correct" has an unambiguous, hand-checkable answer."""
    mgr = ControlPointManager(dim=3)
    n_cp = (pu + 1) * (pv + 1) * (pw + 1)
    for iw in range(pw + 1):
        for iv in range(pv + 1):
            for iu in range(pu + 1):
                mgr.add_point([float(iu), float(iv), float(iw)])
    kv = lambda p: np.array([0.0] * (p + 1) + [1.0] * (p + 1))
    tensor = BSplineVolume(BSpline(pu, kv(pu)), BSpline(pv, kv(pv)), BSpline(pw, kv(pw)))
    return Patch(tensor, mgr, list(range(n_cp)), [pu + 1, pv + 1, pw + 1])


def _vtk_ground_truth_order(pu, pv, pw):
    """index -> (i,j,k) exactly as VTK's own vtkBezierHexahedron assigns it."""
    vtk = pytest.importorskip("vtk")
    order = [pu, pv, pw]
    idx_to_ijk = [None] * ((pu + 1) * (pv + 1) * (pw + 1))
    for i in range(pu + 1):
        for j in range(pv + 1):
            for k in range(pw + 1):
                idx = vtk.vtkBezierHexahedron.PointIndexFromIJK(i, j, k, order)
                idx_to_ijk[idx] = (i, j, k)
    return idx_to_ijk


@pytest.mark.parametrize("pu,pv,pw", [(2, 2, 2), (3, 3, 3), (2, 2, 1), (1, 1, 1), (3, 2, 1)])
def test_exported_box_connectivity_matches_vtk_point_order(tmp_path, pu, pv, pw):
    """
    For a trivial box whose control points sit exactly on the (iu,iv,iw)
    grid, the physical coordinates ARE the grid indices -- so reading the
    exported file's points back in VTK's own cell-point order must give
    exactly the (i,j,k) triples vtkBezierHexahedron itself expects at each
    rank. This is the direct, VTK-level check that catches both the
    face-ordering bug and the file-version bug (either one scrambles this).
    """
    vtk = pytest.importorskip("vtk")
    patch = _build_box_patch(pu, pv, pw)
    out = tmp_path / "box.vtu"
    write_bezier_patch_vtu(patch, out)

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(out))
    reader.Update()
    cell = reader.GetOutput().GetCell(0)
    assert cell.GetNumberOfPoints() == (pu + 1) * (pv + 1) * (pw + 1)

    truth = _vtk_ground_truth_order(pu, pv, pw)
    pts = cell.GetPoints()
    for rank in range(cell.GetNumberOfPoints()):
        got = pts.GetPoint(rank)
        np.testing.assert_allclose(got, truth[rank], err_msg=f"rank={rank}")


def test_exported_anisotropic_hexahedron_order_read_back_correctly(tmp_path):
    """Anisotropic degree (2,2,1): GetOrder() must read back the real,
    per-direction degree (not a duplicated/misread one)."""
    vtk = pytest.importorskip("vtk")
    patch = _build_shell_octant_patch()
    out = tmp_path / "octant.vtu"
    write_bezier_patch_vtu(patch, out)

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(out))
    reader.Update()
    cell = reader.GetOutput().GetCell(0)
    assert [cell.GetOrder(i) for i in range(3)] == [2, 2, 1]


def test_write_bezier_patch_vtu_does_not_mutate_patch(tmp_path):
    patch = _build_shell_octant_patch()
    orig_n_cp = patch.n_cp

    field = np.random.default_rng(0).normal(size=(orig_n_cp, 3))
    write_bezier_patch_vtu(patch, tmp_path / "octant.vtu", field=field, field_name="disp")

    assert patch.n_cp == orig_n_cp


def test_exported_anisotropic_hexahedron_survives_pyvista_decomposition(tmp_path):
    """
    End-to-end regression test for the original VTK crash: write the
    anisotropic-degree NURBS solid directly (no workaround), then run the
    two pyvista operations that used to segfault (warp_by_vector, clip).
    """
    pv = pytest.importorskip("pyvista")

    patch = _build_shell_octant_patch()
    field = np.random.default_rng(0).normal(size=(patch.n_cp, 3))
    out = tmp_path / "octant.vtu"
    write_bezier_patch_vtu(patch, out, field=field, field_name="disp")

    mesh = pv.read(str(out))
    warped = mesh.warp_by_vector("disp", factor=0.3)
    assert warped.n_points == mesh.n_points
    clipped = mesh.clip(normal="z", origin=(0.0, 0.0, 1.0))
    assert clipped.n_cells > 0


def test_exported_box_multi_element_is_watertight(tmp_path):
    """
    Subdivide the trivial degree-2 box into 2x2x2 = 8 elements before
    exporting -- the simplest case that exercises shared connectivity across
    elements (_local_flat_to_global's non-trivial, non-identity path). A
    crack-free shared connectivity implies an exact, predictable number of
    unique points (5 per direction: 2 elements * degree 2 + 1), zero
    non-manifold edges on the rendered surface, zero open boundary edges
    (the box is closed), and the correct total surface area/bounds -- any
    one of these would fail if elements silently duplicated their shared
    boundary points instead of referencing the same global IDs.
    """
    pv = pytest.importorskip("pyvista")
    patch = _build_box_patch(2, 2, 2)
    SubdivisionRefiner(direction=0, n_levels=1).refine(patch)
    SubdivisionRefiner(direction=1, n_levels=1).refine(patch)
    SubdivisionRefiner(direction=2, n_levels=1).refine(patch)

    out = tmp_path / "box_multi.vtu"
    write_bezier_patch_vtu(patch, out)

    mesh = pv.read(str(out))
    assert mesh.n_cells == 8
    assert mesh.n_points == 5 ** 3

    surf = mesh.extract_surface(nonlinear_subdivision=3, algorithm=None)
    non_manifold = surf.extract_feature_edges(
        boundary_edges=False, feature_edges=False,
        manifold_edges=False, non_manifold_edges=True)
    assert non_manifold.n_cells == 0

    open_edges = surf.extract_feature_edges(
        boundary_edges=True, feature_edges=False,
        manifold_edges=False, non_manifold_edges=False)
    assert open_edges.n_cells == 0

    sizes = surf.compute_cell_sizes(area=True)
    np.testing.assert_allclose(sizes["Area"].sum(), 24.0, atol=1.e-10)
    np.testing.assert_allclose(surf.bounds, (0, 2, 0, 2, 0, 2))


def test_exported_nurbs_multi_element_is_watertight(tmp_path):
    """
    Same shared-connectivity check as test_exported_box_multi_element_is_watertight,
    on the real NURBS octant-of-a-hollow-sphere geometry subdivided into
    2x2x2 = 8 elements -- confirms multi-element export also works for a
    rational patch, not only the trivial affine box.
    """
    pv = pytest.importorskip("pyvista")
    patch = _build_shell_octant_patch()
    SubdivisionRefiner(direction=0, n_levels=1).refine(patch)
    SubdivisionRefiner(direction=1, n_levels=1).refine(patch)
    SubdivisionRefiner(direction=2, n_levels=1).refine(patch)

    out = tmp_path / "octant_multi.vtu"
    write_bezier_patch_vtu(patch, out)

    mesh = pv.read(str(out))
    assert mesh.n_cells == 8
    assert mesh.n_points == 5 * 5 * 3  # (2*2+1) * (2*2+1) * (2*1+1)

    surf = mesh.extract_surface(nonlinear_subdivision=4, algorithm=None)
    non_manifold = surf.extract_feature_edges(
        boundary_edges=False, feature_edges=False,
        manifold_edges=False, non_manifold_edges=True)
    assert non_manifold.n_cells == 0


if __name__ == '__main__':
    import pathlib
    tmp = pathlib.Path('/tmp')
    for pu, pv, pw in [(2, 2, 2), (3, 3, 3), (2, 2, 1), (1, 1, 1), (3, 2, 1)]:
        test_exported_box_connectivity_matches_vtk_point_order(tmp, pu, pv, pw)
    test_exported_anisotropic_hexahedron_order_read_back_correctly(tmp)
    test_write_bezier_patch_vtu_does_not_mutate_patch(tmp)
    test_exported_anisotropic_hexahedron_survives_pyvista_decomposition(tmp)
    test_exported_box_multi_element_is_watertight(tmp)
    test_exported_nurbs_multi_element_is_watertight(tmp)
