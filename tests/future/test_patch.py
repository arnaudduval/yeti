"""
Test BSpline patch manipulation :
 - control points manager
 - evaluation
 - iterator over spans
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import BSpline, BSplineSurface, \
    ControlPointManager, Patch, PatchDOFManager, GlobalDOFManager


def test_cp_manager():
    """
    Test control points manager :
     - build a CP manager
     - test getter of CP coordinates
     - build a Patch with 2D param space, CP manager and a connectivity table
     - test view on patch CP coordinates
    """

    mgr = ControlPointManager(dim=2)

    id0 = mgr.add_point([0.0, 0.0])
    id1 = mgr.add_point([1.0, 0.0])
    # id2 and id3 added deliberately in disorder
    id3 = mgr.add_point([1.0, 1.0])
    id2 = mgr.add_point([0.0, 1.0])

    assert mgr.n_points == 4
    assert (id0, id1, id3, id2) == (0, 1, 2, 3)
    # Test view to global coordinates (zero-copy)
    assert np.allclose(mgr.coords_view(),
                       np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]]),
                       rtol=1.e-9)

    # Create surface with degree 1
    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    # u-fastest mapping: flat 0=(u=0,v=0)→0, flat 1=(u=1,v=0)→1, flat 2=(u=0,v=1)→3, flat 3=(u=1,v=1)→2
    mapping = np.array([0, 1, 3, 2], dtype=np.int64)
    local_shape = [2, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    # Test view to local control point (u-fastest: flat 0=(0,0), flat 1=(1,0), flat 2=(0,1), flat 3=(1,1))
    assert np.allclose(patch.control_point(0), np.array([0., 0.]), rtol=1.e-9)
    assert np.allclose(patch.control_point(2), np.array([0., 1.]), rtol=1.e-9)
    # View to CP coordinates using mapping (zero-copy)
    assert np.allclose(patch.local_control_point_view()[mapping, :],
                       np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]]),
                       rtol=1.e-9)


def test_evaluation():
    """
    Test single and multiple evaluation of points on a Patch
    """
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([3.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([3.0, 1.0])
    mgr.add_point([1.5, 0.0])
    mgr.add_point([1.5, 1.0])

    # Create surface with degree 1 / 2
    su = BSpline(2, np.array([0., 0., 0., 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    # u-fastest mapping (nu=3, nv=2): flat = iu + iv*nu
    #   flat 0=(iu=0,iv=0)→(0,0)=mgr[0], flat 1=(iu=1,iv=0)→(1.5,0)=mgr[4]
    #   flat 2=(iu=2,iv=0)→(3,0)=mgr[1], flat 3=(iu=0,iv=1)→(0,1)=mgr[2]
    #   flat 4=(iu=1,iv=1)→(1.5,1)=mgr[5], flat 5=(iu=2,iv=1)→(3,1)=mgr[3]
    mapping = np.array([0, 4, 1, 2, 5, 3], dtype=np.int64)
    local_shape = [3, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    u = np.array([0.25, 0.75], dtype=np.float64)
    spans = surf.find_span_nd(u)
    tensor_basis = surf.basis_funs_nd(spans, u)

    local_pts = patch.local_control_point_view()[mapping, :]
    dim_u, dim_v = tensor_basis.shape

    # u-fastest: reshape as (nv, nu, dim_phys) then transpose to (nu, nv, dim_phys)
    local_pts_reshaped = local_pts.reshape((dim_v, dim_u, mgr.dim_phys)).transpose(1, 0, 2)
    surface_pt = np.tensordot(tensor_basis,
                              local_pts_reshaped,
                              axes=([0, 1], [0, 1]))

    # test single evaluation
    assert np.allclose(surface_pt, u*[3., 1.], rtol=1.e-9)
    assert np.allclose(patch.evaluate_patch_nd([spans], [u]),
                       [u*[3., 1.]],
                       rtol=1.e-9)
    assert np.allclose(patch.evaluate_patch_nd_omp([spans], [u]),
                       [u*[3., 1.]],
                       rtol=1.e-9)

    # Test multiple evaluation
    n_points = 100
    u = np.random.rand(n_points, 2)
    spans = np.array([surf.find_span_nd(pt) for pt in u])
    res_serial = patch.evaluate_patch_nd(spans, u)
    res_omp = patch.evaluate_patch_nd_omp(spans, u)

    assert np.allclose(res_omp, u*[3., 1.], rtol=1.e-9)
    assert np.allclose(res_serial, res_omp, rtol=1.e-9)


def test_evaluation_low_continuity():
    """
    Test multiple evaluations on a Patch when knot vector contains repeated
    inner knots
    """
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([0.75, 0.0])
    mgr.add_point([1.5, 0.0])
    mgr.add_point([2.25, 0.0])
    mgr.add_point([3.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([0.75, 1.0])
    mgr.add_point([1.5, 1.0])
    mgr.add_point([2.25, 1.0])
    mgr.add_point([3.0, 1.0])

    su = BSpline(2, np.array([0., 0., 0., 0.5, 0.5, 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    # u-fastest mapping (nu=5, nv=2): flat = iu + iv*nu
    #   flat 0=(iu=0,iv=0)→mgr[0], flat 1=(iu=1,iv=0)→mgr[1], ...
    #   flat 5=(iu=0,iv=1)→mgr[5], flat 6=(iu=1,iv=1)→mgr[6], ...
    mapping = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.int64)
    local_shape = [5, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    # Test multiple evaluation
    n_points = 100
    u = np.random.rand(n_points, 2)
    spans = np.array([surf.find_span_nd(pt) for pt in u])

    res_serial = patch.evaluate_patch_nd(spans, u)
    res_omp = patch.evaluate_patch_nd_omp(spans, u)

    assert np.allclose(res_omp, u*[3., 1.], rtol=1.e-9)
    assert np.allclose(res_serial, res_omp, rtol=1.e-9)


def test_span_iterator():
    """
    Test iterator on span of a 2D patch
    """
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([0.75, 0.0])
    mgr.add_point([1.5, 0.0])
    mgr.add_point([2.25, 0.0])
    mgr.add_point([3.0, 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([0.75, 1.0])
    mgr.add_point([1.5, 1.0])
    mgr.add_point([2.25, 1.0])
    mgr.add_point([3.0, 1.0])

    su = BSpline(2, np.array([0., 0., 0., 0.5, 0.5, 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    # u-fastest mapping (nu=5, nv=2): flat = iu + iv*nu
    #   flat 0=(iu=0,iv=0)→mgr[0], ..., flat 4=(iu=4,iv=0)→mgr[4]
    #   flat 5=(iu=0,iv=1)→mgr[5], ..., flat 9=(iu=4,iv=1)→mgr[9]
    mapping = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.int64)
    local_shape = [5, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    # In 2D, this patch has only 2 spans (2 in u, 1 in v)
    ref_spans = np.array([[2, 1], [4, 1]])
    for i, span in enumerate(patch.spans()):
        assert (span == ref_spans[i]).all()

    # get CPs for a given span (returned jv-outer, iu-inner)
    pts = patch.control_points_for_span(np.array([4, 1]))
    pts = pts.reshape([2, 3, 2]).transpose(1, 0, 2)

    assert (pts[1, 1] == [2.25, 1.]).all()
    assert (pts[2, 1] == [3., 1.]).all()


def test_coupling_strong():
    """
    2 patchs with strong coupling (common control points)
    1 extra patch without coupling but with different number of DOF/CP
    """
    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])    # index : 0
    mgr.add_point([4.0, 0.0])    # index : 1
    mgr.add_point([10.0, 0.0])    # index : 2
    mgr.add_point([0.0, 4.0])    # index : 3
    mgr.add_point([4.0, 4.0])    # index : 4
    mgr.add_point([10.0, 4.0])    # index : 5

    mgr.add_point([0.0, 0.0])   # index : 6
    mgr.add_point([10.0, 0.0])  # index : 7
    mgr.add_point([0.0, 4.0])   # index : 8
    mgr.add_point([10.0, 4.0])  # index : 9

    # parametric space patch 1
    su1 = BSpline(1, np.array([0., 0., 1., 1.]))
    sv1 = BSpline(1, np.array([0., 0., 1., 1.]))
    surf1 = BSplineSurface(su1, sv1)
    # parametric space patch2
    su2 = BSpline(1, np.array([0., 0., 1., 1.]))
    sv2 = BSpline(1, np.array([0., 0., 1., 1.]))
    surf2 = BSplineSurface(su2, sv2)
    # parametric space patch3
    su3 = BSpline(1, np.array([0., 0., 1., 1.]))
    sv3 = BSpline(1, np.array([0., 0., 1., 1.]))
    surf2 = BSplineSurface(su3, sv3)

    # create DOF manager
    dofs_per_control_point = [2, 2, 2, 2, 2, 2, 1, 1, 1, 1]
    global_dof_manager = GlobalDOFManager(dofs_per_control_point)

    # u-fastest mappings (nu=2, nv=2): flat = iu + iv*nu
    # patch 1: CPs (0,0),(4,0),(0,4),(4,4) → mgr indices 0,1,3,4
    mapping1 = [0, 1, 3, 4]
    local_shape1 = [2, 2]
    dof_manager1 = PatchDOFManager(dofs_per_control_point=2,
                                   control_points=mapping1,
                                   global_dof_manager=global_dof_manager)
    # patch 2: CPs (4,0),(10,0),(4,4),(10,4) → mgr indices 1,2,4,5
    mapping2 = [1, 2, 4, 5]
    local_shape2 = [2, 2]
    dof_manager2 = PatchDOFManager(dofs_per_control_point=2,
                                   control_points=mapping2,
                                   global_dof_manager=global_dof_manager)
    # patch 3: CPs (0,0),(10,0),(0,4),(10,4) → mgr indices 6,7,8,9
    mapping3 = [6, 7, 8, 9]
    local_shape3 = [2, 2]
    dof_manager3 = PatchDOFManager(dofs_per_control_point=1,
                                   control_points=mapping3,
                                   global_dof_manager=global_dof_manager)

    patch1 = Patch(surf1, mgr, mapping1, local_shape1, dof_manager1)
    patch2 = Patch(surf2, mgr, mapping2, local_shape2, dof_manager2)
    patch3 = Patch(surf2, mgr, mapping3, local_shape3, dof_manager3)

    # local CP 1 = global CP 1 (at (4,0)) → DOFs 2,3
    assert patch1.dof_manager.get_global_dof_indices(1) == [2, 3]
    # local CP 2 = global CP 3 (at (0,4)) → DOFs 6,7
    assert patch1.dof_manager.get_global_dof_indices(2) == [6, 7]

    assert patch2.dof_manager.get_global_dof_indices(0) == [2, 3]
    assert patch2.dof_manager.get_global_dof_indices(3) == [10, 11]

    assert patch3.dof_manager.get_global_dof_indices(0) == [12]
    assert patch3.dof_manager.get_global_dof_indices(3) == [15]


if __name__ == '__main__':
    test_cp_manager()
    test_evaluation()
    test_evaluation_low_continuity()
    test_span_iterator()
    test_coupling_strong()
