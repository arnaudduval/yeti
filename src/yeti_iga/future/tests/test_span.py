import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "build"))


import numpy as np
from bspline import BSpline, BSplineSurface, BSplineVolume, ControlPointManager, Patch
from bspline import IGABasis1D, PatchDOFManager, GlobalDOFManager, PatchIntegrator#, IGAAssembler2D
import matplotlib
import matplotlib.pyplot as plt

U = np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])
V = np.array([0., 0., 0., 0.2, 0.4, 0.6, 0.8, 1., 1., 1.])
W = np.array([0., 0., 0., 0., 0.3, 0.4, 0.8, 0.9, 1., 1., 1., 1.])

def test_BSpline_getters():
    """
    Test getters of knot vector and degree for a 1D BSpline parametric space
    """
    b1 = BSpline(2, U)
    assert (b1.knot_vector == np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])).all()
    assert b1.degree == 2

def test_ND_BSpline():
    """
    Test build of BSpline 2D and 3D poarametric space by tensor product of 1D BSplines
    Tets span search and function computation on 2D and 3D space
    """
    b1 = BSpline(2, U)
    b2 = BSpline(2, V)
    b3 = BSpline(3, W)
    surf = BSplineSurface(b1, b2)
    vol = BSplineVolume(b1, b2, b3)

    u = np.array([0.3, 0.45])
    span = surf.find_span_nd(u)
    assert (span == [2, 4]).all()
    ref_2D = np.array([[2.32438017e-03, 5.68181818e-03, 2.58264463e-04],
                       [1.62706612e-01, 3.97727273e-01, 1.80785124e-02],
                       [1.16219008e-01, 2.84090909e-01, 1.29132231e-02]])

    funs = surf.basis_funs_nd(span, u)
    assert np.allclose(np.sum(funs), 1.0, rtol = 1.e-9)
    assert np.allclose(funs, ref_2D, rtol = 1.e-9)


    u = np.array([0.3, 0.72, 0.45])
    span = vol.find_span_nd(u)
    assert( span == [2, 5, 5]).all()
    ref_3D = np.array([[[1.77169421e-04, 4.13739669e-04, 6.95592287e-05, 6.88705234e-07],
                        [1.63881715e-03, 3.82709194e-03, 6.43422865e-04, 6.37052342e-06],
                        [3.98631198e-04, 9.30914256e-04, 1.56508264e-04, 1.54958678e-06]],
                       [[1.24018595e-02, 2.89617769e-02, 4.86914601e-03, 4.82093664e-05],
                        [1.14717200e-01, 2.67896436e-01, 4.50396006e-02, 4.45936639e-04],
                        [2.79041839e-02, 6.51639979e-02, 1.09555785e-02, 1.08471074e-04]],
                       [[8.85847107e-03, 2.06869835e-02, 3.47796143e-03, 3.44352617e-05],
                        [8.19408574e-02, 1.91354597e-01, 3.21711433e-02, 3.18526171e-04],
                        [1.99315599e-02, 4.65457128e-02, 7.82541322e-03, 7.74793388e-05]]])

    funs = vol.basis_funs_nd(span, u)
    # print(funs.shape)
    assert np.allclose(np.sum(funs), 1.0, rtol=1.e-9)
    assert np.allclose(funs, ref_3D, rtol=1.e-9)


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
    # Test view to global coordinates (zero-copy)
    assert np.allclose(mgr.coords_view(), np.array([[0., 0.],[1., 0.],[1., 1.],[0., 1.]]), rtol=1.e-9)

    # Create surface with degree 1
    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = np.array([0, 1, 3, 2], dtype=np.int64)
    local_shape = [2, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    # Test view to local control point
    assert np.allclose(patch.control_point(0), np.array([0., 0.]), rtol=1.e-9)
    assert np.allclose(patch.control_point(2), np.array([0., 1.]), rtol=1.e-9)
    # View to CP coordinates using mapping (zero-copy)
    assert np.allclose(patch.local_control_point_view()[mapping, :],
                       np.array([[0., 0.],[1., 0.],[0., 1.],[1., 1.]]),
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
    mapping = np.array([0, 4, 1, 2, 5, 3], dtype=np.int64)
    local_shape = [3, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    u = np.array([0.25, 0.75], dtype=np.float64)
    spans = surf.find_span_nd(u)
    tensor_basis = surf.basis_funs_nd(spans, u)

    local_pts = patch.local_control_point_view()[mapping, :]
    dim_u, dim_v = tensor_basis.shape

    local_pts_reshaped = local_pts.reshape((dim_v, dim_u, mgr.dim_phys)).transpose(1, 0, 2)
    surface_pt = np.tensordot(tensor_basis, local_pts_reshaped, axes=([0,1],[0,1]))

    # test single evaluation
    assert np.allclose(surface_pt, u*[3., 1.], rtol=1.e-9)
    assert np.allclose(patch.evaluate_patch_nd([spans], [u]),  [u*[3., 1.]], rtol=1.e-9)
    assert np.allclose(patch.evaluate_patch_nd_omp([spans], [u]),  [u*[3., 1.]], rtol=1.e-9)


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
    Test multiple evaluation on a Patch when knot vector contains repeated inner knots
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
    mapping = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.int64)
    local_shape = [5, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    u = np.array([0.5, 1.0])

    # In 2D, this patch has only 2 spans (2 in u, 1 in v)
    ref_spans = np.array([[2, 1], [4, 1]])
    for i, span in enumerate(patch.spans()):
        assert (span == ref_spans[i]).all()

    # get CPs for a given span
    pts = patch.control_points_for_span(np.array([4, 1]))
    # pts = pts.reshape([3, 2, 2])
    # pts is return as stored in memory. It must be set in proper order
    pts = pts.reshape([2, 3, 2]).transpose(1, 0, 2)

    assert (pts[1, 1] == [2.25, 1.]).all()
    assert (pts[2, 1] == [3., 1.]).all()


def test_functions_derivatives():
    """
    Test functions 1st derivative, compared with reference finite differences value
    """
    b = BSpline(3, W)
    u_list = [0.1, 1./3., 0.5, 8./9.]

    eps = 1.e-7
    for u in u_list:
        span = b.find_span(u)
        funs = b.basis_funs(span, u)
        dfuns = b.basis_funs_derivatives(span, u, 2)

        assert np.allclose(dfuns[0, :], funs, rtol = 1.e-9)

        val_minus = b.basis_funs(span, u - eps)
        val_plus = b.basis_funs(span, u + eps)

        assert np.allclose(dfuns[1, :], (val_plus - val_minus)/(2.*eps), rtol = 1.e-6)


def test_coupling_strong():
    """
    2 patchs with string coupling (common control points)
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

    # mapping patch 1
    mapping1 = [0, 1, 3, 4]
    local_shape1 = [2, 2]
    dof_manager1 = PatchDOFManager(dofs_per_control_point=2, control_points=mapping1, global_dof_manager=global_dof_manager)
    # mapping patch 2
    mapping2 = [1, 2, 4, 5]
    local_shape2 = [2, 2]
    dof_manager2 = PatchDOFManager(dofs_per_control_point=2, control_points=mapping2, global_dof_manager=global_dof_manager)
    # mapping patch 3
    mapping3 = [6, 7, 8, 9]
    local_shape3 = [2, 2]
    dof_manager3 = PatchDOFManager(dofs_per_control_point=1, control_points=mapping3, global_dof_manager=global_dof_manager)


    patch1 = Patch(surf1, mgr, mapping1, local_shape1, dof_manager1)
    patch2 = Patch(surf2, mgr, mapping2, local_shape2, dof_manager2)
    patch3 = Patch(surf2, mgr, mapping3, local_shape3, dof_manager3)

    assert patch1.dof_manager.get_global_dof_indices(1) == [2, 3]
    assert patch1.dof_manager.get_global_dof_indices(2) == [6, 7]

    assert patch2.dof_manager.get_global_dof_indices(0) == [2,3]
    assert patch2.dof_manager.get_global_dof_indices(3) == [10, 11]

    assert patch3.dof_manager.get_global_dof_indices(0) == [12]
    assert patch3.dof_manager.get_global_dof_indices(3) == [15]


def test_integration_1elt_lin_square_1():
    """
    Test Gauss integration over a single patch, degree 1, dimension 1x1
    """

        # Compare with legacy YETI overlay
    from yeti_iga.preprocessing.igaparametrization import IGAparametrization
    from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
        as build_stiffmatrix
    import scipy.sparse as sp

    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/1_elt_lin')
    data, row, col, rhs = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())

    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    stiff_tot = stiff_side + stiff_side.transpose()

    print("==========")

    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([1., 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([1.0, 1.0])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    local_shape = [2, 2]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 2)   # 3 Gauss points per span
    basis_v = IGABasis1D.build(sv, 2)   # 2 Gauss points per span

    integrator = PatchIntegrator(patch, basis_u, basis_v)
    stiffness_matrix = integrator.integrate()

    print("stiffness_matrix", stiffness_matrix)
    print(stiff_tot)

    assert np.allclose(stiffness_matrix.toarray(), stiff_tot.toarray(), rtol=1.e-5, atol=1.e-8)

def test_integration_1elt_lin_rect():
    """
    Test Gauss integration over a single patch, degree 1, dimension 1x3
    """

        # Compare with legacy YETI overlay
    from yeti_iga.preprocessing.igaparametrization import IGAparametrization
    from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
        as build_stiffmatrix
    import scipy.sparse as sp

    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/1_elt_lin_rect')
    data, row, col, rhs = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())

    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    stiff_tot = stiff_side + stiff_side.transpose()

    print("==========")

    mgr = ControlPointManager(dim=2)
    mgr.add_point([0.0, 0.0])
    mgr.add_point([3., 0.0])
    mgr.add_point([0.0, 1.0])
    mgr.add_point([3.0, 1.0])

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3]
    local_shape = [2, 2]

    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 2)   # 3 Gauss points per span
    basis_v = IGABasis1D.build(sv, 2)   # 2 Gauss points per span

    integrator = PatchIntegrator(patch, basis_u, basis_v)
    stiffness_matrix = integrator.integrate()

    print("stiffness_matrix", stiffness_matrix)
    print(stiff_tot)

    assert np.allclose(stiffness_matrix.toarray(), stiff_tot.toarray(), rtol=1.e-5, atol=1.e-8)


def test_integration_2_elements_C0():
    """
    Test Gauss integration over a Patch
    WARNING : unfinished
    """

    # Basic demo patch
    # TODO Should be vectorized
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

    dofs_per_control_point = [2 for _ in range(mgr.n_points)]
    dof_manager = GlobalDOFManager(dofs_per_control_point)

    su = BSpline(2, np.array([0., 0., 0., 0.5, 0.5, 1., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    local_shape = [5, 2]


    dof_manager_patch = PatchDOFManager(2, mapping, dof_manager)

    patch = Patch(surf, mgr, mapping, local_shape, dof_manager_patch)

    # Create 1D integration basis
    basis_u = IGABasis1D.build(su, 3)   # 3 Gauss points per span
    basis_v = IGABasis1D.build(sv, 2)   # 2 Gauss points per span

    print(basis_u.gauss_spans)
    print(basis_v.gauss_spans)
    for sp in basis_u.gauss_spans:
        print(f"{sp.u_param = }")
        print(f"{sp.weight}")
        print(f"{sp.N}")
        print(f"{sp.dN}")

    print('---------')

    for sp in basis_v.gauss_spans:
        print(sp.u_param)
        print(sp.weight)
        print(sp.N)
        print(sp.dN)

    # Get global DOF indices of 1st control point
    # print(dof_manager_patch.get_global_dof_indices(0))
    # print(dof_manager_patch.get_global_dof_indices(1))
    # print(dof_manager_patch.get_global_dof_indices(2))
    # print(dof_manager_patch.get_global_dof_indices(3))
    # print(dof_manager_patch.get_global_dof_indices(4))
    # print(dof_manager_patch.get_global_dof_indices(5))
    # print(dof_manager_patch.get_global_dof_indices(6))
    # print(dof_manager_patch.get_global_dof_indices(7))
    # print(dof_manager_patch.get_global_dof_indices(8))
    # print(dof_manager_patch.get_global_dof_indices(9))

    integrator = PatchIntegrator(patch, basis_u, basis_v)
    stiffness_matrix = integrator.integrate()

    print("stiffness_matrix", stiffness_matrix)




def test_test_test():
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
    mapping = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.int64)
    local_shape = [5, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    # Create 1D basis
    basis_u = IGABasis1D.build(su, 3)   # 3 Gauss points per span
    basis_v = IGABasis1D.build(sv, 2)   # 2 Gauss points per span

    print(basis_u.spans)
    print(len(basis_u.spans))
    print(basis_u.spans[0].u_param)
    print(basis_u.spans[1].u_param)
    print(len(basis_v.spans))
    print(basis_v.spans[0].u_param)


    # Assembler
    assembler = IGAAssembler2D(patch, basis_u, basis_v)

    elems = assembler.assemble_stiffness()

    for elem in elems:
        print(elem.global_indices)
        print(elem.get_K_as_numpy())


    patch.test()


if __name__ == '__main__':
    test_integration_1elt_lin_rect()
    # test_integration_2_elements_C0()
    # test_test_test()
    exit()
    test_integration_1elt_lin_square_1()
    test_BSpline_getters()
    test_ND_BSpline()
    test_cp_manager()
    test_evaluation()
    test_evaluation_low_continuity()
    test_span_iterator()
    test_functions_derivatives()
    test_coupling_strong()
    print("All tests finshed !!!")

