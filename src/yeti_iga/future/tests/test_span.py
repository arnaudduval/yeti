import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "build"))


import numpy as np
from bspline import BSpline, BSplineSurface, BSplineVolume, ControlPointManager, Patch
import matplotlib
import matplotlib.pyplot as plt

U = np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])
V = np.array([0., 0., 0., 0.2, 0.4, 0.6, 0.8, 1., 1., 1.])
W = np.array([0., 0., 0., 0., 0.3, 0.4, 0.8, 0.9, 1., 1., 1., 1.])

def test_BSpline_getters():
    b1 = BSpline(2, U)
    assert (b1.knot_vector == np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])).all()
    assert b1.degree == 2

def test_ND_BSpline():
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
    assert np.allclose(patch.local_cp_view(0), np.array([0., 0.]), rtol=1.e-9)
    assert np.allclose(patch.local_cp_view(2), np.array([0., 1.]), rtol=1.e-9)

    # View to CP coordinates using mapping (zero-copy)
    assert np.allclose(patch.local_control_point_view()[mapping, :],
                       np.array([[0., 0.],[1., 0.],[0., 1.],[1., 1.]]),
                       rtol=1.e-9)

def test_evaluation():
    mgr = ControlPointManager(dim=2)

    id0 = mgr.add_point([0.0, 0.0])
    id1 = mgr.add_point([1.0, 0.0])
    id2 = mgr.add_point([0.0, 1.0])
    id3 = mgr.add_point([1.0, 1.0])
    # TODO : les ajouter dans un autre ordre et avoir un mapping non ordonné

    # Create surface with degree 1
    su = BSpline(1, np.array([0., 0., 1., 1.]))
    sv = BSpline(1, np.array([0., 0., 1., 1.]))
    surf = BSplineSurface(su, sv)
    mapping = np.array([0, 1, 2, 3], dtype=np.int64)
    local_shape = [2, 2]
    patch = Patch(surf, mgr, mapping.tolist(), local_shape)

    u = np.array([0.25, 0.75], dtype=np.float64)
    spans = surf.find_span_nd(u)
    tensor_basis = surf.basis_funs_nd(spans, u)

    print(spans)
    print(tensor_basis)


    local_pts = patch.local_control_point_view()[mapping, :]
    dim_u, dim_v = tensor_basis.shape

    print(local_pts)
    print(dim_u, dim_v)

    local_pts_reshaped = local_pts.reshape((dim_u, dim_v, mgr.dim_phys))
    print(f'{local_pts_reshaped[0, 1, :] = }')
    print(f'{local_pts[1, :] = }')
    surface_pt = np.tensordot(tensor_basis, local_pts_reshaped, axes=([1,0],[0,1]))

    print(local_pts_reshaped)
    print(surface_pt)

    # Pb inversion x, y
    # assert np.allclose(surface_pt, u, rtol=1.e-9)





test_BSpline_getters()
test_ND_BSpline()
test_cp_manager()
test_evaluation()


exit()



# Creation d'un espace 2D



print("========")
# ERROR : IL Y A UNE INVERSION DE DIRECTION DANS CE QUI EST FAIT CI-DESSOUS !!!





print(f"Surface evaluated at u={u[0]}, v={u[1]} : {surface_pt}")



print("======")
u = np.array([[0.25, 0.75],[0.7, 0.66]], dtype=np.float64)
spans = np.array([surf.find_span_nd(u[0, :]), surf.find_span_nd(u[1])])
print(spans)

print(patch.evaluate_patch_nd(spans, u))

