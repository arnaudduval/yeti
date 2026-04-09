"""
Basic test for BSpline parameter space and its tensor product
"""

import numpy as np

# pylint: disable=no-name-in-module
from yeti_iga.future.bspline import BSpline, BSplineSurface, BSplineVolume



# reference knot vectors
U = np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])
V = np.array([0., 0., 0., 0.2, 0.4, 0.6, 0.8, 1., 1., 1.])
W = np.array([0., 0., 0., 0., 0.3, 0.4, 0.8, 0.9, 1., 1., 1., 1.])

def test_bspline_getters():
    """
    Test getters of knot vector and degree for a 1D BSpline parametric space
    """
    b1 = BSpline(2, U)
    assert (b1.knot_vector == np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])).all()
    assert b1.degree == 2


def test_ndimension_bspline():
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
    ref_2dim = np.array([[2.32438017e-03, 5.68181818e-03, 2.58264463e-04],
                       [1.62706612e-01, 3.97727273e-01, 1.80785124e-02],
                       [1.16219008e-01, 2.84090909e-01, 1.29132231e-02]])

    funs = surf.basis_funs_nd(span, u)
    assert np.allclose(np.sum(funs), 1.0, rtol = 1.e-9)
    assert np.allclose(funs, ref_2dim, rtol = 1.e-9)


    u = np.array([0.3, 0.72, 0.45])
    span = vol.find_span_nd(u)
    assert( span == [2, 5, 5]).all()
    ref_3dim = np.array([[[1.77169421e-04, 4.13739669e-04, 6.95592287e-05, 6.88705234e-07],
                        [1.63881715e-03, 3.82709194e-03, 6.43422865e-04, 6.37052342e-06],
                        [3.98631198e-04, 9.30914256e-04, 1.56508264e-04, 1.54958678e-06]],
                       [[1.24018595e-02, 2.89617769e-02, 4.86914601e-03, 4.82093664e-05],
                        [1.14717200e-01, 2.67896436e-01, 4.50396006e-02, 4.45936639e-04],
                        [2.79041839e-02, 6.51639979e-02, 1.09555785e-02, 1.08471074e-04]],
                       [[8.85847107e-03, 2.06869835e-02, 3.47796143e-03, 3.44352617e-05],
                        [8.19408574e-02, 1.91354597e-01, 3.21711433e-02, 3.18526171e-04],
                        [1.99315599e-02, 4.65457128e-02, 7.82541322e-03, 7.74793388e-05]]])

    funs = vol.basis_funs_nd(span, u)
    assert np.allclose(np.sum(funs), 1.0, rtol=1.e-9)
    assert np.allclose(funs, ref_3dim, rtol=1.e-9)


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

        assert np.allclose(dfuns[1, :],
                           (val_plus - val_minus)/(2.*eps),
                           rtol = 1.e-6)


if __name__ == "__main__":
    test_bspline_getters()
    test_ndimension_bspline()
    test_functions_derivatives()
