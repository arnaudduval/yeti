"""
Cross-check WeightedQuadrature1D (C++ port) against pymfiga's Python
reference implementation -- both "method 1" and "method 2" (method "2" is
the only variant exercised by pymfiga's own validated regression
benchmarks, see benchs/pymfiga/, but both are ported here).
"""

import numpy as np
import pytest

from yeti_iga.future.bspline import BSpline, WeightedQuadrature1D
from yeti_iga.pymfiga.common.numerics.quadrature_rules.weighted_quadrature import (
    WeightedQuadrature,
)


@pytest.mark.parametrize(
    "degree,knot_vector",
    [
        (2, [0, 0, 0, 1, 1, 1]),                                  # single element
        (2, [0, 0, 0, 0.25, 0.5, 0.5, 0.75, 1, 1, 1]),             # multi-element, reduced C1
        (3, [0, 0, 0, 0, 1, 1, 1, 1]),                             # single element
        (3, [0, 0, 0, 0, 0.2, 0.4, 0.6, 0.8, 1, 1, 1, 1]),         # 5 elements, C2 interior
        (3, [0, 0, 0, 0, 0.3, 0.3, 0.3, 0.7, 1, 1, 1, 1]),         # interior knot multiplicity 3
    ],
    ids=["deg2-1elt", "deg2-multi-C1", "deg3-1elt", "deg3-multi-C2", "deg3-mult3"],
)
@pytest.mark.parametrize("quadtype", ["1", "2"])
def test_weighted_quadrature_matches_pymfiga(degree, knot_vector, quadtype):
    kv = np.array(knot_vector, dtype=float)

    bsp = BSpline(degree, kv)
    wq_cpp = WeightedQuadrature1D.build(bsp, quadtype)

    wq_py = WeightedQuadrature(degree, kv, quadtype=quadtype)
    wq_py.export_quadrature_rules()

    assert wq_cpp.quadtype == quadtype
    assert np.allclose(np.array(wq_cpp.quadpts), wq_py.quadpts, atol=1e-10)

    B0_cpp = np.array(wq_cpp.B0.todense())
    B1_cpp = np.array(wq_cpp.B1.todense())
    assert np.allclose(B0_cpp, wq_py.basis[0].toarray(), atol=1e-8)
    assert np.allclose(B1_cpp, wq_py.basis[1].toarray(), atol=1e-8)

    W00_cpp = np.array(wq_cpp.W00.todense())
    W01_cpp = np.array(wq_cpp.W01.todense())
    W10_cpp = np.array(wq_cpp.W10.todense())
    W11_cpp = np.array(wq_cpp.W11.todense())

    if quadtype == "2":
        # Method "2": W00==W01 and W10==W11, both on the C++ side and the
        # pymfiga reference -- checked directly rather than assumed.
        assert np.allclose(W00_cpp, W01_cpp)
        assert np.allclose(W10_cpp, W11_cpp)
    else:
        # Method "1": all four weight matrices are genuinely distinct.
        assert not np.allclose(W00_cpp, W01_cpp)
        assert not np.allclose(W10_cpp, W11_cpp)

    assert np.allclose(W00_cpp, wq_py.weights[0].toarray(), atol=1e-6)
    assert np.allclose(W01_cpp, wq_py.weights[1].toarray(), atol=1e-6)
    assert np.allclose(W10_cpp, wq_py.weights[2].toarray(), atol=1e-6)
    assert np.allclose(W11_cpp, wq_py.weights[3].toarray(), atol=1e-6)


def test_weighted_quadrature_fewer_points_than_gauss():
    """
    The headline benefit of weighted quadrature: fewer points than standard
    Gauss for the same spline space, once there are enough elements for the
    per-span vs. per-knot-break scaling to matter.
    """
    degree = 3
    kv = np.array([0, 0, 0, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1, 1, 1],
                  dtype=float)
    bsp = BSpline(degree, kv)
    wq = WeightedQuadrature1D.build(bsp)

    nbelem = len(np.unique(kv)) - 1
    n_gauss = nbelem * (degree + 1)

    assert len(wq.quadpts) < n_gauss
