import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "build"))


import numpy as np
from bspline import BSpline, BSplineSurface, BSplineVolume, ControlPointManager, Patch
from bspline import IGABasis1D, PatchDOFManager, GlobalDOFManager, PatchIntegrator, MaterialProperties
import matplotlib
import matplotlib.pyplot as plt

U = np.array([0., 0., 0.,0.33, 0.66, 1., 1., 1.])
V = np.array([0., 0., 0., 0.2, 0.4, 0.6, 0.8, 1., 1., 1.])
W = np.array([0., 0., 0., 0., 0.3, 0.4, 0.8, 0.9, 1., 1., 1., 1.])















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


    # test_test_test()
    # exit()

    test_BSpline_getters()
    test_ND_BSpline()
    test_cp_manager()
    test_evaluation()
    test_evaluation_low_continuity()
    test_span_iterator()
    test_functions_derivatives()
    test_coupling_strong()
    test_integration_1elt_lin_rect()
    test_integration_1elt_lin_square_1()
    test_integration_1elt_d2_square_1()
    test_integration_1elt_d2_rect()
    test_integration_2_elements_C0()
    test_integration_2_elements_C1()
    print("All tests finshed !!!")

