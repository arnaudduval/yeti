# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test point inversion on a trivariate model
"""

import os
import numpy as np


from yeti_iga import Patch, QuarterPipe
from yeti_iga import IgaModel, ElasticMaterial
from yeti_iga.exchange import load_geomdl_json
from yeti_iga.utils import trivariate_point_inversion
from yeti_iga.preprocessing.igaparametrization import IGAparametrization



def test_trivariate_pt_inversion_linear(tmp_path):
    """
    Define a trivariate patch and perform point inversion
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/beam-dist-field')

    # Refine modele
    iga_model.refine(nb_refinementByDirection=np.array([2, 2, 5]),
                     nb_degreeElevationByDirection=np.array([1, 1, 1]))

    params = iga_model.get_inputs4point_inversion(np.array([10., 2., 0.7]),
                                                  eps1=1.e-5)
    param, status, distance = trivariate_point_inversion(**params)

    assert np.allclose(param, np.array([2./3., 0.7, 0.5]), rtol=1.e-10)
    assert status == 0
    assert distance < 1.e-5


def test_trivariate_pt_inversion_nonlinear(tmp_path):
    """
    Define trivariate patches and perform point inversion
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    filename = 'svol_face.json'

    model = IgaModel("3D solid")
    patchs = load_geomdl_json(f'{script_dir}/{filename}')
    print(len(patchs))
    for patch in patchs:
        model.add_patch(patch)
    qpipe = QuarterPipe(10., 100., 120., ElasticMaterial(210000., 0.3))
    qpipe.translate([0., 0., -50.])
    model.add_patch(qpipe)

    point = np.array([-5.0, -38.0, 0.0])
    params = model.iga_param.get_inputs4point_inversion(point, i_patch=0, eps1=1.e-5)
    print(trivariate_point_inversion(**params))

    model.write_solution_vtu(np.zeros(shape=(model.iga_param.nb_dof_free)),
                              'test.vtu')



    return



    patch = Patch(element_type='U1',
                  degrees=np.array([1, 1, 1]),
                  knot_vectors=[np.array([0., 0., 1., 1.]),
                                np.array([0., 0., 1., 1.]),
                                np.array([0., 0., 1., 1.])],
                  control_points=np.array([[0., 0., 0.],
                                           [0., 0.1, 0.],
                                           [0., 0., 0.03],
                                           [0., 0.1, 0.03],
                                           [10., 0., 0.],
                                           [10., 0.1, 0.],
                                           [10., 0., 0.03],
                                           [10., 0.1, 0.03]]),
                  weights=np.array([1., 1., 1., 1., 1., 1., 1., 1.]),
                  connectivity=np.array([[7, 6, 5, 4, 3, 2, 1, 0]]),
                  spans=np.array([[1, 1, 1]]),
                  material=ElasticMaterial(210000., 0.3)
                  )
    model.add_patch(patch)

    point = np.array([1., 2., 3.])

    print(len(model.iga_param._IEN))
    for ie in model.iga_param._IEN:
        print(len(ie))
    print(len(model.iga_param._IEN_flat))
    print(f'{model.iga_param._COORDS.shape = }')

    trivariate_point_inversion(**model.iga_param.get_inputs4point_inversion(point, i_patch=0))

    return

    # model.write_solution_vtu(np.zeros(shape=(model.iga_param.nb_dof_free)),
    #                          'test.vtu')

    return

    patch = patchs[0]
    patch.invert_point(np.array([5., 0.05, 0.01]))
    return

    print(patch.knot_vectors[0])
    print(len(patch.knot_vectors[0]))
    # print(patch.spans[0, :])

    from yeti_iga import tools

    print(tools.find_span(0.5, patch.knot_vectors[0]))
    print(tools.find_span(0.49, patch.knot_vectors[0]))
    print(tools.find_span(0.51, patch.knot_vectors[0]))



    return

    patch = Patch(element_type='U1',
                  degrees=np.array([1, 1, 1]),
                  knot_vectors=[np.array([0., 0., 1., 1.]),
                                np.array([0., 0., 1., 1.]),
                                np.array([0., 0., 1., 1.])],
                  control_points=np.array([[0., 0., 0.],
                                           [0., 0.1, 0.],
                                           [0., 0., 0.03],
                                           [0., 0.1, 0.03],
                                           [10., 0., 0.],
                                           [10., 0.1, 0.],
                                           [10., 0., 0.03],
                                           [10., 0.1, 0.03]]),
                  weights=np.array([1., 1., 1., 1., 1., 1., 1., 1.]),
                  connectivity=np.array([[7, 6, 5, 4, 3, 2, 1, 0]]),
                  spans=np.array([[1, 1, 1]]),
                  material=None
                  )

    print(type(patch.knot_vectors))
    patch.invert_point(np.array([5., 0.05, 0.01]))


if __name__ == '__main__':
    # test_trivariate_pt_inversion_linear('.')
    test_trivariate_pt_inversion_nonlinear('.')
