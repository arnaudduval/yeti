# Copyright 2022 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Verify gradient coputation of a solid structure subjected to centrifugal body
force.
Results are compared with reference numerical results.
"""

import os

import numpy as np

from yeti_iga.preprocessing.igaparametrization import IGAparametrization, \
                                                      OPTmodelling


def test_grad_centrif_coupling_u1_c0():
    """
    Test computed gradients
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/centrif_U1_C0')

    iga_model.refine(nb_refinementByDirection=np.array([0, 2, 0]),
                     nb_degreeElevationByDirection=np.array([0, 1, 0]),
                     additional_knots={"patches": np.array([0]),
                                       "1": np.array([]),
                                       "2": np.array([0.2, 0.2]),
                                       "3": np.array([])})

    icp = np.where(iga_model.coords[1, :] > 1.2)[0]
    nb_var = np.unique(iga_model.coords[1, icp]).size
    print(nb_var)

    # Define shape change from design variables
    # min and max dimensions, assuming that design variables are in [0,1]
    min_dim = 0.1
    max_dim = 3.5

    def shapemodif(coords0, igapara, var):
        """"
        Shape parametrization
        Design variables drive coordinates of control point located at y > 1.2
        """
        igapara.coords[:, :] = coords0[:, :]
        # shape change is made on points with y coord higher than 3
        i = 0
        for y_loc in np.unique(igapara.coords[1, icp]):
            # WARNING exact real value comparison is unsafe
            jcp = np.where(igapara.coords[1, :] == y_loc)[0]

            igapara.coords[2, jcp[0]] = -(min_dim +
                                          var[i]*(max_dim-min_dim))/2.
            igapara.coords[2, jcp[1]] = -(min_dim +
                                          var[i]*(max_dim-min_dim))/2.
            igapara.coords[2, jcp[2]] = (min_dim +
                                         var[i]*(max_dim-min_dim))/2.
            igapara.coords[2, jcp[3]] = (min_dim +
                                         var[i]*(max_dim-min_dim))/2.

            i += 1

    nb_deg_analysis = np.array([1, 0, 1])
    nb_ref_analysis = np.array([2, 2, 2])

    optim_problem = OPTmodelling(iga_model, nb_var, shapemodif,
                                 nb_degreeElevationByDirection=nb_deg_analysis,
                                 nb_refinementByDirection=nb_ref_analysis)

    xk_test = np.array([0.93079, 0.69162, 0.13132, 0.01605, 0.16685])
    ref_compliance = 1.431973e-05
    ref_volume = 13.5
    ref_grad_comp = np.array([8.59656e-06, 1.53569e-05, 1.88451e-05,
                              1.51298e-05, 8.63750e-06])
    ref_grad_vol = np.array([4.08, 6.12, 6.12, 4.08, 2.04])

    assert np.allclose(optim_problem.compute_volume(xk_test),
                       ref_volume, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_gradVolume_AN(xk_test),
                       ref_grad_vol, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_compliance_discrete(xk_test),
        ref_compliance, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_gradCompliance_AN(xk_test),
        ref_grad_comp, rtol=1.e-5)


if __name__ == "__main__":
    test_grad_centrif_coupling_u1_c0()
