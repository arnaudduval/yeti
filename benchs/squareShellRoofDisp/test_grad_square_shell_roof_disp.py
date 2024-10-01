# Copyright 2020 Thibaut Hirschler
# Copyright 2020 Arnaud Duval

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
This cas is described in the following publication :
Hirschler, T., Bouclier, R., Duval, A. et al.
A New Lighting on Analytical Discrete Sensitivities in the Context of
IsoGeometric Shape Optimization. Arch Computat Methods Eng (2020).
https://doi.org/10.1007/s11831-020-09458-6

The shape of a Kirchhoff-Love shell roof is optimized versus its maximal
displacement.
Volume is kept constant
"""

import os
import numpy as np

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, \
                                             IGAmanip as manip
from yeti_iga.preprocessing.igaparametrization import OPTmodelling

P_NORM = 20


def test_grad_square_shell_roof_disp():
    """
    Test computed gradients values
    """

    # Creation of the IGA object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/squareShellPlate')

    iga_model.refine(np.array([2, 2, 0]), np.array([1, 1, 0]))

    # Shape Parametrization
    vertex = manip.get_vertexCPindice(iga_model.n_kv, iga_model.j_pqr,
                                      iga_model.dim)[:4]
    freecp = np.setxor1d(np.arange(0, iga_model.nb_cp), vertex)
    nb_var = freecp.size

    def altitude(coords0, igapara, var):
        """
        Change z coordinate of control points as a function of design variables
        """
        igapara.coords[:, :] = coords0[:, :]
        igapara.coords[2, freecp] = coords0[2, freecp] + 5.*var[:]

    # Build the optimization pb
    nb_deg = np.maximum(np.array([0, 0, 0]), 0)
    nb_ref = np.maximum(np.array([3, 3, 0]), 0)

    optim_problem = OPTmodelling(iga_model, nb_var, altitude,
                                 nb_degreeElevationByDirection=nb_deg,
                                 nb_refinementByDirection=nb_ref)

    xk_test = np.array([0.04693029, 0.20588302, 0.20588302, 0.04693029,
                        0.04693029, 0.04693029, 0.47372698, 0.47372698,
                        0.04693029, 0.04693029, 0.20588302, 0.47372698,
                        0.57762737, 0.57762737, 0.47372698, 0.20588302,
                        0.20588302, 0.47372698, 0.57762737, 0.57762737,
                        0.47372698, 0.20588302, 0.04693029, 0.04693029,
                        0.47372698, 0.47372698, 0.04693029, 0.04693029,
                        0.04693029, 0.20588302, 0.20588302, 0.04693029])

    ref_volume = 116.19997424272571
    ref_disp_aggreg = 0.00024629470705899787
    ref_gradvol = np.array([-3.23237683, -5.29094997, -5.29094997, -3.23237683,
                            -3.23237683, 0.19737303, 5.51194684, 5.51194684,
                            0.19737303, -3.23237683, -5.29094997, 5.51194684,
                            7.93954514, 7.93954514, 5.51194684, -5.29094997,
                            -5.29094997, 5.51194684, 7.93954514, 7.93954514,
                            5.51194684, -5.29094997, -3.23237683, 0.19737303,
                            5.51194684, 5.51194684, 0.19737303, -3.23237683,
                            -3.23237683, -5.29094997, -5.29094997, -3.23237683
                            ])
    ref_graddisp_aggreg = \
        np.array([-5.63355389e-04, 5.89090961e-05, 5.89090961e-05,
                  -5.63355389e-04, -5.63355389e-04, -9.40645528e-05,
                  -1.02028135e-05, -1.02028135e-05, -9.40645528e-05,
                  -5.63355389e-04,  5.89090961e-05, -1.02028135e-05,
                  2.38128879e-05,  2.38128879e-05, -1.02028135e-05,
                  5.89090961e-05, 5.89090961e-05, -1.02028135e-05,
                  2.38128879e-05,  2.38128879e-05, -1.02028135e-05,
                  5.89090961e-05, -5.63355389e-04, -9.40645528e-05,
                  -1.02028135e-05, -1.02028135e-05, -9.40645528e-05,
                  -5.63355389e-04, -5.63355389e-04,  5.89090961e-05,
                  5.89090961e-05, -5.63355389e-04])

    assert np.allclose(optim_problem.compute_volume(xk_test),
                       ref_volume, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_displacementAggreg(xk_test, pnorm=P_NORM),
        ref_disp_aggreg, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_gradVolume_AN(xk_test),
                       ref_gradvol, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_gradDisplacementAggreg_AN(xk_test, pnorm=P_NORM),
        ref_graddisp_aggreg, rtol=1.e-5)


if __name__ == '__main__':
    test_grad_square_shell_roof_disp()
