# Copyright 2020 Thibaut Hirschler
# Copyright 2020 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
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
IsoGeometric Shape Optimization.
Arch Computat Methods Eng (2020). https://doi.org/10.1007/s11831-020-09458-6

The shape of Kirchhoff-Love shell roof is optimized versus its compliance
Volume is kept constant
Resuting shape is compared to reference numerical results

"""

import os
import numpy as np

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, \
      IGAmanip as manip
from yeti_iga.preprocessing.igaparametrization import OPTmodelling


def test_grad_square_shell_roof():
    """
    Test computed gradients values
    """

    # Creation of the IGA object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/squareShellPlate')

    iga_model.refine(nb_refinementByDirection=np.array([[2, 2, 0]]).T,
                     nb_degreeElevationByDirection=np.array([[1, 1, 0]]).T)

    # Shape Parametrization
    vertex = manip.get_vertexCPindice(iga_model.n_kv, iga_model.j_pqr,
                                      iga_model.dim)[:4]
    freecp = np.setxor1d(np.arange(0, iga_model.nb_cp), vertex)
    nb_var = freecp.size

    def altitude(coords0, igapara, var):
        igapara.coords[:, :] = coords0[:, :]
        igapara.coords[2, freecp] = coords0[2, freecp] + 5.*var[:]

    optim_problem = \
        OPTmodelling(iga_model, nb_var, altitude,
                     nb_degreeElevationByDirection=np.array([0, 0, 0]).T,
                     nb_refinementByDirection=np.array([3, 3, 3]).T)

    xk_test = np.array([0.11579141, 0.47655367, 0.47655367, 0.11579141,
                        0.11579141, 0.21496752, 0.35017255, 0.35017255,
                        0.21496752, 0.11579141, 0.47655367, 0.35017255,
                        0.52252416, 0.52252416, 0.35017255, 0.47655367,
                        0.47655367, 0.35017255, 0.52252416, 0.52252416,
                        0.35017255, 0.47655367, 0.11579141, 0.21496752,
                        0.35017255, 0.35017255, 0.21496752, 0.11579141,
                        0.11579141, 0.47655367, 0.47655367, 0.11579141])

    ref_volume = 106.95664994656539
    ref_compliance = 0.2658143469250488
    ref_gradvol = np.array([-1.4684551, 2.55633407, 2.55633407, -1.4684551,
                            -1.4684551, -1.1877157, -0.99198918, -0.99198918,
                            -1.1877157, -1.4684551, 2.55633407, -0.99198918,
                            4.14112959, 4.14112959, -0.99198918, 2.55633407,
                            2.55633407, -0.99198918, 4.14112959, 4.14112959,
                            -0.99198918, 2.55633407, -1.4684551, -1.1877157,
                            -0.99198918, -0.99198918, -1.1877157, -1.4684551,
                            -1.4684551, 2.55633407, 2.55633407, -1.4684551])

    ref_gradcomp = np.array([0.42266874, -0.00427181, -0.00427181, 0.42266874,
                             0.42266874, 0.21505752, -0.16105879, -0.16105879,
                             0.21505752, 0.42266874, -0.00427181, -0.16105879,
                             -0.19598612, -0.19598612, -0.16105879,
                             -0.00427181, -0.00427181, -0.16105879,
                             -0.19598612, -0.19598612, -0.16105879,
                             -0.00427181, 0.42266874, 0.21505752, -0.16105879,
                             -0.16105879, 0.21505752, 0.42266874, 0.42266874,
                             -0.00427181, -0.00427181, 0.42266874])

    assert np.allclose(optim_problem.compute_volume(xk_test),
                       ref_volume, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_compliance_discrete(xk_test),
                       ref_compliance, rtol=1.e-5)

    assert np.allclose(optim_problem.compute_gradVolume_AN(xk_test),
                       ref_gradvol, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_gradCompliance_AN(xk_test),
                       ref_gradcomp, rtol=1.e-5)


if __name__ == '__main__':
    test_grad_square_shell_roof()
