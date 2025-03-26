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
# FOR A PARTICULAR
# PURPOSE. See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This case is described in the following publication :
Hirschler, T., Bouclier, R., Duval, A. et al.
A New Lighting on Analytical Discrete Sensitivities in the Context of
IsoGeometric Shape Optimization.
Arch Computat Methods Eng (2020). https://doi.org/10.1007/s11831-020-09458-6

The shape of a solid 3D elephant trunk is optimized versus its eigenfrequencies
Volume is kept constant
"""

import os
import numpy as np

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, \
                                                      IGAmanip as manip
from yeti_iga.preprocessing.igaparametrization import OPTmodelling


def test_grad_elephant_trunk():
    """
    Test computed gradients
    """
    # Creation of the IGA object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/elephantTrunk')

    nb_deg_design = np.array([0, 0, 1])
    nb_ref_design = np.array([0, 0, 2])

    iga_model.refine(nb_ref_design, nb_deg_design)
    nb_var = int((iga_model.n_kv[2, 0] - iga_model.j_pqr[2, 0] - 1)*2)

    def dilatation(coords0, igapara, var):
        igapara.coords[:, :] = coords0[:, :]
        for i in np.arange(0, var.size/2, dtype=np.intp):
            icps = np.intp(manip.get_directionCP(igapara, 5, 0, i) - 1)
            igapara.coords[0, icps] *= var[i]
            igapara.coords[1, icps] *= var[int(i+var.size/2)]

    nb_deg_analysis = np.maximum(np.array([0, 0, 1]) - nb_deg_design, 0)
    nb_ref_analysis = np.maximum(np.array([2, 2, 4]) - nb_ref_design, 0)

    optim_problem = OPTmodelling(iga_model, nb_var, dilatation,
                                 nb_degreeElevationByDirection=nb_deg_analysis,
                                 nb_refinementByDirection=nb_ref_analysis)

    nb_freq = 5
    xk_test = np.array([1.62707153, 1.54842318, 1.2854958, 0.77428138,
                        0.25, 0.25, 1.80428529, 1.85106827, 1.35227388,
                        0.62181301, 0.25, 0.25])
    ref_volume = 39.27003
    ref_grad_vol = np.array([4.70730, 8.75269, 10.18088, 5.55013, 1.98951,
                             0.703168, 4.1428859, 7.69870, 9.63348, 6.07530,
                             2.23462, 0.72312])
    ref_eigen = np.array([53.94943, 57.40211, 170.11873, 186.88762, 622.07787])
    ref_grad_eigen = np.array(
        [[19.67357, 5.17106, 17.16848, 42.10873, 45.12770],
         [27.29204, 7.29723, 11.52960, 22.96115, -4.84958],
         [23.97798, 6.24620, -23.81477, -3.10598, -4.52177],
         [3.55793, -7.50407, 27.63472, 247.38143595, -19.55357],
         [-20.87234, -38.03131, 17.62084, 242.09844, 36.04903],
         [-32.5637, -51.4918, -169.18556, -187.67891, -255.1387],
         [5.93584, 13.87351, 45.06018, 13.30093, 92.47524],
         [7.62753, 18.37585, 28.37321, 5.64129, 5.94001],
         [2.91972, 27.01660, -11.94378, -20.10759, 185.69639],
         [-18.79006, 33.89090, 207.90137, 67.95655, 353.49149],
         [-35.45262, -12.62862, 240.17743, 20.99953, 707.75118],
         [-34.59689, -49.42659, -139.29835, -225.96793, -75.76490]])

    assert np.allclose(optim_problem.compute_volume(xk_test),
                       ref_volume, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_gradVolume_AN(xk_test),
                       ref_grad_vol, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_vibrationMode(xk_test, nb_freq)[0],
        ref_eigen, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_gradVibration_AN(xk_test, nb_frq=nb_freq),
        ref_grad_eigen, rtol=1.e-5)


if __name__ == '__main__':
    test_grad_elephant_trunk()
