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

The shape of a solid 2D beam is optimized versus its maximal displacement
Volume is kept constant
"""

import os
import numpy as np

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, \
    IGAmanip as manip
from yeti_iga.preprocessing.igaparametrization import OPTmodelling


def test_grad_cantilever_beam():
    """
    Test computed gradients values
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/cantileverBeam')

    iga_model.refine(np.array([2, 0, 0]), np.array([1, 0, 0]))

    topcps = manip.get_directionCP(iga_model, 4, 0, 0)-1

    nb_var = topcps.size
    h_max = 10.0
    h_min = 1.5

    def beam_height(coords0, igapara, var):
        """
        Compute beam height depending on control points coordinates
        """
        igapara.coords[:, :] = coords0[:, :]
        igapara.coords[1, topcps] = (h_max - h_min)*var[:] + h_min

    # Build the optimization pb
    optim_problem = OPTmodelling(
        iga_model, nb_var, beam_height,
        nb_degreeElevationByDirection=np.array([0, 1, 0]),
        nb_refinementByDirection=np.array([3, 3, 0]))

    xi4disp = np.array([0.75, 0., 0.])

    xk_test = np.array([1.35290147, 1.48566331, 0.96101777,
                        0.07349478, 0.03275368, 0.27143277])

    ref_volume = 210.
    ref_disp = np.array([-0.00420967, -0.00946925, 0.])
    ref_gradvol = np.array([21.25, 42.5, 63.75, 63.75, 42.5, 21.25])
    ref_graddisp = np.array([[1.43443470e-04, 8.39615330e-04, 0.],
                             [1.80603238e-04, 7.93166153e-04, 0.],
                             [2.04967108e-03, 8.44682388e-03, 0.],
                             [1.38449992e-02, 3.92369606e-02, 0.],
                             [5.78761198e-03, 1.08681108e-02, 0.],
                             [7.36744690e-07, -1.02986420e-05, 0.]])

    assert np.allclose(optim_problem.compute_gradVolume_AN(xk_test),
                       ref_gradvol, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_volume(xk_test),
                       ref_volume, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_displacement(xk_test, xi4disp),
                       ref_disp, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_gradDisplacement_AN(xk_test, xi4disp),
        ref_graddisp, rtol=1.e-5)


if __name__ == '__main__':
    test_grad_cantilever_beam()
