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
This case is described in the following publication :
Hirschler, T., Bouclier, R., Duval, A. et al.
A New Lighting on Analytical Discrete Sensitivities in the Context of
IsoGeometric Shape Optimization.
Arch Computat Methods Eng (2020). https://doi.org/10.1007/s11831-020-09458-6

The shape of a 2D solid plate with a hole is optimized versus its compliance
Volume is kept constant
Resutling hole shape is compared to a circle

"""

import os
import numpy as np

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.preprocessing.igaparametrization import OPTmodelling


def test_grad_plate_with_hole():
    """
    Test computed gradients values
    """

    # Creation of the IGA object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/plateWithHole')

    # Shape Parametrization
    nb_var = 6

    def holeshape(coords0, igapara, var):
        """
        Define plate hole shape by moving controls points as a function of
        design variables.
        """
        move_x = np.arange(0, 3, dtype=np.intp)
        move_y = np.arange(1, 4, dtype=np.intp)

        igapara.coords[:, :] = coords0[:, :]
        igapara.coords[0, move_x] += var[:move_x.size]
        igapara.coords[1, move_y] += var[move_x.size:]

    # Build the optimization pb
    nb_deg = np.array([1, 1, 0])
    nb_ref = np.array([3, 4, 0])

    optim_problem = OPTmodelling(iga_model, nb_var, holeshape,
                                 nb_degreeElevationByDirection=nb_deg,
                                 nb_refinementByDirection=nb_ref)

    xk_test = np.array([0.19723271, -0.07043156, -0.06654755, 0.06654755,
                        0.07043156, -0.19723271])

    ref_volume = 15.500081105173459
    ref_compliance = 0.00926406727182702
    ref_gradvol = np.array([0.21908539,  0.31977025,  0.21193217, -0.21193217,
                            -0.31977025, -0.21908539])
    ref_gradcomp = np.array([-0.0004263, -0.00062404, -0.00042611, 0.00042611,
                             0.00062404, 0.0004263])

    assert np.allclose(optim_problem.compute_volume(xk_test),
                       ref_volume, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_compliance_discrete(xk_test),
                       ref_compliance, rtol=1.e-5)

    assert np.allclose(optim_problem.compute_gradVolume_AN(xk_test),
                       ref_gradvol, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_gradCompliance_AN(xk_test),
                       ref_gradcomp, rtol=1.e-5)


if __name__ == '__main__':
    test_grad_plate_with_hole()
