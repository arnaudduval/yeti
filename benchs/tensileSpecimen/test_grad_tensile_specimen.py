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

The shape of a 2D solid tensile specimen is optimized versus its maximal Von
Mises stress
"""

import os
import numpy as np

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, \
    IGAmanip as manip
from yeti_iga.preprocessing.igaparametrization import OPTmodelling

PNORM = 40


def test_grad_tensile_specimen():
    """
    Test computed gradients values
    """

    # Creation of the IGA object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/tensileSpecimen')

    nb_deg_design = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref_design = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_deg_design[0, 1] = 2
    nb_ref_design[0, 1] = 2

    iga_model.refine(nb_ref_design, nb_deg_design)

    # --
    # Shape Parametrization

    top_cps_mid = manip.get_directionCP(iga_model, 4, 1, 0) - 1

    nb_var = top_cps_mid.size - 2

    def height_mid_patch(coords0, igapara, var):
        igapara.coords[:, :] = coords0[:, :]
        igapara.coords[1, top_cps_mid[1:-1]] += var[:]

    # Build the optimization problem
    nb_deg_an = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref_an = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

    nb_deg_an[:2, :] = 2
    nb_ref_an[:2, :] = 4
    nb_deg_an = np.maximum(nb_deg_an-nb_deg_design, 0)
    nb_ref_an = np.maximum(nb_ref_an-nb_ref_design, 0)

    optim_problem = OPTmodelling(iga_model, nb_var, height_mid_patch,
                                 nb_degreeElevationByDirection=nb_deg_an,
                                 nb_refinementByDirection=nb_ref_an)

    xk_test = np.array([-0.00891814, -0.04949665, -0.57123583,
                        -1.04720008, -0.34795149])

    ref_vonmises = 28.246541380870248
    ref_gradvonmises = np.array([0.02645387, 0.21844495, 2.22742017,
                                 -3.76106857, 1.3700102])

    assert np.allclose(
        optim_problem.compute_vonmisesAggreg(xk_test, pnorm=PNORM),
        ref_vonmises, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_gradVonMisesAggreg_AN(xk_test, pnorm=PNORM),
        ref_gradvonmises, rtol=1.e-5)


if __name__ == '__main__':
    test_grad_tensile_specimen()
