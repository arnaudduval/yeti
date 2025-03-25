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

The shape of a Kirchhoff-Love shell catenary is optimized versus its maximal
bending moment
"""

import os
import numpy as np

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, \
      IGAmanip as manip
from yeti_iga.preprocessing.igaparametrization import OPTmodelling

PNORM = 40


def test_grad_catenary():
    """
    Test computed gradients values
    """
    # # Creation of the IGA object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/shellArch')

    nb_deg_design = np.array([1, 0, 0])
    nb_ref_design = np.array([2, 0, 0])

    iga_model.refine(nb_ref_design, nb_deg_design)

    nb_var = int(iga_model.n_kv[0, 0] - iga_model.j_pqr[0, 0] - 3)

    def altitude(coords0, igapara, var):
        igapara.coords[:, :] = coords0[:, :]
        count = 1
        for v in var:
            icps = manip.get_boundCPindice_wEdges(igapara.n_kv, igapara.j_pqr,
                                                  igapara.dim, 1,
                                                  num_patch=0, offset=count)
            igapara.coords[2, igapara.ind_cp_by_patch[0][icps]-1] += v*10.
            count += 1

    # Build the optimization pb
    nb_deg_an = np.array([2, 2, 0]) - nb_deg_design
    nb_ref_an = np.array([5, 1, 0]) - nb_ref_design

    optim_problem = OPTmodelling(iga_model, nb_var, altitude,
                                 nb_degreeElevationByDirection=nb_deg_an,
                                 nb_refinementByDirection=nb_ref_an)

    xk_test = np.array([0.31499875, 0.41392725, 0.41392725, 0.31499877])

    ref_volume = 7.08008736685128
    ref_stress = np.array([82.68791113, 0., 0., 13.25017662, 0., 0.])
    ref_gradvol = np.array([1.90690505, 2.35165181, 2.35165179, 1.90690522])
    ref_gradstress = np.array([74.79151097, -58.90890841,
                               -58.90889942, 74.79151713])

    assert np.allclose(optim_problem.compute_volume(xk_test),
                       ref_volume, rtol=1.e-5)
    assert np.allclose(
        optim_problem.compute_stressAggreg(xk_test, pnorm=PNORM),
        ref_stress, rtol=1.e-5)
    assert np.allclose(optim_problem.compute_gradVolume_AN(xk_test),
                       ref_gradvol, rtol=1.e-5)
    assert np.allclose(  # Some non sliced components have NaN values
        optim_problem.compute_gradStressAggreg_AN(xk_test, pnorm=PNORM)[:, 3],
        ref_gradstress, rtol=1.e-5)


if __name__ == '__main__':
    test_grad_catenary()
