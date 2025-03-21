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

# -*- coding: utf-8 -*-

"""
A test to validate volume gradient computation for embedded solid element.
Structure is an embedded cube with non linear parametrization for embedded
patch.
Gradient of volume is computed with 2 methods:
 - finite differences
 - analytical
"""

import os

import numpy as np

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.preprocessing.igaparametrization import OPTmodelling


def test_gradvol_embdsol_cube_nonli_embd():
    """
    Compute gradients and compare values
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/embd_cube_nonli_embd_press')

    initcoords = iga_model.coords

    # Build optim model
    def shapemodif(coords0, igapara, var):
        """
        A shape change that apply on all coordinates of all control points
        """
        igapara.coords[:, :] = coords0[:, :]
        dim = coords0.shape[0]
        n_cp = coords0.shape[1]
        i = 0

        for i_cp in range(n_cp):
            for i_dim in range(dim):
                igapara.coords[i_dim, i_cp] = initcoords[i_dim, i_cp] + var[i]
                i += 1

    # Refinement from optim model to analysis model
    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

    nb_var = iga_model.coords.size

    optim_problem = OPTmodelling(iga_model, nb_var, shapemodif,
                                 nb_degreeElevationByDirection=nb_deg,
                                 nb_refinementByDirection=nb_ref)

    x0 = np.zeros((nb_var))
    v0 = optim_problem.compute_volume(x0, listpatch=[0, 1])

    error = np.linalg.norm(optim_problem.compute_gradVolume_DF(x0) -
                           optim_problem.compute_gradVolume_AN(x0)) / v0

    assert error < 1.e-6


if __name__ == "__main__":
    test_gradvol_embdsol_cube_nonli_embd()
