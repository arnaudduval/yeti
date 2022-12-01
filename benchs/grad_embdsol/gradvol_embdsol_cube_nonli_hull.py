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
"""

import numpy as np

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import OPTmodelling

if __name__ == "__main__":

    modeleIGA = IGAparametrization(filename='embd_cube_nonli_hull')

    # Set arguments for model refinement
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    additional_knots = {"patches": np.array([]),
                        "1": np.array([]),
                        "2": np.array([]),
                        "3": np.array([])}

    initcoords = modeleIGA.coords

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
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)

    nb_var = modeleIGA.coords.size

    optPB = OPTmodelling(modeleIGA, nb_var, shapemodif,
                         nb_degreeElevationByDirection=nb_deg,
                         nb_refinementByDirection=nb_ref)

    x0 = np.zeros((nb_var))
    v0 = optPB.compute_volume(x0, listpatch=[0, 1])
    gradV_DF = optPB.compute_gradVolume_DF(x0)
    gradV_AN = optPB.compute_gradVolume_AN(x0)

    error = np.linalg.norm(gradV_DF - gradV_AN) / v0
    print(f"Error : {error:.02E}")

    assert error < 1.e-5
