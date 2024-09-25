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
A test to validate compliance gradient computation for embedded solid element.
A degree 1 cube is embedded in a degree 1 cube
Loading is a constant pressure applied on face 6
symmetry boundary condition are applied on faces 1,3 and 5
Analytical compliance gradient is compared with the one computed with finite
differences
Design variables drive the control points of embedded entity except the loaded
face.
"""

import numpy as np

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.preprocessing.igaparametrization import OPTmodelling

if __name__ == "__main__":

    modeleIGA = IGAparametrization(filename='embd_cube_lin_press')

    # Set arguments for model refinement
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    additional_knots = {"patches": np.array([]),
                        "1": np.array([]),
                        "2": np.array([]),
                        "3": np.array([])}

    modeleIGA.refine(nb_ref, nb_deg, additional_knots)

    initcoords = modeleIGA.coords

    # Build optim model
    def shapemodif(coords0, igapara, var):
        """
        A shape change that apply on all coordinates of all control points
        """
        igapara.coords[:, :] = coords0[:, :]
        dim = coords0.shape[0]
        i = 0

        for i_cp in range(8, 12):
            for i_dim in range(dim):
                igapara.coords[i_dim, i_cp] = initcoords[i_dim, i_cp] + var[i]
                i += 1

    # Refinement from optim model to analysis model
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)

    NB_VAR = 4 * 3      # 4 control points x 3 directions
    optPB = OPTmodelling(modeleIGA, NB_VAR, shapemodif,
                         nb_degreeElevationByDirection=nb_deg,
                         nb_refinementByDirection=nb_ref)

    x0 = np.zeros((NB_VAR))

    c0 = optPB.compute_compliance_discrete(x0)
    gradC_DF = optPB.compute_gradCompliance_FD(x0)
    gradC_AN = optPB.compute_gradCompliance_AN(x0)

    error = np.linalg.norm(gradC_DF - gradC_AN) / c0
    print(f"Compliance : {c0:.02E}")
    print(f"Error : {error:.02E}")

    assert error < 1.e-6
