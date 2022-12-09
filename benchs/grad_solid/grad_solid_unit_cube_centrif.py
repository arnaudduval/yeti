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
Test on a unit cube made with solid element U1 of degree 1
Structure is loaded with a centrifugal body force
Gradients of compliance and volume are computed with 3 methods:
 - finite differences
 - semi analytical
 - analytical
Results must be equal within a tolerance
"""

import numpy as np

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import OPTmodelling

if __name__ == "__main__":

    modeleIGA = IGAparametrization(filename='unit_cube_U1')

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
    v0 = optPB.compute_volume(x0)
    c0 = optPB.compute_compliance_discrete(x0)
    print(f"Volume : {v0:2.5E}")
    print(f"Compliance : {c0:2.5E}")

    gradV_AN = optPB.compute_gradVolume_AN(x0)
    gradV_FD = optPB.compute_gradVolume_DF(x0)

    gradC_FD = optPB.compute_gradCompliance_FD(x0)
    gradC_AN = optPB.compute_gradCompliance_AN(x0)
    gradC_SAN = optPB.compute_gradCompliance_semiAN(x0)

    print("gradV FD : ", gradV_FD)
    print("gradV AN : ", gradV_AN)

    print("gradC FD : ", gradC_FD)
    print("gradC AN : ", gradC_AN)
    print("gradC SAN : ", gradC_SAN)

    errorV = np.linalg.norm(gradV_AN - gradV_FD)/v0

    errorC_AN_FD = np.linalg.norm(gradC_AN - gradC_FD)/c0
    errorC_SAN_FD = np.linalg.norm(gradC_SAN - gradC_FD)/c0
    errorC_AN_SAN = np.linalg.norm(gradC_AN - gradC_SAN)/c0

    print(f"Error on grad V : {errorV:2.5E}")
    print(f"Error on grad C (AN/FD) : {errorC_AN_FD:2.5E}")
    print(f"Error on grad C (SAN/FD) : {errorC_SAN_FD:2.5E}")
    print(f"Error on grad C (AN/SAN) : {errorC_AN_SAN:2.5E}")

    assert errorV < 1.e6
    assert errorC_AN_FD < 1.e5
    assert errorC_SAN_FD < 1.e5
    assert errorC_AN_SAN < 1.e5
