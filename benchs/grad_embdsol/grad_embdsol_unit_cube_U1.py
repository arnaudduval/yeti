# Copyright 2022 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the terms 
# of the GNU Lesser General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
# PURPOSE. See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along 
# with Yeti. If not, see <https://www.gnu.org/licenses/>

# -*- coding: utf-8 -*-

import sys
import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import OPTmodelling
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp


modeleIGA = IGAparametrization(filename='unit_cube_U1')

# Set arguments for model refinement
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
additional_knots = {"patches": np.array([]),
                    "1": np.array([]),
                    "2": np.array([]),
                    "3": np.array([])}

nb_ref[:, 0] = [0, 0, 0]

modeleIGA.refine(nb_ref, nb_deg, additional_knots)
initcoords = modeleIGA._COORDS

# Build optim model
def shapemodif(coords0, igapara, var):
    igapara._COORDS[:,:] = coords0[:,:]
    # shape change is made on nodes 2 and 4
    igapara._COORDS[0, 1] = initcoords[0, 1] + var[0]
    igapara._COORDS[1, 1] = initcoords[1, 1] + var[1]
    igapara._COORDS[2, 1] = initcoords[2, 1] + var[2]
    igapara._COORDS[0, 3] = initcoords[0, 3] + var[3]
    igapara._COORDS[1, 3] = initcoords[1, 3] + var[4]
    igapara._COORDS[2, 3] = initcoords[2, 3] + var[5]

    return None

# Refinement from optim model to analysis model
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

nb_var = 6

optPB = OPTmodelling(modeleIGA, nb_var, shapemodif,
                     nb_degreeElevationByDirection=nb_deg,
                     nb_refinementByDirection=nb_ref)

x0 = np.zeros((6))
print(optPB.compute_volume(x0))                     
print(optPB.compute_compliance_discrete(x0))

print("FD : ", optPB.compute_gradCompliance_FD(x0))
print("AN : ", optPB.compute_gradCompliance_AN(x0))

#print(optPB.compute_gradVolume_AN(x0))

exit()

# Static study
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof] - 1
print(idof)
print(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_free)
data, row, col, Fb = build_stiffmatrix(
                        *modeleIGA.get_inputs4system_elemStorage())
Kside = sp.coo_matrix((data, (row, col)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside + Kside.transpose()
K2solve = Ktot[idof, :][:, idof]
del Kside, data, row, col

x = sp.linalg.spsolve(K2solve, Fb[idof])

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'test_embd_cube',
    SOL.transpose(), nb_ref=[3, 3, 3], Flag=[True, False, False]))
modeleIGA.generate_vtk4controlMeshVisu(
    'test_embd_cube_map_cp', 0)
modeleIGA.generate_vtk4controlMeshVisu(
    'test_embd_cube_embded_cp', 1)