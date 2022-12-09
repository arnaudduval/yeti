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

import sys
import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import OPTmodelling
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp


modeleIGA = IGAparametrization(filename='embd_cube_nonli')

# Set arguments for model refinement
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
additional_knots = {"patches": np.array([]),
                    "1": np.array([]),
                    "2": np.array([]),
                    "3": np.array([])}

nb_ref[:, 1] = [0, 0, 0]

modeleIGA.refine(nb_ref, nb_deg, additional_knots)
initcoords = modeleIGA._COORDS


# Build optim model
def shapemodif(coords0, igapara, var):
    igapara._COORDS[:, :] = coords0[:, :]
    dim = coords0.shape[0]
    n_cp = coords0.shape[1]
    i = 0
    # SHape change on all coords of all control points
    for i_cp in range(n_cp):
        for i_dim in range(dim):
            igapara._COORDS[i_dim, i_cp] = initcoords[i_dim, i_cp] + var[i]
            i += 1

    return None


# Refinement from optim model to analysis model
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

nb_var = modeleIGA._COORDS.size

optPB = OPTmodelling(modeleIGA, nb_var, shapemodif,
                     nb_degreeElevationByDirection=nb_deg,
                     nb_refinementByDirection=nb_ref)

x0 = np.zeros((nb_var))
v0 = optPB.compute_volume(x0, listpatch=[0, 1])
gradV_DF = optPB.compute_gradVolume_DF(x0)
gradV_AN = optPB.compute_gradVolume_AN(x0)
# dim = modeleIGA._COORDS.shape[0]
# n_cp = modeleIGA._COORDS.shape[1]
# i = 0
# for i_cp in range(n_cp):
#     for i_dim in range(dim):
#         print(i_cp, i_dim, gradV_DF[i], gradV_AN[i])
#         i += 1

error = np.linalg.norm(gradV_DF - gradV_AN) / v0


assert(error < 1.e-5)

exit(0)

print(optPB.compute_volume(x0, listpatch=[0, 1]))

print(optPB.compute_compliance_discrete(x0))

print(optPB.compute_gradCompliance_AN(x0))

print(optPB.compute_gradVolume_AN(x0, listpatch=[0, 1]))


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
