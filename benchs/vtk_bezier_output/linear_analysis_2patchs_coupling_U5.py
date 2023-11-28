# Copyright 2023 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>


# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generate VTU output using Bezier en assembly of 2 patchs with Mortar
coupling
No test is made on the content of output files.
"""

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# modeleIGA = IGAparametrization(filename='BentPipeQuarter')
modeleIGA = IGAparametrization(filename='inputs/BentPipe_2patchs_coupling_U5')

# Model refinement
nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)


nb_ref[:, 0] = np.array([1, 1, 2])
nb_ref[:, 1] = np.array([1, 1, 2])

nb_deg[:, 0] = np.array([1, 2, 1])  # Degree 3 in all directions
nb_deg[:, 1] = np.array([1, 2, 1])  # Degree 3 in all directions

# Coupling patches
nb_ref[:, 2] = np.array([1, 2, 0])
nb_ref[:, 3] = np.array([1, 2, 0])

nb_deg[:, 2] = np.array([2, 2, 0])  # Degree 2
nb_deg[:, 3] = np.array([2, 2, 0])  # Degree 2

print(modeleIGA.nb_patch)
modeleIGA.refine(nb_ref, nb_deg)

# Matrix assembly
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

print("Build stiffmatrix")

data, row, col, Fb = build_stiffmatrix(
                        *modeleIGA.get_inputs4system_elemStorage())

Kside = sp.coo_matrix((data, (row, col)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside + Kside.transpose()
del Kside, data, row, col


print("Build coupling matrix")
Cdata, Crow, Ccol = cplg_matrixU5( *modeleIGA.get_inputs4cplgmatrixU5() )
Cside = sp.coo_matrix((Cdata,(Crow,Ccol)), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ctot  = Cside + Cside.transpose()
del Cdata,Crow,Ccol,Cside


# Resolution
K2solve = Ktot[idof,:][:,idof]
C2solve = Ctot[idof,:][:,idof] * K2solve.max()
x = sp.linalg.spsolve(K2solve + C2solve, Fb[idof])

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'couplingU5',SOL.transpose(),nb_ref=np.array([3,3,3]),
    Flag=np.array([True,True,False])))

exit()

# VTK output using Bezier elements
for i_patch in range(modeleIGA._nb_patch):
    pp.generate_vtu_bezier(**modeleIGA.get_inputs4postproc_bezier(
        i_patch+1,
        f'BentPipeBezier_2patchs_P{i_patch+1}',
        SOL.transpose(),
        ))

