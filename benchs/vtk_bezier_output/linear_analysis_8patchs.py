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
Generate VTU output using Bezier for a strong (common control points)
assembly of 2 patchs.
No test is made on the content of output files.
"""

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# modeleIGA = IGAparametrization(filename='BentPipeQuarter')
modeleIGA = IGAparametrization(filename='inputs/BentPipe_8patchs')

# Model refinement
nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)


nb_ref[:, 0] = np.array([3, 2, 3])
nb_ref[:, 1] = np.array([3, 2, 3])
nb_ref[:, 2] = np.array([3, 2, 3])
nb_ref[:, 3] = np.array([3, 2, 3])
nb_ref[:, 4] = np.array([3, 2, 3])
nb_ref[:, 5] = np.array([3, 2, 3])
nb_ref[:, 6] = np.array([3, 2, 3])
nb_ref[:, 7] = np.array([3, 2, 3])

nb_deg[:, 0] = np.array([0, 1, 0])
nb_deg[:, 1] = np.array([0, 1, 0])
nb_deg[:, 2] = np.array([0, 1, 0])
nb_deg[:, 3] = np.array([0, 1, 0])
nb_deg[:, 4] = np.array([0, 1, 0])
nb_deg[:, 5] = np.array([0, 1, 0])
nb_deg[:, 6] = np.array([0, 1, 0])
nb_deg[:, 7] = np.array([0, 1, 0])

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
K2solve = Ktot[idof, :][:, idof]
del Kside, data, row, col

# Resolution
x = sp.linalg.spsolve(K2solve, Fb[idof])

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))


# VTK output using Bezier elements
for i_patch in range(modeleIGA._nb_patch):
    pp.generate_vtu_bezier(**modeleIGA.get_inputs4postproc_bezier(
        i_patch+1,
        f'BentPipeBezier_8patchs_P{i_patch+1}',
        SOL.transpose(),
        ))

