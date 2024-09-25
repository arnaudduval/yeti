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
A single shell patch clamped at corner
Shell is divided in 9 elements (directly in input files) and only central
element is loaded with pressure.
Normal diosplacement at center is compared with reference result computed in
Abaqus.
"""

import sys

import numpy as np
import scipy.sparse as sp

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp

# Reference normal displacement at center computed with Abaqus
# (S8R elements, mesh size 0.1)
REF_ABA = -2.140E-4

# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename='9_elements_degree_3')

# Model refinement
nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)

nb_ref[:, 0] = np.array([2, 2, 0])
modeleIGA.refine(nb_ref, nb_deg)

# Matrix assembly
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1


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
# print(Fb[idof])

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'Shell_d3_9e',
    SOL.transpose(),
    nb_ref=np.array([4, 4, 0]),
    Flag=np.array([True, False, False])))

disp_center = pp.evaldisp(*modeleIGA.get_inputs4evaldisp(
    SOL.transpose(), np.array([0.5, 0.5, 0.]), numpatch=1))

print('yeti solution : ', disp_center[2])
print('reference solution : ', REF_ABA)
print('error : ', np.abs((REF_ABA-disp_center[2])/REF_ABA))

if np.abs((REF_ABA-disp_center[2])/REF_ABA) < 0.04:
    sys.exit(0)
else:
    sys.exit(-1)
