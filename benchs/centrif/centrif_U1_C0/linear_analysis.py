# Copyright 2020 Arnaud Duval

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

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Compute deflexion due to distributed pressure field on a beam
Compare with analytic result
"""

import sys
import time

import numpy as np
import scipy.sparse as sp

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.preprocessing.igaparametrization import IGAmanip as manip
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp

# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename='centrif_U1_C0')

# Refine modele
nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

nb_ref[:,0] = np.array([2,2,2])
nb_deg[:,0] = np.array([1,1,1])
additional_knots = {"patches": np.array([0]),
                    "1": np.array([]), "2": np.array([]), "3": np.array([0.6,0.6])}


# Initial refinement (none)
modeleIGA.refine(nb_ref, nb_deg,additional_knots)

# Matrix assembly
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

t1 = time.time()

data, row, col, Fb = build_stiffmatrix(
                        *modeleIGA.get_inputs4system_elemStorage())

Kside = sp.coo_matrix((data, (row, col)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside + Kside.transpose()
K2solve = Ktot[idof, :][:, idof]
del Kside, data, row, col
print(('\n Time to build stiffness matrix : %.2f s' % (time.time() - t1)))


# RESOLUTION : direct solver
t2 = time.time()
x  = sp.linalg.spsolve(K2solve, Fb[idof])
# LU = sp.linalg.splu(K2solve)
# x = LU.solve(Fb[idof])
print(('\n Time for direct solving :    %.2f s\n\n' % (time.time() - t2)))

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'centrif_U1_C0',
    SOL.transpose(),
    nb_ref=np.array([2, 2, 2]),
    Flag=np.array([True, False, False])))

# Get solution at interpolating control points and compare with Abaqus reference result
# Tolerance = max 2% on quais radial component
pt_coords = np.array([12.0, 1.5, 0.0])
ref_aba = np.array([1.60785E-6, -1.94110E-8, 3.15129E-7])

found = False

for idx in range(np.shape(modeleIGA._COORDS)[1]):
    if np.all(modeleIGA._COORDS[:,idx] == pt_coords):
        found = True
        break

if not found:
    sys.exit(-1)

print("yeti solution : ", SOL[idx,:])
print("reference solution : ", ref_aba)
print("error / component : ", np.abs((SOL[idx,:]-ref_aba)/ref_aba))

if np.abs((SOL[idx,:]-ref_aba)/ref_aba)[0] < 0.02:
    sys.exit(0)
else:
    sys.exit(-1)
