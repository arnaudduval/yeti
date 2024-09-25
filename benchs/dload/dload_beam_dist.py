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
modeleIGA = IGAparametrization(filename='beam-dist')

# Refine modele
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

nb_deg[0, 0] = 1
nb_deg[1, 0] = 1
nb_deg[2, 0] = 1

nb_ref[0, 0] = 2
nb_ref[1, 0] = 2
nb_ref[2, 0] = 5

# Initial refinement (none)
modeleIGA.refine(nb_ref, nb_deg)

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
# x  = sp.linalg.spsolve(K2solve, Fb[idof])
LU = sp.linalg.splu(K2solve)
x = LU.solve(Fb[idof])
print(('\n Time for direct solving :    %.2f s\n\n' % (time.time() - t2)))

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
# Recherche d'un point pour verifier la solution
# Méthode pas optimale
pt_coords = np.array([20., 0., 0.])

found = False
for idx in range(np.shape(modeleIGA._COORDS)[1]):
    if(np.all(modeleIGA._COORDS[:, idx] == pt_coords)):
        found = True
        break

if not found:
    sys.exit(-1)

print(SOL[idx, 2])

# Analytic solution

b = 3.          # Base of section
h = 1.          # Height of section
E = 210000.     # Young Modulus
L = 20.         # Beam length
q0 = 1000.      # Load
I = b * h**3. / 12.     # Inertia of section

# max deflexion
f = - q0 * b * L**4. / (8. * E * I)
print(f)

error = ( SOL[idx,2] - f ) / f
print(error)

if abs(error) > 0.02:
    sys.exit(-1)
else:
    sys.exit(0)
