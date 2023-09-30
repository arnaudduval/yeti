# Copyright 2023 Arnaud Duval

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

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Compute deflection due to 2 triangular pressure fields applied on 2 different
faces of a beam
Compare with the combination of analytic results
"""

import sys
import time

import numpy as np
import scipy.sparse as sp

import reconstructionSOL as rsol
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix

# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename='beam-dloadF2F4')

# Refine modele
nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)

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

# RESOLUTION : direct solver
t2 = time.time()
# x  = sp.linalg.spsolve(K2solve, Fb[idof])
LU = sp.linalg.splu(K2solve)
x = LU.solve(Fb[idof])

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
# Recherche d'un point pour verifier la solution
# Méthode pas optimale
pt_coords = np.array([20., 0., 0.])

found = False
for idx in range(np.shape(modeleIGA.coords)[1]):
    if np.all(modeleIGA.coords[:, idx] == pt_coords):
        found = True
        break

if not found:
    sys.exit(-1)

print('disp z : ', SOL[idx, 2])
print('disp y : ', SOL[idx, 1])

# Analytic solution for deflection along z axis

b = 3.          # Base of section
h = 1.          # Height of section
E = 210000.     # Yound Modulus
L = 20.         # Beam length
q0 = 1000.      # Load
inertia = b * h**3. / 12.     # Inertia of section

# max deflexion
fz = - q0 * b * L**4. / (30. * E * inertia)
print('disp z anaytic : ', fz)

# Analytic solution for defection along y axis
h = 3.          # Base of section
b = 1.          # Height of section
E = 210000.     # Yound Modulus
L = 20.         # Beam length
q0 = 666.      # Load
inertia = b * h**3. / 12.     # Inertia of section

# max deflexion
fy = - 11. * q0 * b * L**4. / (120. * E * inertia)
print('disp y anaylitic : ', fy)

error_z = (SOL[idx, 2] - fz) / fz
print('Error on z displacement: ', error_z)
error_y = (SOL[idx, 1] - fy) / fy
print('Error on y displacement: ', error_y)

if abs(error_z) > 0.02 or abs(error_y) > 0.02:
    sys.exit(-1)
else:
    sys.exit(0)
