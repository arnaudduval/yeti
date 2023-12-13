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
Test computation time for different problem sizes
"""

import sys
import time
import os

import numpy as np
import scipy.sparse as sp

from copy import deepcopy

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import IGAmanip as manip
print('R1')
from stiffmtrx_elemstorage_omp import sys_linmat_lindef_static_omp as build_stiffmatrix_omp
print('R2')
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# Read data and create initial model
iga_model_ini = IGAparametrization(filename='beam-dist')
modeleIGA = deepcopy(iga_model_ini)

# Refine modele
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_deg[:, 0] = np.array([0, 0, 0])
nb_ref[:, 0] = np.array([1, 1, 0])
modeleIGA.refine(nb_ref, nb_deg)

# Matrix assembly
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

# Sequential
t0 = time.time()
data, row, col, Fb = build_stiffmatrix(
                    *modeleIGA.get_inputs4system_elemStorage())
t1 = time.time()
Kside = sp.coo_matrix((data, (row, col)),
                shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                dtype='float64').tocsc()

t2 = time.time()
# Multithreaded
modeleIGA2 = deepcopy(modeleIGA)
t3 = time.time()
dataMP, rowMP, colMP, FbMP = build_stiffmatrix_omp(
                    *modeleIGA2.get_inputs4system_elemStorage())
t4 = time.time()

KsideMP = sp.coo_matrix((dataMP, (rowMP, colMP)),
                shape=(modeleIGA2._nb_dof_tot, modeleIGA2._nb_dof_tot),
                dtype='float64').tocsc()

t5 = time.time()

print(np.linalg.norm(dataMP))

print("Erreur calcul K : ", sp.linalg.norm(Kside - KsideMP))
print("Acceleration factor : ", (t1-t0)/(t4-t3))

# print(data, dataMP)


exit()




# Ktot = Kside + Kside.transpose()
# t3 = time.time()
# K2solve = Ktot[idof, :][:, idof]
# t4 = time.time()
# del Kside, data, row, col
# t5 = time.time()
# x = sp.linalg.spsolve(K2solve, Fb[idof], use_umfpack=umfpack)




ref = {}
ref['initial'] = ([0, 0, 0],[0, 0, 0])
# ref['p2'] = ([1, 1, 1],[0, 0, 0])
# ref['p3'] = ([2, 2, 2],[0, 0, 0])
# ref['p4'] = ([3, 3, 3],[0, 0, 0])

# ref['p2r01'] = ([1, 1, 1],[1, 0, 0])
# ref['p2r02'] = ([1, 1, 1],[1, 1, 0])
ref['p2r03'] = ([1, 1, 1],[1, 1, 1])
# ref['p2r04'] = ([1, 1, 1],[2, 1, 1])
# ref['p2r05'] = ([1, 1, 1],[2, 2, 1])
# ref['p2r06'] = ([1, 1, 1],[2, 2, 2])
# ref['p2r07'] = ([1, 1, 1],[3, 2, 2])
ref['p2r08'] = ([1, 1, 1],[3, 3, 2])
# ref['p2r09'] = ([1, 1, 1],[3, 3, 3])
# ref['p2r10'] = ([1, 1, 1],[4, 3, 3])
# ref['p2r11'] = ([1, 1, 1],[4, 4, 3])
ref['p2r12'] = ([1, 1, 1],[4, 4, 4])
# ref['p2r13'] = ([1, 1, 1],[5, 4, 4])
# ref['p2r14'] = ([1, 1, 1],[5, 5, 4])
ref['p2r15'] = ([1, 1, 1],[5, 5, 5])

for umfpack in [True, False]:
    for r in ref.keys():
        # Copy model
        modeleIGA = deepcopy(iga_model_ini)
        # Refine modele
        nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
        nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

        print(ref[r][0], ref[r][1])

        nb_deg[:, 0] = ref[r][0]
        nb_ref[:, 0] = ref[r][1]



        # Initial refinement (none)
        modeleIGA.refine(nb_ref, nb_deg)

        # Matrix assembly
        ndof = modeleIGA._nb_dof_free
        idof = modeleIGA._ind_dof_free[:ndof]-1

        t0 = time.time()
        data, row, col, Fb = build_stiffmatrix(
                            *modeleIGA.get_inputs4system_elemStorage())

        t1 = time.time()
        Kside = sp.coo_matrix((data, (row, col)),
                        shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                        dtype='float64').tocsc()

        t2 = time.time()
        Ktot = Kside + Kside.transpose()
        t3 = time.time()
        K2solve = Ktot[idof, :][:, idof]
        t4 = time.time()
        del Kside, data, row, col
        t5 = time.time()
        x = sp.linalg.spsolve(K2solve, Fb[idof], use_umfpack=umfpack)
        t6 = time.time()
        with open('time.txt', 'a') as f:
            print(t1 - t0, t2 - t1, t3 - t2, t4 - t3, t5 - t4, t6 - t5, ndof, os.environ["OMP_NUM_THREADS"], umfpack)
            f.write(f'{t1 - t0}\t{t2 - t1}\t{t3 - t2}\t{t4 - t3}\t{t5 - t4}\t{t6 - t5}\t{ndof}\t{os.environ["OMP_NUM_THREADS"]}\t{umfpack}\n')

exit()


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
