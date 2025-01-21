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
Compute deflexion due to volume force
Compare with reference Abaqus result
"""

import sys
import time

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import IGAmanip as manip
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename='volume_force_3D')

# Refine modele
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

nb_deg[0, 0] = 1
nb_deg[1, 0] = 1
nb_deg[2, 0] = 1

nb_ref[0, 0] = 3
nb_ref[1, 0] = 4
nb_ref[2, 0] = 3

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

pp.generate_vtu_bezier(**modeleIGA.get_inputs4postproc_bezier(
        1,
        'fvol',
        SOL.transpose(),
        ))



disp = pp.evaldisp(*modeleIGA.get_inputs4evaldisp(
    SOL.transpose(), np.array([1., 1., 0.]), numpatch=1))


print(disp)

REF_SOLUTION = [-9.96639E-5, -1.01212, -6.42538]

error = np.linalg.norm(disp-REF_SOLUTION)/np.linalg.norm(REF_SOLUTION)

assert error < 2.e-3


