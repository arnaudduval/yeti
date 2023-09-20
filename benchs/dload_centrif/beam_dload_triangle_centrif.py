# Copyright 2023 Arnaud Duval

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
A beam subjected to centrifugal force and distributed triangular pressure
"""

import sys
import time
import copy

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import IGAmanip as manip
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# Read data and create IGAparametrization objects
# Distributed pressure + centrifugal force
model_iga_pc = IGAparametrization(filename='beam-triangle-centrif')
# Distributed pressure only
model_iga_p = IGAparametrization(filename='beam-triangle')
# Centrifugal force only
model_iga_c = IGAparametrization(filename='beam-centrif')

# Refine models (single patch and 3D for all)
nb_deg = np.zeros((3, 1), dtype=np.intp)
nb_ref = np.zeros((3, 1), dtype=np.intp)

nb_deg[:, 0] = np.array([1, 1, 1])
nb_ref[:, 0] = np.array([2, 2, 5])

# Initial refinement
model_iga_pc.refine(nb_ref, nb_deg)
model_iga_p.refine(nb_ref, nb_deg)
model_iga_c.refine(nb_ref, nb_deg)


def matrix_assembly(model_iga):
    """
    Assemble stiffness matrix and force vector for a given IGA model
    """
    ndof = model_iga._nb_dof_free
    idof = model_iga._ind_dof_free[:ndof]-1
    data, row, col, Fb = build_stiffmatrix(
                        *model_iga.get_inputs4system_elemStorage())
    Kside = sp.coo_matrix((data, (row, col)),
                          shape=(model_iga._nb_dof_tot,
                          model_iga._nb_dof_tot),
                          dtype='float64').tocsc()
    Ktot = Kside + Kside.transpose()

    del Kside, data, row, col

    return copy.copy(Ktot), copy.copy(Fb), copy.copy(idof)


# Assemble stiffness matrices and force vectors
Ktot_pc, Fb_pc, idof_pc = matrix_assembly(model_iga_pc)
Ktot_p, Fb_p, idof_p = matrix_assembly(model_iga_p)
Ktot_c, Fb_c, idof_c = matrix_assembly(model_iga_c)

print(np.linalg.norm(Fb_pc - Fb_p - Fb_c))


def resolution_and_output(model_iga, Ktot, Fb, idof, name):
    """
    Solve system and save results in VTU format
    """

    K2solve = Ktot[idof, :][:, idof]
    LU = sp.linalg.splu(K2solve)
    x = LU.solve(Fb[idof])
    SOL, u = rsol.reconstruction(**model_iga.get_inputs4solution(x))
    pp.generatevtu(*model_iga.get_inputs4postprocVTU(
        name, SOL.transpose(), nb_ref=np.array([1, 1, 1]),
        Flag=np.array([True, True, False])))


resolution_and_output(model_iga_pc, Ktot_pc, Fb_pc, idof_pc,
                      'presfield+centrif')
resolution_and_output(model_iga_p, Ktot_p, Fb_p, idof_p,
                      'presfield')
resolution_and_output(model_iga_c, Ktot_c, Fb_c, idof_c,
                      'centrif')

exit()


# RESOLUTION : direct solver
t2 = time.time()
# x  = sp.linalg.spsolve(K2solve, Fb[idof])


print(('\n Time for direct solving :    %.2f s\n\n' % (time.time() - t2)))

# Solution reconstruction








exit()

# Recherche d'un point pour verifier la solution
# Méthode pas optimale
pt_coords = np.array([20.,0.,0.])

found = False
for idx in range(np.shape(modeleIGA._COORDS)[1]):
    if(np.all(modeleIGA._COORDS[:,idx] == pt_coords)):
        found = True
        break

if not found:
    sys.exit(-1)

print(SOL[idx,2])

# Analytic solution

b = 3.      # Base of section
h = 1.      # height of section
E = 210000. # Yound Modulus
L = 20.     # beam length
q0 = 1000.  # Load
I = b * h**3. / 12. #Inertia of section

# max deflexion
f = - q0 * b * L**4. / (30. * E * I)
print(f)

error = ( SOL[idx,2] - f ) / f
print(error)

if abs(error) > 0.02:
    sys.exit(-1)
else:
    sys.exit(0)
