# Copyright 2021 Marie Guerder

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

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 16:27:57 2021

@author: mguerder
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

# =============================================================================
# Creation of the IGA object
# =============================================================================
# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename='test_embd_cube')

# Set arguments for model refinement
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
additional_knots = {"patches": np.array([]),
                    "1": np.array([]),
                    "2": np.array([]),
                    "3": np.array([])}

nb_ref[:, 1] = [1, 1, 1]    # Knot insertion, [xi, eta, zeta]

# Refine model
modeleIGA.refine(nb_ref, nb_deg, additional_knots)

# =============================================================================
# Static study
# =============================================================================
# Matrix assembly
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof] - 1

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
LU = sp.linalg.splu(K2solve)
x = LU.solve(Fb[idof])

print(('\n Time for direct solving :    %.2f s\n\n' % (time.time() - t2)))

# =============================================================================
# Postprocessing
# =============================================================================
# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

# VTU file generation
nb_ref_visu = np.array([1, 1, 1])       # Refinement for visu: [xi, eta, zeta]
output = np.array([True, True, False])  # Output type: [disp, stress, VM]
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'test_embd_cube',
    SOL.transpose(), nb_ref=3*nb_ref_visu, Flag=output))
modeleIGA.generate_vtk4controlMeshVisu(
    'test_embd_cube_map_cp', 0)
modeleIGA.generate_vtk4controlMeshVisu(
    'test_embd_cube_embded_cp', 1)

# =============================================================================
# Check results
# =============================================================================
ref_sol = np.loadtxt('ref_sol.txt')

error = np.sqrt(sum((x - ref_sol)**2.))

print(error)

if error > 1.e-3:
    sys.exit(-1)
else:
    sys.exit(0)
