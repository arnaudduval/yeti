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

# -*- coding: utf-8 -*-

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
# Selection of .INP and .NB file
# =============================================================================
DIRECTORY = 'inputFiles/'
CASE = []
CASE.append('AIRBUSCROR')  # --------------------- 0
CASE.append('hemisphere_shell')  # --------------- 1
CASE.append('HorseshoeD2S335')  # ---------------- 2
CASE.append('multipatch')  # --------------------- 3
CASE.append('Ring3DD221S550')  # ----------------- 4
CASE.append('beam')  # --------------------------- 5
CASE.append('QuarterRing2DD2S20')  # ------------- 6
CASE.append('beam2patch')  # --------------------- 7
CASE.append('hemisphere_solid')  # --------------- 8
CASE.append('scordelis-Lo')  # ------------------- 9
CASE.append('beamKL')  # ------------------------- 10
CASE.append('Lshaped-plate')  # ------------------ 11
CASE.append('plate')  # -------------------------- 12
CASE.append('arc_bot')  # ------------------------ 13
CASE.append('beamShear')  # ---------------------- 14
CASE.append('stiffenedPanel')  # ----------------- 15
CASE.append('TubeRectU1')  # --------------------- 16
CASE.append('beamEmbded')  # --------------------- 17
CASE.append('plaqueTrouee_embded')  # ------------ 18
CASE.append('plaqueTrouee')  # ------------------- 19
CASE.append('Contour_v1')  # --------------------- 20
CASE.append('test')  # --------------------------- 21
CASE.append('curvedwallstrangeStiff')  # --------- 22
CASE.append('twoplatesDD')  # -------------------- 23
CASE.append('wingraw2')  # ----------------------- 24
CASE.append('aero_bot1')  # ---------------------- 25
CASE.append('plateVolume')  # -------------------- 26
CASE.append('NR4L_20200115')  # ------------------ 27
CASE.append('NR4L_20200120')  # ------------------ 28
CASE.append('volblade')  # ----------------------- 29
CASE.append('DA36_2020_02_12')  # ---------------- 30
CASE.append('NR4L_20200114')  # ------------------ 31
CASE.append('pavecentrif')  # -------------------- 32
CASE.append('DA36_2020_03_10')  # ---------------- 33
CASE.append('beam_2')  # ------------------------- 34
CASE.append('beam2D')  # ------------------------- 35
CASE.append('beam_3')  # ------------------------- 36
CASE.append('NR4L_2020_06_25_pres')  # ----------- 37
CASE.append('verif_pres')  # --------------------- 38
CASE.append('NR4L_2020_09_21_pres_centrif')  # --- 39
CASE.append('pres_uni_face')  # ------------------ 40
CASE.append('pres_uni_field')  # ----------------- 41
CASE.append('pres_non_uni_field')  # ------------- 42
CASE.append('pres_uni_field_new')  # ------------- 43
CASE.append('pres_sin_field')  # ------------- 44
CASE.append('pres_sin_field_enc')  # ------------- 45
CASE.append('poutre-triangle-1')  # ------------- 46


EXEMPLE_NB = 46

FILENAME = DIRECTORY + CASE[EXEMPLE_NB]

# =============================================================================
# Creation of the IGA object
# =============================================================================
# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename=FILENAME)


# Set arguments for model refinement
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
additional_knots = {"patches": np.array([]),
                    "1": np.array([]),
                    "2": np.array([]),
                    "3": np.array([])}

nb_deg[:, 0] = [1, 1, 1]    # Degree elevation, [xi, eta, zeta]
nb_ref[:, 0] = [3, 2, 5]    # Knot insertion, [xi, eta, zeta]
# additional_knots['patches'] = np.array([0])    # Specific knot insertion
# additional_knots['2'] = np.array([1./3.,2./3.],dtype=np.float64)
# additional_knots['3'] = np.array([1./3.,2./3.],dtype=np.float64)

# Refine model
modeleIGA.refine(nb_ref, nb_deg, additional_knots)




# =============================================================================
# Static study
# =============================================================================
# Matrix assembly
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

t1 = time.time()

data, row, col, Fb = build_stiffmatrix(
                        *modeleIGA.get_inputs4system_elemStorage())
#exit()
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

# RESOLUTION : conjugate gradient
# Ddata  = 1./K2solve.diagonal()
# jacobi = sp.dia_matrix((Ddata[np.newaxis, :], np.array([0])),
#                        shape=K2solve.shape,
#                        dtype=np.float64)
# t2 = time.time()
# xI = sp.linalg.cg(K2solve, Fb[idof],tol=1.e-6, M=jacobi)
# print ('\n Time for iterative solving : %.2f s\n\n' % (time.time() - t2))

# =============================================================================
# Postprocessing
# =============================================================================
# Solution reconstruction
SOL, u = rsol.reconstruction(*modeleIGA.get_inputs4solution(x))

# Recherche d'un point pour verifier la solution
# Méthode pas optimale
pt_coords = np.array([3.,1.,20.])

for idx in range(np.shape(modeleIGA._COORDS)[1]):
    if(np.all(modeleIGA._COORDS[:,idx] == pt_coords)):
        break

print(SOL[idx,:])
print(modeleIGA._COORDS[:,idx])



# VTU file generation
nb_ref_visu = np.array([1, 0, 1])       # Refinement for visu: [xi, eta, zeta]
output = np.array([True, True, False])  # Output type: [disp, stress, VM]
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'linear_analysis_' + CASE[EXEMPLE_NB],
    SOL.transpose(), nb_ref=3*nb_ref_visu, Flag=output))
#modeleIGA.generate_vtk4controlMeshVisu(
#    'linear_analysis_' + CASE[EXEMPLE_NB], 0)

# pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
#     'pres_field_' + CASE[EXEMPLE_NB],
#     [modeleIGA._pres_field]*3, nb_ref=3*nb_ref_visu, Flag=output))
