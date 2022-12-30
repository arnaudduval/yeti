# Copyright 2021 Arnaud Duval

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

# Python module
import numpy as np
import scipy.sparse as sp
import time

# IGA module
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip

# Selection of .INP and .NB file
# ------------------------------
FILENAME = 'chamfer'


# Creation of the IGA object
# --------------------------
modeleIGA = IGAparametrization(filename=FILENAME)

ti = time.time()

# REFINEMENT for Analysis
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

# domain 1 : chamfer
nb_deg[:, 0] = [0, 0, 0]
nb_ref[:, 0] = [3, 4, 3]
# domain 2 : tube
nb_deg[:, 1] = [1, 0, 1]
nb_ref[:, 1] = [2, 4, 3]
# domain 3 : base
nb_deg[:, 2] = [1, 0, 1]
nb_ref[:, 2] = [4, 4, 3]

additional_knots = {"patches": np.array([1]),
                    "1": np.array([]),
                    "2": np.array([]),
                    "3": np.array([0.25])
                    }

# Interfaces
nb_deg[:2, (2, 3, 5, 6, 7, 8)] = 0
nb_ref[:2, (3, 4)] = ([2, 4], [2, 4])
nb_ref[:2, (5, 6)] = ([2, 4], [2, 4])
nb_ref[:2, (7, 8)] = ([2, 4], [2, 4])

# Lgrge
nb_deg[:2, (9, 10, 11)] = 1
nb_ref[:2, 9] = [2, 4]
nb_ref[:2, 10] = [2, 4]
nb_ref[:2, 11] = [2, 4]

modeleIGA.refine(nb_ref, nb_deg, additional_knots)

modeleIGA._NBPINT[np.where(modeleIGA._ELT_TYPE == 'U00')] = \
    6**modeleIGA._dim.min()

# --
# STATIC STUDY

# MATRIX Assembly

ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

t1 = time.time()
data, row, col, Fb = \
    build_stiffmatrix(*modeleIGA.get_inputs4system_elemStorage())
Kside = sp.coo_matrix((data, (row, col)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside + Kside.transpose()
del Kside, data, row, col
print(('\n Time to build stiffness matrix : %.2f s' % (time.time() - t1)))

t1 = time.time()
Cdata, Crow, Ccol = cplg_matrix(*modeleIGA.get_inputs4cplgmatrix())
Cside = sp.coo_matrix((Cdata, (Crow, Ccol)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ctot = Cside + Cside.transpose()
del Cdata, Crow, Ccol, Cside
print(('\n Time to build coupling matrix  : %.2f s\n' % (time.time() - t1)))


# RESOLUTION Monolithique
t2 = time.time()
K2solve = Ktot[idof, :][:, idof]
C2solve = Ctot[idof, :][:, idof] * K2solve.max()
LU = sp.linalg.splu(K2solve + C2solve)
x = LU.solve(Fb[idof])
# x = np.zeros(idof.size)
# x = sp.linalg.spsolve(K2solve + C2solve,Fb[idof])
print(('\n Time for monolithique solving : %.2f s\n\n' % (time.time() - t2)))

# Postprocessing
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'coupling_mono', SOL.transpose(), nb_ref=np.array([3, 3, 3]),
    flag=np.array([True, False, False])))

print(('Total time for Mono analysis : %.2f s' % (time.time() - ti)))

# Postproc at interfaces
pp.generate_coupling_vtu(*modeleIGA.get_inputs4postproc_cplg_vtu(
    'interfaces_10', 10, SOL.transpose(), nb_ref=np.array([3, 3])))
pp.generate_coupling_vtu(*modeleIGA.get_inputs4postproc_cplg_vtu(
    'interfaces_11', 11, SOL.transpose(), nb_ref=np.array([3, 3])))
pp.generate_coupling_vtu(*modeleIGA.get_inputs4postproc_cplg_vtu(
    'interfaces_12', 12, SOL.transpose(), nb_ref=np.array([3, 3])))
