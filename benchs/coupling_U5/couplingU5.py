# Copyright 2021 Arnaud Duval

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

# Python module
import numpy as np
from numpy.lib.arraysetops import intersect1d
import scipy.sparse as sp
import sys
import time

#IGA module
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5
from utils import gausspts
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip

FILENAME='testU5'

modeleIGA = IGAparametrization(filename=FILENAME)

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

nb_ref[:,0] = np.array([3,3,3])
nb_deg[:,0] = np.array([1,1,1])

nb_ref[:,1] = np.array([3,4,4])
nb_deg[:,1] = np.array([1,1,1])

nb_ref[:,2] = np.array([3,3,0])
nb_deg[:,2] = np.array([1,1,0])


modeleIGA.refine(nb_ref,nb_deg,additional_knots)

# STIFFNESS MATRIX
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

data, row, col, Fb = build_stiffmatrix( *modeleIGA.get_inputs4system_elemStorage() )
Kside = sp.coo_matrix((data, (row,col)), shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                        dtype='float64').tocsc()
Ktot = Kside+Kside.transpose()
del Kside, data, row, col

print("stiffness matrix assembly done")

# COUPLING MATRIX
Cdata, Crow, Ccol = cplg_matrixU5( *modeleIGA.get_inputs4cplgmatrixU5() )
Cside = sp.coo_matrix((Cdata,(Crow,Ccol)), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ctot  = Cside + Cside.transpose()
del Cdata,Crow,Ccol,Cside

print("Coupling matrix assembly done")

# Monolithic SOLVE
t2 = time.time()
K2solve = Ktot[idof,:][:,idof]
C2solve = Ctot[idof,:][:,idof] * K2solve.max()
LU = sp.linalg.splu(K2solve + C2solve)
x  = LU.solve(Fb[idof])
t3 = time.time()

print("Resolution done ", t3-t2, " seconds")

# Postprocessing
SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'couplingU5',SOL.transpose(),nb_ref=np.array([1,1,1]),
    Flag=np.array([True,True,False])))

# Verify result
xtest = np.array([6.,0.,2.])
search = np.arange(np.shape(modeleIGA._COORDS)[1])
for i in range(3):
    search = intersect1d(search, np.where(modeleIGA._COORDS[i,:] == xtest[i]))

# Reference FE solution computed with Abaqus
ref = -1.136E-2
print("Reference solution : ", ref)
print("Computed solution : ", SOL[search[0],2])

error = ( SOL[search[0],2] - ref ) / ref

print(f"{error = }")

if abs(error) > 0.02:
    sys.exit(-1)
else:
    sys.exit(0)
