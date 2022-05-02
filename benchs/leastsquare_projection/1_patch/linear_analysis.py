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

# Single domain with 1 patch

# Python module
import numpy as np
import scipy.sparse as sp
import sys
import time

#IGA module
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip

FILENAME = 'plateWithRoundHole'

# Creation of the IGA object
# --------------------------
modeleIGA = IGAparametrization(filename=FILENAME)

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

nb_deg[:,0] = [0,0,0]
nb_ref[:,0] = [3,3,0]

modeleIGA.refine(nb_ref,nb_deg,additional_knots)


# --
# STATIC STUDY

# MATRIX Assembly

ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

t1 = time.time()
data,row,col,Fb = build_stiffmatrix( *modeleIGA.get_inputs4system_elemStorage() )
Kside = sp.coo_matrix((data,(row,col)), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot), 
                      dtype='float64').tocsc()
Ktot  = Kside + Kside.transpose()
K2solve = Ktot[idof,:][:,idof]
del Kside,data,row,col
print(('\n Time to build stiffness matrix : %.2f s' % (time.time() - t1)))


# SOLVING : direct solver
t2 = time.time()
LU = sp.linalg.splu(K2solve)
x  = LU.solve(Fb[idof])
print(('\n Time for direct solving :    %.2f s\n\n' % (time.time() - t2)))

# Standard postprocessing
# -----------------------
SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
# Stress and strain post processing is deactivated due to a singular point
# leading to NaN values
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'linear_analysis',SOL.transpose(),nb_ref=3*np.array([0,0,0]),
    Flag=np.array([True,False,False])))

# Postprocessing with least square projection
# -------------------------------------------
# Compute gram matrix
data, row, col = pp.build_cgrammatrix( *modeleIGA.get_inputs4grammat())
# RHS vector
rhs = pp.compute_svars_solid_rhs( *modeleIGA.get_inputs4svarsrhs(SOL.transpose()))

size = max(col)+1
Gside = sp.coo_matrix((data,(row,col)), shape=(size,size), 
                      dtype='float64').tocsc()
Gtot  = Gside + Gside.transpose()

del Gside,data,row,col
LU = sp.linalg.splu(Gtot)

projected_svars = LU.solve(rhs.transpose())

pp.generate_proj_vtu(*modeleIGA.get_inputs4proj_vtu('least_square', SOL.transpose(), projected_svars.transpose(), nb_ref=3*np.array([1,1,0])))
