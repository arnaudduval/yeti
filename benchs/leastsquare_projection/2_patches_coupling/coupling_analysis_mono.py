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

# Coupled domain with monolithic solving

# Python module
import numpy as np
import scipy.sparse as sp
import sys
import time

#IGA module
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip

FILENAME = 'plateWithRoundHoleCPLG'


# Creation of the IGA object
# --------------------------
modeleIGA = IGAparametrization(filename=FILENAME)

ti = time.time()

# REFINEMENT for Analysis
nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

p = 2   # analysis degree
r = 2   # refinement level

# domain 1 (initially at degree 2)
nb_deg[:, 0] = [p - 2, p - 2, p - 2]
nb_ref[:, 0] = [r, r, 0]
# domain 2 (initially at degree 2)
nb_deg[:, 1] = [p - 2, p - 2, p - 2]
nb_ref[:, 1] = [r+1, r+1, 0]
additional_knots = {"patches":np.array([1]),"1":np.array([]),
                    "2":np.array([0.3]),"3":np.array([])} 
# interface
nb_ref[:1, 2] = [r+2]
nb_ref[:1, 3] = [r+2]
# Lagrange field (initially at degree 0)
nb_deg[:1, 4] = p - 1
nb_ref[:1, 4] = r


modeleIGA.refine(nb_ref,nb_deg,additional_knots)

modeleIGA._NBPINT[ np.where(modeleIGA._ELT_TYPE == 'U00') ] = 6**modeleIGA._dim.min()
print(modeleIGA._dim.min())
print(modeleIGA._NBPINT)

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
del Kside,data,row,col
print(('\n Time to build stiffness matrix : %.2f s' % (time.time() - t1)))

t1 = time.time()
Cdata,Crow,Ccol = cplg_matrix( *modeleIGA.get_inputs4cplgmatrix() )
Cside = sp.coo_matrix((Cdata,(Crow,Ccol)), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ctot  = Cside + Cside.transpose()
del Cdata,Crow,Ccol,Cside
print(('\n Time to build coupling matrix  : %.2f s\n' % (time.time() - t1)))


# Monolithic solving
t2 = time.time()
K2solve = Ktot[idof,:][:,idof]
C2solve = Ctot[idof,:][:,idof] * K2solve.max()
LU = sp.linalg.splu(K2solve + C2solve)
x  = LU.solve(Fb[idof])
print(('\n Time for monolithique solving : %.2f s\n\n' % (time.time() - t2)))

# Standard Postprocessing
# -----------------------
SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'coupling_mono',SOL.transpose(),nb_ref=np.array([3,3,3]),
        flag=np.array([True,False,False])))

print(('Total time for Mono analysis : %.2f s' % (time.time() - ti)))
    
# Least square projection post processing
# ---------------------------------------
# Gram matrix
data, row, col = pp.build_cgrammatrix( *modeleIGA.get_inputs4grammat())
# RHS vector
rhs = pp.compute_svars_solid_rhs( *modeleIGA.get_inputs4svarsrhs(SOL.transpose()))

size = modeleIGA._nb_cp
        
# get max node index for U1 elements
# WARNING : Assumption : U1 elements are given first
max_ind = 0
for i in np.where(modeleIGA._ELT_TYPE == 'U1')[0]:
    max_ind = max(max_ind, np.max(modeleIGA._IEN[i]))
        
# Add identity for Gram matrix coefficients not related to U1 elements
data = np.concatenate((data, 0.5*np.ones(size-max_ind)), axis=None)
row = np.concatenate((row, np.arange(max_ind, size)), axis=None)
col = np.concatenate((col, np.arange(max_ind, size)), axis=None)
        
Gside = sp.coo_matrix((data,(row,col)), shape=(size,size), 
                      dtype='float64').tocsc()
Gtot  = Gside + Gside.transpose()

del Gside,data,row,col
LU = sp.linalg.splu(Gtot)

projected_svars = LU.solve(rhs.transpose())

pp.generate_proj_vtu(*modeleIGA.get_inputs4proj_vtu('coupling_mono_least_square', SOL.transpose(), projected_svars.transpose(), nb_ref=3*np.array([1,1,0])))
