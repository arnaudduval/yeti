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
Compute eigenfrequencies of a plate modeled with 3D solid elements
Results a compared with numerical reference (Abaqus)
"""

# Python module
import numpy as np
import scipy.sparse as sp
import sys

# IGA module
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from massmtrx import build_cmassmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# Base path for .inp and .NB files
FILENAME = 'plateVolume'

# Create IGA object
modeleIGA = IGAparametrization(filename=FILENAME)

# Refine modele
nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

nb_deg[0,0] = 1
nb_deg[1,0] = 1
nb_deg[2,0] = 1

nb_ref[0,0] = 5
nb_ref[1,0] = 5
nb_ref[2,0] = 1

# Initial refinement (none)
modeleIGA.refine(nb_ref,nb_deg)

# Build stiffness matrix
data,row,col,Fb = build_stiffmatrix( *modeleIGA.get_inputs4system_elemStorage() )
Kside = sp.coo_matrix((data,(row,col)), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot), 
                      dtype='float64').tocsc()
Ktot  = Kside + Kside.transpose()
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1
K2solve = Ktot[idof,:][:,idof]
del Kside,data,row,col,Ktot

# Build mass matrix
data,row,col = build_cmassmatrix( *modeleIGA.get_inputs4massmat() )
Mside = sp.coo_matrix((data,(row,col)), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot), 
                          dtype='float64').tocsc()
Mtot  = Mside + Mside.transpose()
M2solve = Mtot[idof,:][:,idof]
del Mside,data,row,col,Mtot

# Compute eigenvalues and eigenfrequencies
nb_frq = 10
vals, vecs = sp.linalg.eigsh(K2solve,k=nb_frq,M=M2solve,sigma=0.)
frq = np.sqrt(vals[:])/2./np.pi

# Save results in VTU file
Output = np.array([True, False, False])
for i in range(nb_frq):
    SOL,U = rsol.reconstruction( **modeleIGA.get_inputs4solution(vecs[:,i]) )
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
            'vib%0.2d' % i,SOL.transpose(),nb_ref=np.array([2,2,2]),flag=Output) )

# Reference eigenfrequencies computed in Abaqus
ref_frq = np.array([0.23088, 0.56354, 1.4185, 1.7967, 2.0466, 3.5630, 4.0492, 4.2456, 4.6715, 6.0879])

maxerror = max(abs(ref_frq[:]-frq[:])/ref_frq[:])
print(maxerror)

if maxerror > 0.01:
    sys.exit(-1)
else:
    sys.exit(0)
