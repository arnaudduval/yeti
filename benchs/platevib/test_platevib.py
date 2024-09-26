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
Results a compared with numerical reference
"""

import os
import numpy as np
import scipy.sparse as sp
import sys

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from yeti_iga.massmtrx import build_cmassmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp

# Creation of the IGA object
script_dir = os.path.dirname(os.path.realpath(__file__))
modeleIGA = IGAparametrization(filename=f'{script_dir}/plateVolume')

# Refine modele
nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

nb_deg[0,0] = 1
nb_deg[1,0] = 1
nb_deg[2,0] = 1

nb_ref[0,0] = 2
nb_ref[1,0] = 2
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
nb_frq = 5
vals, _ = sp.linalg.eigsh(K2solve, k=nb_frq, M=M2solve, sigma=0.)
frq = np.sqrt(vals[:])/2./np.pi


# Reference eigenfrequencies
ref_frq = np.array([0.2548, 0.6217, 1.9546, 2.3870, 2.8585])

assert np.allclose(frq, ref_frq, rtol=1.e-3)
