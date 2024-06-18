"""
Comparaison of U4 and U5 coupling methods in 2D case
Compare computed coupling matrices
"""


# Python module
import sys

import numpy as np
import scipy.sparse as sp

#IGA module
from preprocessing.igaparametrization import IGAparametrization
from coupling.cplgmatrix import cplg_matrix, cplg_matrixu5

modeleIGA_U4 = IGAparametrization(filename='twoplates_U4')
modeleIGA_U5 = IGAparametrization(filename='twoplates_U5')

# Models refinement parameters
P = 0
R = 0

# Refine U4 model
nb_deg = np.zeros((3,modeleIGA_U4.nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA_U4.nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

# domains
nb_deg[:2,:2] = P
nb_ref[:2, 0] = R
nb_ref[:2, 1] = R
# curves
nb_ref[0,(2,3)] = R
# lgrge
nb_deg[0,4] = P
nb_ref[0,4] = R

modeleIGA_U4.refine(nb_ref,nb_deg,additional_knots)
modeleIGA_U4._NBPINT[ np.where(modeleIGA_U4._ELT_TYPE == 'U00') ] = 3


# refine U5 model

nb_deg = np.zeros((3,modeleIGA_U5.nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA_U5.nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

# domains
nb_deg[:2,:2] = P
nb_ref[:2, 0] = R
nb_ref[:2, 1] = R
# lgrge
nb_ref[:,2] = np.array([R,0,0])
nb_deg[:,2] = np.array([R,0,0])

modeleIGA_U5.refine(nb_ref,nb_deg,additional_knots)

# Build coupling matrix for U4 model
Cdata,Crow,Ccol = cplg_matrix( *modeleIGA_U4.get_inputs4cplgmatrix() )
Cside_U4 = sp.coo_matrix((Cdata,(Crow,Ccol)),
                         shape=(modeleIGA_U4.nb_dof_tot, modeleIGA_U4.nb_dof_tot),
                         dtype='float64').tocsc()

# Build coupling matrix for U5 model
Cdata, Crow, Ccol = cplg_matrixu5( *modeleIGA_U5.get_inputs4cplgmatrixU5(integrationOrder=3) )

Cside_U5 = sp.coo_matrix((Cdata,(Crow,Ccol)),
                         shape=(modeleIGA_U5.nb_dof_tot, modeleIGA_U5.nb_dof_tot),
                         dtype='float64').tocsc()


# Slice to get components used in linear solver
ndof_U4 = modeleIGA_U4.nb_dof_free
idof_U4 = modeleIGA_U4.ind_dof_free[:ndof_U4]-1

ndof_U5 = modeleIGA_U5.nb_dof_free
idof_U5 = modeleIGA_U5.ind_dof_free[:ndof_U5]-1

error = np.linalg.norm((Cside_U4[idof_U4,:][:,idof_U4] - \
                        Cside_U5[idof_U5,:][:,idof_U5]).todense()) / \
    np.linalg.norm(Cside_U4[idof_U4,:][:,idof_U4].todense())

print(f"{error = }")

if error > 1.E-15:
    sys.exit(-1)
