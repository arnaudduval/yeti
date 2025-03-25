"""
Linear analysis of a clamped square in 2D
Computation is made for different degrees and refinement levels
"""


import time

import numpy as np
import scipy.sparse as sp

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp



# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename='square')

# Refine modele
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

nb_ref[:, 0] = np.array([0, 0, 0])
nb_deg[:, 0] = np.array([0, 0, 0])

# Initial refinement (none)
modeleIGA.refine(nb_ref, nb_deg)

# Matrix assembly
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

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

print(K2solve.toarray())
print(K2solve.toarray().shape)
print(Fb[idof])

 # RESOLUTION : direct solver
t2 = time.time()
x  = sp.linalg.spsolve(K2solve, Fb[idof])
# LU = sp.linalg.splu(K2solve)
# x = LU.solve(Fb[idof])
print(('\n Time for direct solving :    %.2f s\n\n' % (time.time() - t2)))

print(x)

from solver import pseudoLUstep
import numpy as np
import scipy.sparse as sp

# A = sp.coo_matrix((np.array([1,2,1]),(np.array([0,1,3]),(np.array([0,1,3])))), shape=(4,4)).tocsc()
#A = sp.coo_matrix((np.array([1,2e-9,2,1e-9,1]),(np.array([0,1,2,3,4]),(np.array([0,1,2,3,4])))), shape=(5,5)).tocsc()

A = Ktot

LU_ = sp.linalg.splu(A,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0.,
                            options=dict(SymmetricMode=True))


print(('\n Matrice U post permutation :'))
print(LU_.U.toarray())
print(('\n Matrice L post permutation :'))
print(LU_.L.toarray())
print("\n")
print(type(LU_.L))

tol = 1e-8
Ukk = np.abs(LU_.U.diagonal())
ind = np.where(Ukk<tol*Ukk.max())[0]


print(('\n Ukk :'))
print(Ukk)
print(('\n ind :'))
print(ind)
print("\n")

LU = pseudoLUstep(A,tol=1.e-8)


print(('\n Matrice R :'))
print(LU.R)
print(LU.R.shape)
print(('\n Matrice U :'))
print(LU._LU.U)
print(LU._LU.U.shape)
print(('\n Matrice L :'))
print(LU._LU.L)
print(LU._LU.L.shape)


print('\n')
print(Ktot.toarray())
print(Ktot.toarray().shape)