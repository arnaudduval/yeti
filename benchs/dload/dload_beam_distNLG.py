"""
Compute deflexion due to distributed pressure field on a beam with nonlinear geomtric 
"""

import sys
import time

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import IGAmanip as manip
from tanmtrx_elemStorageNLG import build_tanmatrix 
from extvec_elemStorageNLG import build_extvector
from intvec_elemStorageNLG import compute_intvector
import reconstructionSOL as rsol
import postprocessing.postproc as pp

eps = 1E-08
nb_itMAX = 15
nb_loadF = 1

# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename='beam-dist')

# Refine modele
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

nb_deg[0, 0] = 0
nb_deg[1, 0] = 0
nb_deg[2, 0] = 0

nb_ref[0, 0] = 2
nb_ref[1, 0] = 2
nb_ref[2, 0] = 5

# Initial refinement (none)
modeleIGA.refine(nb_ref, nb_deg)

# Matrix assembly
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

# Compute loda vector 
Fext_tot = build_extvector(*modeleIGA.get_inputs4system_elemStorage())

DltF = Fext_tot.copy()
DltF[:] = DltF[:]/nb_loadF

# Initial guess and residual
Uk = np.zeros(modeleIGA._nb_dof_tot,dtype='float64')
Uk_new = np.zeros(modeleIGA._nb_dof_tot,dtype='float64')
Fext = np.zeros(modeleIGA._nb_dof_tot,dtype='float64')

for load in range(nb_loadF):
    bk = np.zeros(modeleIGA._nb_dof_tot,dtype='float64')
    Fext[:] += DltF[:]
    DltF0 = DltF.copy()
    normF = np.linalg.norm(Fext)
    Rk = np.zeros(modeleIGA._nb_dof_tot,dtype='float64')

    print(('Load n° %i' % (load)))
    cond = True
    it = 1
    while cond:

        t2 = time.time()
        # Compute tangent matrix 
        data, row, col = build_tanmatrix(
                                *modeleIGA.get_inputs4tanmtrx(Uk))

        Kside = sp.coo_matrix((data, (row, col)),
                            shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                            dtype='float64').tocsc()
        Ktot = Kside + Kside.transpose()
        KT = Ktot[idof, :][:, idof]
        del Kside, data, row, col
        # print(('Itération n° %i ; Time to build stiffness matrix : %.2f s' % (it,time.time() - t2)))
        bk[:] = Rk[:] + DltF0[:]
        
        # RESOLUTION : direct solver
        t3 = time.time()
        LU = sp.linalg.splu(KT)
        Dlt_Uk = LU.solve(bk[idof])
        # print(('Time for direct solving :    %.2f s' % (time.time() - t3)))
        Uk_new[idof] = Dlt_Uk[:]
        Uk[:] = Uk[:] + Uk_new[:]

        # Compute internal efforts and update residual
        Fint = compute_intvector(*modeleIGA.get_inputs4tanmtrx(Uk))
        Rk[idof] = Fext[idof]-(Fint[idof])
        
        it = it +1
        residual = np.linalg.norm(Rk)/(1+normF)
        print(('Iteration n° %i ; Residual :    %.2E \n' % (it,residual)))
        cond = (residual>eps and it< nb_itMAX)
        DltF0[:] = 0.
    

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(Uk[idof]))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU("dload_1",SOL.T))
