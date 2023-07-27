# Python module
import numpy as np
import scipy.sparse as sp
import sys
import time

#IGA module
from solver import PCPGortho
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

modeleIGA = IGAparametrization(filename='twoplatesDD')

ti = time.time()

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

p = 2
r = 2
# domains
nb_deg[:2,:2] = p
nb_ref[:2, 0] = r+1
nb_ref[:2, 1] = r
# curves
nb_ref[0,(2,3)] = r+2
# lgrge
nb_deg[0,4] = p
nb_ref[0,4] = r

modeleIGA.refine(nb_ref,nb_deg,additional_knots)
modeleIGA._NBPINT[ np.where(modeleIGA._ELT_TYPE == 'U00') ] = 6

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


# RESOLUTION Monolithique
t2 = time.time()
K2solve = Ktot[idof,:][:,idof]
C2solve = Ctot[idof,:][:,idof] * K2solve.max()
LU = sp.linalg.splu(K2solve + C2solve)
x  = LU.solve(Fb[idof])
#x = np.zeros(idof.size)
#x = sp.linalg.spsolve(K2solve + C2solve,Fb[idof])
print(('\n Time for monolithique solving : %.2f s\n\n' % (time.time() - t2)))

# Postprocessing
print(modeleIGA.get_inputs4solution(x))
# exit()
SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'coupling_mono',SOL.transpose(),nb_ref=np.array([3,3,3]),
    Flag=np.array([True,False,False])))

print(('Total time for Mono analysis : %.2f s' % (time.time() - ti)))