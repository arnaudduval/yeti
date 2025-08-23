

# Python module
import numpy as np
import sys
import time

import scipy.sparse as sp
from scipy.sparse import linalg as sla

from preprocessing.igaparametrization import IGAmanip as manip
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

#IGA module
from preprocessing.igaparametrization import IGAparametrization, IGAmanip as manip
#from stiffmtrx_embedded import sys_linmat_lindef_static as build_stiffmatrix_emb
#import coupling.penlty_embeddedstruct_rot as cplg
#import bandstorage as bs
import postprocessing.postproc as pp
#import entities2param as tt
import reconstructionSOL as rsol

# Selection of .INP and .NB file
# ------------------------------
DIRECTORY='inputFiles/'
CASE=[]
CASE.append('stiffenedPanelEmbded') # ------ 0

EXEMPLE_NB = 0

FILENAME = DIRECTORY+CASE[EXEMPLE_NB]


# Creation of the IGA object
# --------------------------
modeleIGA = IGAparametrization(filename=FILENAME)


ti = time.time()
nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

p = 3
r = 5
# domains
nb_deg[:2,0] = [p,p]
nb_deg[:2,1] = [p,p]
nb_ref[:2,1] = [r,r]
nb_deg[:2,2] = [p,p]
nb_ref[:2,2] = [r+1,r-2]
# curves
nb_deg[0,(3,4,5,6)] = 1
nb_ref[0,(3,4,5,6)] = r+2
# lgrge
nb_deg[0,(7,8)] = [p,p-1]
nb_ref[0,(7,8)] = r

modeleIGA.refine(nb_ref,nb_deg)

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
np.savetxt('/home/agagnaire/yeti/temp/CsideU4', (Cside[idof,:][:,idof])[:,-4:].todense())

with open("/home/agagnaire/yeti/temp/CsideU4.txt","w") as file:
    file.write(str(Cside[idof,:][:,idof])) 
Ctot  = Cside + Cside.transpose()
del Cdata,Crow,Ccol,Cside
print(('\n Time to build coupling matrix  : %.2f s\n' % (time.time() - t1)))

print(idof.shape)

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
    'coupling_mono',SOL.transpose(),nb_ref=np.array([3,3,0]),
    Flag=np.array([True,False,False])))

print(('Total time for Mono analysis : %.2f s' % (time.time() - ti)))







