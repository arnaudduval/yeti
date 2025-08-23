# Python module
import numpy as np
import scipy.sparse as sp
from scipy.linalg import null_space, svd, eig, pinv
import sys
import time

#IGA module
from solver import PCPGortho
from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import IGAmanip as manip
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrix_collocation as cplg_matrix_collocation
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from utils import rigid_motion

# modeleIGA = IGAparametrization(filename='rectangle')
modeleIGA = IGAparametrization(filename='plate_disk_monopatch')

ti = time.time()

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}


p = 1 #p-1
r = 3
# domains
nb_deg[:2,0] = p
nb_ref[0, 0] = r+1  #r+1
nb_ref[1, 0] = r
modeleIGA.refine(nb_ref,nb_deg,additional_knots)

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}
modeleIGA.refine(nb_ref,nb_deg,additional_knots)

additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}





# modeleIGA.refine(nb_ref,nb_deg,additional_knots)
# --
# STATIC STUDY

# MATRIX Assembly

ndof = modeleIGA._nb_dof_free
print(f'{ndof = }')
idof = modeleIGA._ind_dof_free[:ndof]-1


t1 = time.time()


data, row, col, Fb = build_stiffmatrix( *modeleIGA.get_inputs4system_elemStorage() )
Kside = sp.coo_matrix((data, (row,col)), shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                        dtype='float64').tocsc()
Ktot = Kside+Kside.transpose()
del Kside, data, row, col

print("stiffness matrix assembly done")
print(f'{Ktot.shape = }')
print(f'{modeleIGA._nnode = }')
print(f'{modeleIGA._nb_dof_tot = }')


# Monolithic SOLVE
t2 = time.time()
K2solve = Ktot[idof,:][:,idof]

# print(null_space(K2solve.todense()+ C2solve.todense(),1e-2))

# kernel = null_space(K2solve.todense()+ C2solve.todense(),1e-3)
# print(f'{K2solve.max()=} ')

# eigenvalues, eigenvectors = eig(Ktot[:8,:8].todense())
# print('-------')
# print(eigenvalues)
# print('-------')

# eigenvalues, eigenvectors = eig(K2solve[:12,:12].todense())
# print('-------')
# print(eigenvalues)
# print('-------')

# for imode in range (kernel.shape[1]):
#     print(kernel[-6:, imode])
#     print('-------')



# U, s, Vh = svd(K2solve.todense()+ C2solve.todense(),1e-3)

# print(s)

LU = sp.linalg.splu(K2solve)
x  = LU.solve(Fb[idof])

# x = sp.linalg.spsolve(K2solve + C2solve, Fb[idof])



print(('\n Time for monolithique solving : %.2f s\n\n' % (time.time() - t2)))

# Postprocessing
# print(modeleIGA.get_inputs4solution(x))
# exit()

# ur = rigid_motion(*modeleIGA.get_inputs4rigid_motion(1, np.array([0.,0.]), np.array([5.,2.5]), 1.))
# ur = ur + rigid_motion(*modeleIGA.get_inputs4rigid_motion(2, np.array([0.,0.]), np.array([5.,2.5]), 1.))
# SOLr = np.reshape(ur,(modeleIGA._nb_dof_tot//2, 2))

# np.savetxt('/home/agagnaire/yeti/temp/SOLr', SOLr)

SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'mono',SOL.transpose(),nb_ref=np.array([3,3,3]),
    Flag=np.array([True,False,False])))

print(('Total time for Mono analysis : %.2f s' % (time.time() - ti)))

print(np.max(SOL))

# # - factorization
# from solver import pseudoLUstep

# # print(Ktot[:8,:8].todense())
# np.savetxt("/home/agagnaire/yeti/temp/K1K2.txt",  Ktot.todense(), fmt='%.4e')


# LU1 = pseudoLUstep(Ktot, tol=1.e-8)

# LU = pseudoLUstep(K2solve,tol=1.e-8)

# moore = (pinv(K2solve.todense()))
# mo = sp.csc_matrix(moore)







