# Python module
import numpy as np
import scipy.sparse as sp
from scipy.linalg import null_space, svd, eig, pinv, svdvals
import sys
import time
import matplotlib.pyplot as plt

#IGA module
from solver import PCPGortho
from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import IGAmanip as manip
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrix_collocation as cplg_matrix_collocation
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from utils import rigid_motion

# modeleIGA = IGAparametrization(filename='plate_disk_collocation')
modeleIGA = IGAparametrization(filename='twoplatesDDcas1_collocation')

ti = time.time()

diametre = 0.0
epsi1 = 0.05
epsi2 = 0.05
u_inter = 15
v_inter = 15

count = 0
xi_inter = np.zeros((u_inter*v_inter,3))
xi_master = np.zeros((u_inter*v_inter,3))
xi_slave = np.zeros((u_inter*v_inter,3))

for i in range(u_inter):
    for j in range(v_inter):
        print(f'{count = }')
        print(f'{i = }')
        print(epsi1 + ((1-2*epsi1)/(u_inter-1))*i)
        xi_inter[count, 0] = epsi1 + ((1-2*epsi1)/(u_inter-1))*i
        print(f'{j = }')
        print(epsi2 + ((1-2*epsi2)/(v_inter-1))*j)
        xi_inter[count, 1] = epsi2 + ((1-2*epsi2)/(v_inter-1))*j
        count = count + 1
print(f'{xi_inter = }')
print(f'{xi_inter.shape[0] = }')

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

p = 1 #p-1
r = 3
# domains
nb_deg[:2,0] = p 
nb_deg[:2,1] = p
nb_deg[:2,2] = 0
modeleIGA.refine(nb_ref,nb_deg,additional_knots)

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}
modeleIGA.refine(nb_ref,nb_deg,additional_knots)

additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

nb_ref[:2, 0] = r 
nb_ref[:2, 1] = r
nb_ref[0, 2] = r - 1
nb_ref[1, 2] = r 

# p = 1 #p-1
# r = 3
# # domains
# nb_deg[:2,0] = p 
# nb_deg[:2,1] = 0
# nb_deg[:2,2] = 0
# modeleIGA.refine(nb_ref,nb_deg,additional_knots)

# nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
# nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

# additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}
# modeleIGA.refine(nb_ref,nb_deg,additional_knots)

# additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

# nb_ref[:2, 0] = r 
# nb_ref[:2, 1] = r
# nb_ref[0, 2] = r - 1
# nb_ref[1, 2] = r 

modeleIGA.refine(nb_ref,nb_deg,additional_knots)
# --
# STATIC STUDY

# MATRIX Assembly

ndof = modeleIGA._nb_dof_free
print(f'{ndof = }')
idof = modeleIGA._ind_dof_free[:ndof]-1


t1 = time.time()


# COUPLING MATRIX
# Cdata, Crow, Ccol = cplg_matrixU5( *modeleIGA.get_inputs4cplgmatrix_collocation(integrationOrder=6) )
# print(Cdata)

# Cdata1, Crow1, Ccol1, Cdata2, Crow2, Ccol2 = cplg_matrix_collocation( **modeleIGA.get_inputs4cplgmatrix_collocation(xi_master=np.array([0.9,0.333333,0]), 
#                                                                       xi_slave= np.array([0.166666,0.4,0]),
#                                                                       xi_inter= np.array([0.5,0.333333,0]),
#                                                                       i_master=1, i_slave=2, i_inter=3, d=0.1) )

# xi_master = np.array([[0.5,1./3.,0],[0.5,2./3.,0]]) 
# xi_slave = np.array([[0.5,1./3.,0],[0.5,2./3.,0]]) 
# xi_inter = np.array([[0.5,1./3.,0],[0.5,2./3.,0]]) 

# xi_master = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0],[2./3.,2./3.,0]]) 
# xi_slave = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0],[2./3.,2./3.,0]]) 
# xi_inter = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0],[2./3.,2./3.,0]])  


# xi_master = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0],[2./3.,2./3.,0],[0.5,0.5,0]]) 
# xi_slave = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0],[2./3.,2./3.,0],[0.5,0.5,0]]) 
# xi_inter = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0],[2./3.,2./3.,0],[0.5,0.5,0]])  

# xi_master = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0]]) 
# xi_slave = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0]]) 
# xi_inter = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0]]) 


# patchs 1/3 recouvrants


# xi_master = np.array([[0.55,1./3.,0],[0.75,1./3.,0],[0.95,1./3.,0],[0.55,0.5,0],[0.75,0.5,0],[0.95,0.5,0],[0.55,2./3.,0],[0.75,2./3.,0],[0.95,2./3.,0]]) 
# xi_slave = np.array([[0.05,1./3.,0],[0.25,1./3.,0],[0.45,1./3.,0],[0.05,0.5,0],[0.25,0.5,0],[0.45,0.5,0],[0.05,2./3.,0],[0.25,2./3.,0],[0.45,2./3.,0]])
# xi_inter = np.array([[0.36666,1./3.,0],[0.5,1./3.,0],[0.63333,1./3.,0],[0.36666,0.5,0],[0.5,0.5,0],[0.63333,0.5,0],[0.36666,2./3.,0],[0.5,2./3.,0],[0.63333,2./3.,0]])  

# xi_master = np.array([[0.55,0.1,0],[0.75,0.1,0],[0.95,0.1,0],[0.55,0.5,0],[0.75,0.5,0],[0.95,0.5,0],[0.55,0.9,0],[0.75,0.9,0],[0.95,0.9,0]]) 
# xi_slave = np.array([[0.05,0.1,0],[0.25,0.1,0],[0.45,0.1,0],[0.05,0.5,0],[0.25,0.5,0],[0.45,0.5,0],[0.05,0.9,0],[0.25,0.9,0],[0.45,0.9,0]])
# xi_inter = np.array([[0.36666,0.1,0],[0.5,0.1,0],[0.63333,0.1,0],[0.36666,0.5,0],[0.5,0.5,0],[0.63333,0.5,0],[0.36666,0.9,0],[0.5,0.9,0],[0.63333,0.9,0]])  

# xi_master = np.array([[0.55,0.1,0],[0.75,0.1,0],[0.95,0.1,0],[0.55,0.3,0],[0.75,0.3,0],[0.95,0.3,0],[0.55,0.7,0],[0.75,0.7,0],[0.95,0.7,0],[0.55,0.9,0],[0.75,0.9,0],[0.95,0.9,0]]) 
# xi_slave = np.array([[0.05,0.1,0],[0.25,0.1,0],[0.45,0.1,0],[0.05,0.3,0],[0.25,0.3,0],[0.45,0.3,0],[0.05,0.7,0],[0.25,0.7,0],[0.45,0.7,0],[0.05,0.9,0],[0.25,0.9,0],[0.45,0.9,0]])
# xi_inter = np.array([[0.36666,0.1,0],[0.5,0.1,0],[0.63333,0.1,0],[0.36666,0.3,0],[0.5,0.3,0],[0.63333,0.3,0],[0.36666,0.7,0],[0.5,0.7,0],[0.63333,0.7,0],[0.36666,0.9,0],[0.5,0.9,0],[0.63333,0.9,0]])  


# xi_master = np.array([[0.5,1./3.,0],[0.5,2./3.,0],[1,1./3.,0],[1,2./3.,0],[0.75,0.5,0]]) 
# xi_slave = np.array([[0,1./3.,0],[0,2./3.,0],[0.5,1./3.,0],[0.5,2./3.,0],[0.25,0.5,0]]) 
# xi_inter = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0],[2./3.,2./3.,0],[0.5,0.5,0]])  

# xi_master = np.array([[0.5,1./3.,0],[0.5,2./3.,0],[1,1./3.,0],[1,2./3.,0]]) 
# xi_slave = np.array([[0,1./3.,0],[0,2./3.,0],[0.5,1./3.,0],[0.5,2./3.,0]]) 
# xi_inter = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0],[2./3.,2./3.,0]])  

# xi_master = np.array([[0.5,1./3.,0],[0.5,2./3.,0],[1,1./3.,0]]) 
# xi_slave = np.array([[0,1./3.,0],[0,2./3.,0],[0.5,1./3.,0]]) 
# xi_inter = np.array([[1./3.,1./3.,0],[1./3.,2./3.,0],[2./3.,1./3.,0]])  

# xi_master = np.array([[0.75,1./3.,0],[0.75,2./3.,0]]) 
# xi_slave = np.array([[0.25,1./3.,0],[0.25,2./3.,0]]) 
# xi_inter = np.array([[0.5,1./3.,0],[0.5,2./3.,0]]) 



x_phys_m = np.zeros((xi_inter.shape[0],3))
Ctot = sp.coo_matrix(([],([],[])), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot),
                        dtype='float64').tocsc()

for ipoint in range(xi_inter.shape[0]):

    Cdata1, Crow1, Ccol1, Cdata2, Crow2, Ccol2, x_phys_m[ipoint,:] = cplg_matrix_collocation( **modeleIGA.get_inputs4cplgmatrix_collocation(xi_master= xi_master[ipoint,:], 
                                                                        xi_slave= xi_slave[ipoint,:],
                                                                        xi_inter= xi_inter[ipoint,:],
                                                                        i_master=1, i_slave=2, i_inter=3, d=diametre) )



    Cside1 = sp.coo_matrix((Cdata1,(Crow1,Ccol1)), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot),
                        dtype='float64').tocsc()

    Cside2 = sp.coo_matrix((Cdata2,(Crow2,Ccol2)), shape=(modeleIGA._nb_dof_tot,modeleIGA._nb_dof_tot),
                        dtype='float64').tocsc()

    # with open("/home/agagnaire/yeti/temp/Cside1.txt","w") as file:
    #     file.write(str(Cside1)) 

    # with open("/home/agagnaire/yeti/temp/Cside2.txt","w") as file:
    #     file.write(str(Cside2)) 

    # write txt

    Ctot += Cside1 + Cside1.transpose() + Cside2 + Cside2.transpose()

    # del Cdata1,Crow1,Ccol1,Cside1,Cdata2,Crow2,Ccol2,Cside2


print("Coupling matrix assembly done")

from pyevtk.hl import pointsToVTK

x = np.array(x_phys_m[:,0])
y = np.array(x_phys_m[:,1])
z = np.array(x_phys_m[:,2])

pointsToVTK("/home/agagnaire/yeti/examples/feti/results/colloc_points", x, y, z)

print(idof.shape)

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
C2solve = Ctot[idof,:][:,idof] * K2solve.max()

# plt.matshow(abs((Ctot[idof,:][:,idof]).todense()))
# plt.matshow(abs((Ctot[idof,:][:,idof]).todense()))
# plt.show()

# plt.spy(abs((Ctot[idof,:][:,idof]).todense()))
# plt.spy(abs(K2solve).todense())
# plt.show()

# print(null_space(K2solve.todense()+ C2solve.todense(),1e-2))

# print(Ctot[36:,:35].todense())
# valsing = svdvals(Ctot[36:,:35].todense())
# print(f'{valsing = }')

# singKC2solve = svdvals(K2solve.todense()+C2solve.todense())
# print(f'{singKC2solve = }')

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


idofl = modeleIGA._indCPbyPatch[2].size*2
print(idofl)
print('valeurs singulières')

singv = svdvals(Ctot[-idofl:,:-idofl].todense())
print(singv)

print(np.sum(singv < 1e-14))


# U, s, Vh = svd(C2solve.todense(),1e-3)

# print(s)

LU = sp.linalg.splu(K2solve + C2solve)
x  = LU.solve(Fb[idof])

# x = sp.linalg.spsolve(K2solve + C2solve, Fb[idof])

# print(f'{x[-18:] = }')

# for imode in range (kernel.shape[1]):
#     SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(kernel[:,imode]))
#     pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
#         f'mode_sol{imode}',SOL.transpose(),nb_ref=np.array([3,3,3]),
#         Flag=np.array([True,False,False])))

# print(kernel.shape[1])


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
    'coupling_mono_coloc',SOL.transpose(),nb_ref=np.array([3,3,3]),
    Flag=np.array([True,False,False])))

print(('Total time for Mono analysis : %.2f s' % (time.time() - ti)))

# print(f'{SOL = }')
# print(f'{np.max(SOL) = }')
# print(f'{u[-18:] = }')

# 0.30403442580078477
# u[-18:] = 0.0


# np.savetxt('/home/agagnaire/yeti/temp/Ctotd0.txt', Ctot.todense())
# np.savetxt("/home/agagnaire/yeti/temp/Ud0.txt",  u)
ud0 = np.loadtxt("/home/agagnaire/yeti/temp/Ud0.txt")

Ctotd0 = np.loadtxt("/home/agagnaire/yeti/temp/Ctotd0.txt")

# np.savetxt('/home/agagnaire/yeti/temp/Ktotd03.txt', Ktot.todense())


# print(f'{Ctot@u = }')

# print(f'{Ctot@ud0 = }')

# print(f'{np.max(Ctot - Ctotd0) = }')

# print(f'{np.max(Ctotd0@ud0) = }')

# print(f'{np.max(Ctotd0@u) = }')

# print(f'{np.max(Ctot@ud0) = }')

# print(f'{np.max(Ctot@u) = }')

# print(f'{np.max(u-ud0) = }')

# print(f'{(0.5*(u.T@Ktot@u) - u.T@Fb) = }')

# print(f'{(0.5*(ud0.T@Ktot@ud0) - ud0.T@Fb) = }')

# print(f'{(ud0.T@Ktot@ud0) = }')


# print(Ctot[-18:,:-18].shape)

# print(svdvals(Ctot[64:,:64].todense()))

# vecteuru = np.zeros(44)
# vecteuru[4] = 750.
# vecteuru[10] = 750.
# vecteuru[16] = 750.
# vecteuru[22] = 750.
# vecteuru[28] = 750.
# vecteuru[34] = 750.

# vecteuru[2] = 750./2.
# vecteuru[8] = 750./2.
# vecteuru[14] = 750./2.
# vecteuru[20] = 750./2.
# vecteuru[26] = 750./2.
# vecteuru[32] = 750./2.


# vecteuru = np.zeros(44)
# vecteuru[4] = 750.*(2./3.)
# vecteuru[10] = 750.*(2./3.)
# vecteuru[16] = 750.*(2./3.)

# vecteuru[22] = 750.
# vecteuru[28] = 750.
# vecteuru[34] = 750.

# vecteuru[2] = 750./3.
# vecteuru[8] = 750./3.
# vecteuru[14] = 750./3.

# vecteuru[18] = 750./3.
# vecteuru[24] = 750./3.
# vecteuru[30] = 750./3.

# vecteuru[20] = 750.*(2./3.)
# vecteuru[26] = 750.*(2./3.)
# vecteuru[32] = 750.*(2./3.)

# # vecteuru[-8:] = x[-8:]

# SOL2,u2 = rsol.reconstruction(**modeleIGA.get_inputs4solution(vecteuru[idof]))
# pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
#     'Solref',SOL2.transpose(),nb_ref=np.array([3,3,3]),
#     Flag=np.array([True,False,False])))


# Force  = (Ktot+Ctot*K2solve.max())@vecteuru
# Force2 = (Ktot+Ctot*K2solve.max())@u


# print(f'{Ktot[0,0] = }')
# np.savetxt("/home/agagnaire/yeti/temp/Ktot.txt",  Ktot.todense())

# print(f'{Ctot.max() = }')

# print(f'{vecteuru = }')

# print(f'{Fb = }')

# print(f'{Force = }')

# print(f'{Force2 = }')

# print(f'{Force2[0] = } {Force2[22] =}')
# print(f'{Force2[6] = } {Force2[28] =}')
# print(f'{Force2[12] = } {Force2[34] =}')


# # VTK output using Bezier elements
# for i_patch in range(modeleIGA._nb_patch):
#     pp.generate_vtu_bezier(**modeleIGA.get_inputs4postproc_bezier(
#         i_patch+1,
#         f'BentPipeBezier_2patchs_P{i_patch+1}',
#         SOL.transpose(),
#         ))





LU = sp.linalg.splu(K2solve + C2solve)
x  = LU.solve(Fb[idof])


if True:
    # resolution with FETI
    print('\n--\nResolution with FETI solver')

    # - dof infos
    ndof = modeleIGA._nb_dof_free
    idof = modeleIGA._ind_dof_free[:ndof]-1

    list_patch = np.where(np.isin(modeleIGA._ELT_TYPE,np.array(['U1','U3','U30', 'U10'])))[0] + 1
    list_curve = np.where(np.isin(modeleIGA._ELT_TYPE,np.array(['U00'])))[0] + 1
    list_lgrge = np.where(np.isin(modeleIGA._ELT_TYPE,np.array([ 'U4', 'U5', 'U23'])))[0] + 1

    idof_internal = np.array([],dtype=np.intp)
    mcrd = modeleIGA._mcrd
    for p in list_patch-1:
        idof_internal = np.concatenate(
            (idof_internal,np.repeat(np.unique(modeleIGA._IEN[p])-1,mcrd)*mcrd
             + np.tile(np.arange(0,mcrd),modeleIGA._indCPbyPatch[p].size)))
    idof_internal = np.intersect1d(idof_internal,idof)

    idof_lgrge = np.array([],dtype=np.intp)
    idof_lgrge_list = []
    for p in list_lgrge-1:
        temp = np.repeat(np.unique(modeleIGA._IEN[p])-1,mcrd)*mcrd \
               + np.tile(np.arange(0,mcrd),modeleIGA._indCPbyPatch[p].size)
        idof_lgrge = np.concatenate((idof_lgrge,temp))
        idof_lgrge_list.append(np.intersect1d(temp,idof))

    idof_lgrge_tot = np.intersect1d(idof_lgrge,idof)
    idof_lgrge     = idof_lgrge_list

    # - operators
    K2solve = Ktot[idof_internal,:][:,idof_internal]
    f2solve =   Fb[idof_internal]

    C2solve = Ctot[idof_lgrge_tot,:][:,idof_internal]

    # - factorization
    from solver import pseudoLUstep
    from solver import pseudoDense

    # print(Ktot[:8,:8].todense())
    np.savetxt("/home/agagnaire/yeti/temp/K1K2.txt",  Ktot[:16,:16].todense(), fmt='%.4e')
    np.savetxt("/home/agagnaire/yeti/temp/K22solve.txt", K2solve[4:12,4:12].todense(), fmt='%.4e')

    # LU = sp.linalg.splu(K2solve[4:12,4:12],permc_spec='MMD_AT_PLUS_A',
    #                             options=dict(SymmetricMode=True))
    # np.savetxt("/home/agagnaire/yeti/temp/LU2.txt", LU.U.todense(), fmt='%.4e')

    # np.savetxt("/home/agagnaire/yeti/temp/LU_U.txt", (LU.L@LU.U).todense(), fmt='%.4e')

    # print(f'{LU.perm_c = }')
    # print(f'{LU.perm_r = }')

    # LU = pseudoLU(K2solve,tol=1.e-8)
    

    LU = pseudoDense(K2solve,tol=1.e-8)

    moore = (pinv(K2solve.todense()))
    mo = sp.csc_matrix(moore)
    
    moore1 = (pinv(K2solve.todense()))

    null = null_space(K2solve.todense())
    

    np.savetxt("/home/agagnaire/yeti/temp/pinv1.txt",  moore1, fmt='%.4e')
    np.savetxt("/home/agagnaire/yeti/temp/null_space.txt",  null, fmt='%.4e')

  
    # LU1solve = pseudoLUstep(K2solve[:4,:4], tol=1.e-8)
    
    # LU = pseudoLUstep(K2solve,tol=1.e-8)
    # LU2solve = pseudoLUstep(K2solve[4:12,4:12], tol=1.e-8)


    # LU1 = pseudoLUstep(Ktot[:8,:8], tol=1.e-8)
    # LU2 = pseudoLUstep(Ktot[8:16,8:16], tol=1.e-8)

    def shur(x):
        y = np.zeros_like(x)
        y[:] = C2solve * LU.solve(C2solve.T * x)
        return y
    dofl = C2solve.shape[0]
    Sd = sp.linalg.LinearOperator((dofl,dofl),matvec=shur)
    t  = C2solve*LU.solve(f2solve)

    # - projector
    R = LU.R
    G =-C2solve*R
    e =-R.T * f2solve

    import scipy.linalg as la
    invGG= la.inv( (G.T * G).toarray() )
    print(la.eig( (G.T * G).toarray() ))
    def project(v):
        y = np.zeros_like(v)
        y[:] = v[:] - G * invGG.dot( G.T*v )
        return y
    Proj = sp.linalg.LinearOperator((dofl,dofl),matvec=project)

    lmbda0 = G.dot(invGG.dot(e))

    # - preconditionner

    # - resolution
    from solver import PCPGortho
    t1 = time.time()
    lmbda,nbiter= PCPGortho(Sd,t, x0=lmbda0, P=Proj,
                            M=None, tol=1.e-8, maxiter=1000,savetxt=True)
    print(' Resolution with PCPG algo')
    print(' (duration : %.2f s, nb iter : %i).\n\n' % (time.time() - t1,nbiter))


    # - postproc
    alpha = invGG.dot( G.T * (t - Sd(lmbda) ) )
    utot = np.zeros(modeleIGA._nb_dof_tot)
    utot[idof_internal] = LU.solve(f2solve - C2solve.T*lmbda) + R*alpha

    SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(utot[idof]))
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'coupling_pcpg_colloc',SOL.transpose(),nb_ref=np.array([3,3,3]),
        Flag=np.array([True,True,True])))
    
exit()
n_sample = 1000

# Attention : enlever le file name 'patch_1'
x_sample_P1, u_sample_P1, dudx_sample_P1, norm_sample_P1, tan_sample_P1, dudxi_sample_P1 = \
    pp.postproc_curve_2d(**modeleIGA.get_inputs4post_curve_2D(1, 2, n_sample, SOL.transpose())) 
x_sample_P2, u_sample_P2, dudx_sample_P2, norm_sample_P2, tan_sample_P2, dudxi_sample_P2 = \
    pp.postproc_curve_2d(**modeleIGA.get_inputs4post_curve_2D(2, 1, n_sample, SOL.transpose()))

# print(x_sample_P1)

delta = np.zeros((2, n_sample))
for i_sample in range(n_sample):
    for i in range(2):
        delta[i, i_sample] = (dudx_sample_P1[i, :, i_sample] @ norm_sample_P1[:, i_sample]) - \
                                (dudx_sample_P2[i, :, i_sample] @ norm_sample_P2[:, i_sample])


print(np.shape(u_sample_P1))



deltaT = np.zeros((2, n_sample))
for i_sample in range(n_sample):
    for i in range(2):
        deltaT[i, i_sample] = (dudx_sample_P1[i, :, i_sample] @ tan_sample_P1[:, i_sample]) - \
                                (dudx_sample_P2[i, :, i_sample] @ tan_sample_P2[:, i_sample])
    
deltau = np.zeros((2, n_sample))
for i_sample in range(n_sample):
    for i in range(2):
        deltau[i, i_sample] = u_sample_P1[i, i_sample] - u_sample_P2[i, i_sample]


dudx_t1 = np.zeros((2, n_sample))
for i_sample in range(n_sample):
    dudx_t1[0, i_sample] = dudx_sample_P1[0, 0, i_sample] * tan_sample_P1[0, i_sample] + dudx_sample_P1[0, 1, i_sample] * tan_sample_P1[1, i_sample]
    dudx_t1[1, i_sample] = dudx_sample_P1[1, 0, i_sample] * tan_sample_P1[0, i_sample] + dudx_sample_P1[1, 1, i_sample] * tan_sample_P1[1, i_sample]

dudx_t2 = np.zeros((2, n_sample))
for i_sample in range(n_sample):
    for i in range(2):
        dudx_t2[i, i_sample] = (dudx_sample_P2[i, :, i_sample] @ tan_sample_P2[:, i_sample]) 

dudx_n1 = np.zeros((2, n_sample))
for i_sample in range(n_sample):
    for i in range(2):
        dudx_n1[i, i_sample] = (dudx_sample_P1[i, :, i_sample] @ norm_sample_P1[:, i_sample]) 

dudx_n2 = np.zeros((2, n_sample))
for i_sample in range(n_sample):
    for i in range(2):
        dudx_n2[i, i_sample] = (dudx_sample_P2[i, :, i_sample] @ norm_sample_P2[:, i_sample]) 
                                    



# print(delta[0, :]) 


# print(dudx_sample_P1[0,0])
# print(dudx_n1[0,:])




# cpi1 = manip.get_boundCPindice_wEdges(modeleIGA._Nkv,modeleIGA._Jpqr,modeleIGA._dim, 2, num_patch=0, offset=0,num_orientation=0)
# print(cpi1)

# print(np.max(cpi1))

# cpi2 = manip.get_boundCPindice(modeleIGA._Nkv,modeleIGA._Jpqr, 1, num_patch=1, offset=0) + np.max(cpi1) + 1
# print(cpi2)

# SOLi1 = SOL[cpi1]

# SOLi2 = SOL[cpi2]

# print(SOL[cpi1])
# print(SOL[cpi2])
# print(modeleIGA._dim)


# pp.generatecplginterfacetxt(*modeleIGA.get_inputs4postprocCPLG(
#     'coupling_txt',SOL.transpose(), nb_ref=5,
#     Flag=np.array([True, False, False])))


np.save('/home/agagnaire/yeti/temp/x_sample_P1', x_sample_P1)
np.save('/home/agagnaire/yeti/temp/x_sample_P2', x_sample_P2)
np.save('/home/agagnaire/yeti/temp/u_sample_P1', u_sample_P1)
np.save('/home/agagnaire/yeti/temp/u_sample_P2', u_sample_P2)
np.save('/home/agagnaire/yeti/temp/dudx_sample_P1', dudx_sample_P1)
np.save('/home/agagnaire/yeti/temp/dudx_sample_P2', dudx_sample_P2)

np.savetxt('/home/agagnaire/yeti/temp/delta_u', np.abs(u_sample_P1-u_sample_P2))

print(np.max(delta))
print(np.max(delta)/np.max(dudx_sample_P1))
print(np.max(delta)/np.max(dudx_sample_P2))

print(f'{np.max(np.abs(u_sample_P1-u_sample_P2)) = }')


print(np.array([[1,4],[7,6]])@np.array([3,5]))


plt.subplot(221)
plt.plot(range(n_sample), dudx_sample_P1[0,0], label='dUx/dx sous domaine 1') # dUx/dx
plt.plot(range(n_sample), dudx_sample_P2[0,0], label='dUx/dx sous domaine 2')
plt.ylim(-1.5e-6,1.5e-6)
plt.legend()
plt.title('dUx/dx')

plt.subplot(222)
plt.plot(range(n_sample), dudx_sample_P1[0,1], label='dUx/dy sous domaine 1') # dUx/dy
plt.plot(range(n_sample), dudx_sample_P2[0,1], label='dUx/dy sous domaine 2')
plt.ylim(-1.25e-6,1.25e-6)
plt.legend()
plt.title('dUx/dy')

plt.subplot(223)
plt.plot(range(n_sample), dudx_sample_P1[1,0], label='dUy/dx sous domaine 1') # dUy/dx
plt.plot(range(n_sample), dudx_sample_P2[1,0], label='dUy/dx sous domaine 2')
plt.ylim(-1.5e-6,1.5e-6)
plt.legend()
plt.title('dUy/dx')

plt.subplot(224)
plt.plot(range(n_sample), dudx_sample_P1[1,1], label='dUy/dy sous domaine 1') # dUy/dy
plt.plot(range(n_sample), dudx_sample_P2[1,1], label='dUy/dy sous domaine 2')
plt.ylim(-1.25e-6,1.25e-6)
plt.legend()
plt.title('dUy/dy')
plt.show()

plt.subplot(221)
plt.plot(range(n_sample), dudxi_sample_P1[0,0], label='dUx/dXi sous domaine 1')
plt.plot(range(n_sample), dudxi_sample_P2[0,0], label='dUx/dXi sous domaine 2')
plt.autoscale()
plt.legend()
plt.title('dUx/dXi')

plt.subplot(222)
plt.plot(range(n_sample), dudxi_sample_P1[0,1], label='dUx/dEta sous domaine 1')
plt.plot(range(n_sample), dudxi_sample_P2[0,1], label='dUx/dEta sous domaine 2')
plt.autoscale()
plt.legend()
plt.title('dUx/dEta')

plt.subplot(223)
plt.plot(range(n_sample), dudxi_sample_P1[1,0], label='dUy/dXi sous domaine 1')
plt.plot(range(n_sample), dudxi_sample_P2[1,0], label='dUy/dXi sous domaine 2')
plt.autoscale()
plt.legend()
plt.title('dUy/dXi')

plt.subplot(224)
plt.plot(range(n_sample), dudxi_sample_P1[1,1], label='dUy/dEta sous domaine 1')
plt.plot(range(n_sample), dudxi_sample_P2[1,1], label='dUy/dEta sous domaine 2')
plt.autoscale()
plt.legend()
plt.title('dUy/dEta')
plt.show()



plt.subplot(221)
plt.plot(range(n_sample), dudx_n1[0,:], label='dUx/dx*nx + dUx/dy*ny sous domaine 1') # dUx/dx
plt.plot(range(n_sample), dudx_n2[0,:], label='dUx/dx*nx + dUx/dy*ny sous domaine 2')
plt.autoscale()
plt.legend()
plt.title('dUx/dx*nx + dUx/dy*ny')



plt.subplot(222)
plt.plot(range(n_sample), dudx_n1[1,:], label='dUy/dx*nx + dUy/dy*ny sous domaine 1') # dUx/dy
plt.plot(range(n_sample), dudx_n2[1,:], label='dUy/dx*nx + dUy/dy*ny sous domaine 2')
plt.autoscale()
plt.legend()
plt.title('dUy/dx*nx + dUy/dy*ny')


# 0 dUx 1 dUy, 0 dx, 1 dy     Ex :  dudx_sample_P1[1,1] dUy/dy interface 1


plt.subplot(223)
plt.plot(range(n_sample), dudx_t1[0,:], label='dUx/dx*tx + dUx/dy*ty sous domaine 1') # dUy/dx
plt.plot(range(n_sample), dudx_t2[0,:], label='dUx/dx*tx + dUx/dy*ty sous domaine 2')
plt.autoscale()
plt.legend()
plt.title('dUx/dx*tx + dUx/dy*ty ')


plt.subplot(224)
plt.plot(range(n_sample), dudx_t1[1,:], label='dUy/dx*tx + dUy/dy*ty sous domaine 1') # dUy/dy
plt.plot(range(n_sample), dudx_t2[1,:], label='dUy/dx*tx + dUy/dy*ty sous domaine 2')
plt.autoscale()
plt.legend()
plt.title('dUy/dx*tx + dUy/dy*ty')
plt.show()

delta_t = np.abs(dudx_t1-dudx_t2)
delta_n = np.abs(dudx_t1-dudx_t2)


plt.subplot(221)
plt.plot(range(n_sample), delta_n[0,:]/np.max((np.abs(dudx_n1[0,:]),np.abs(dudx_n2[0,:])))) 
plt.autoscale()
plt.title('delta norm x')


plt.subplot(222)
plt.plot(range(n_sample), delta_n[1,:]/np.max((np.abs(dudx_n1[1,:]),np.abs(dudx_n2[1,:])))) 
plt.autoscale()
plt.title('delta norm y')

plt.subplot(223)
plt.plot(range(n_sample), delta_t[0,:]/np.max((np.abs(dudx_t1[0,:]),np.abs(dudx_t2[0,:])))) 
plt.autoscale()
plt.title('delta tan x')


plt.subplot(224)
plt.plot(range(n_sample), delta_t[1,:]/np.max((np.abs(dudx_t1[1,:]),np.abs(dudx_t2[1,:])))) 
plt.autoscale()
plt.title('delta tan y')
plt.show()

exit()

# plt.figure()

# plt.subplot(221)
# plt.plot(range(n_sample), dudx_sample_P1[0,0], label='sous domaine 1') # dUx/dx
# plt.plot(range(n_sample), dudx_sample_P2[0,0], label='sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUx/dx')



# plt.subplot(222)
# plt.plot(range(n_sample), dudx_sample_P1[0,1], label='sous domaine 1') # dUx/dy
# plt.plot(range(n_sample), dudx_sample_P2[0,1], label='sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUx/dy')


# # 0 dUx 1 dUy, 0 dx, 1 dy     Ex :  dudx_sample_P1[1,1] dUy/dy interface 1


# plt.subplot(223)
# plt.plot(range(n_sample), dudx_sample_P1[1,0], label='sous domaine 1') # dUy/dx
# plt.plot(range(n_sample), dudx_sample_P2[1,0], label='sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUy/dx')


# plt.subplot(224)
# plt.plot(range(n_sample), dudx_sample_P1[1,1], label='sous domaine 1') # dUy/dy
# plt.plot(range(n_sample), dudx_sample_P2[1,1], label='sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUy/dy')
# plt.tight_layout()
# plt.savefig('dUdX.png')
# # plt.show()
# plt.clf()

# plt.subplot(221)
# plt.plot(range(n_sample), dudx_n1[0,:], label='sous domaine 1') 
# plt.plot(range(n_sample), dudx_n2[0,:], label='sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUx/dx*nx + dUx/dy*ny')


# plt.subplot(222)
# plt.plot(range(n_sample), dudx_n1[1,:], label='sous domaine 1') 
# plt.plot(range(n_sample), dudx_n2[1,:], label='sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUy/dx*nx + dUy/dy*ny')


# plt.subplot(223)
# plt.plot(range(n_sample), dudx_t1[0,:], label='sous domaine 1') 
# plt.plot(range(n_sample), dudx_t2[0,:], label='sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUx/dx*tx + dUx/dy*ty ')


# plt.subplot(224)
# plt.plot(range(n_sample), dudx_t1[1,:], label='sous domaine 1') 
# plt.plot(range(n_sample), dudx_t2[1,:], label='sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUy/dx*tx + dUy/dy*ty')
# plt.tight_layout()
# plt.savefig('dUdX_n_t.png')
# # plt.show()
# plt.clf()

# delta_t = np.abs(dudx_t1-dudx_t2)
# delta_n = np.abs(dudx_t1-dudx_t2)
# print(np.max((np.abs(dudx_n1),np.abs(dudx_n2))))
# print(np.max((np.abs(dudx_t1),np.abs(dudx_t2))))
# print(np.max((delta_n)))
# print(np.max((delta_t)))

# plt.subplot(221)
# plt.plot(range(n_sample), delta_n[0,:]/np.max((np.abs(dudx_n1[0,:]),np.abs(dudx_n2[0,:])))) 
# plt.autoscale()
# plt.title('delta norm x')


# plt.subplot(222)
# plt.plot(range(n_sample), delta_n[1,:]/np.max((np.abs(dudx_n1[1,:]),np.abs(dudx_n2[1,:])))) 
# plt.autoscale()
# plt.title('delta norm y')

# plt.subplot(223)
# plt.plot(range(n_sample), delta_t[0,:]/np.max((np.abs(dudx_t1[0,:]),np.abs(dudx_t2[0,:])))) 
# plt.autoscale()
# plt.title('delta tan x')


# plt.subplot(224)
# plt.plot(range(n_sample), delta_t[1,:]/np.max((np.abs(dudx_t1[1,:]),np.abs(dudx_t2[1,:])))) 
# plt.autoscale()
# plt.title('delta tan y')
# plt.tight_layout()
# plt.savefig('delta_n_t.png')
# # plt.show()
# plt.clf()




# plt.subplot(221)
# plt.plot(range(n_sample), dudx_sample_P1[0,0], label='dUx/dx sous domaine 1') # dUx/dx
# plt.plot(range(n_sample), dudx_sample_P2[0,0], label='dUx/dx sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUx/dx')



# plt.subplot(222)
# plt.plot(range(n_sample), dudx_sample_P1[0,1], label='dUx/dy sous domaine 1') # dUx/dy
# plt.plot(range(n_sample), dudx_sample_P2[0,1], label='dUx/dy sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUx/dy')


# # 0 dUx 1 dUy, 0 dx, 1 dy     Ex :  dudx_sample_P1[1,1] dUy/dy interface 1


# plt.subplot(223)
# plt.plot(range(n_sample), dudx_sample_P1[1,0], label='dUy/dx sous domaine 1') # dUy/dx
# plt.plot(range(n_sample), dudx_sample_P2[1,0], label='dUy/dx sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.title('dUy/dx')


# plt.subplot(224)
# plt.plot(range(n_sample), dudx_sample_P1[1,1], label='dUy/dy sous domaine 1') # dUy/dy
# plt.plot(range(n_sample), dudx_sample_P2[1,1], label='dUy/dy sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.show()
# plt.title('dUy/dy')

# plt.subplot(221)
# plt.plot(range(n_sample), dudxi_sample_P1[0,0], label='dUx/dXi sous domaine 1')
# plt.plot(range(n_sample), dudxi_sample_P2[0,0], label='dUx/dXi sous domaine 2')
# plt.autoscale()
# plt.legend()

# plt.subplot(222)
# plt.plot(range(n_sample), dudxi_sample_P1[0,1], label='dUx/dEta sous domaine 1')
# plt.plot(range(n_sample), dudxi_sample_P2[0,1], label='dUx/dEta sous domaine 2')
# plt.autoscale()
# plt.legend()

# plt.subplot(223)
# plt.plot(range(n_sample), dudxi_sample_P1[1,0], label='dUy/dXi sous domaine 1')
# plt.plot(range(n_sample), dudxi_sample_P2[1,0], label='dUy/dXi sous domaine 2')
# plt.autoscale()
# plt.legend()

# plt.subplot(224)
# plt.plot(range(n_sample), dudxi_sample_P1[1,1], label='dUy/dEta sous domaine 1')
# plt.plot(range(n_sample), dudxi_sample_P2[1,1], label='dUy/dEta sous domaine 2')
# plt.autoscale()
# plt.legend()
# plt.show()

# exit()

# RESOLUTION Decomposition de Domaine
time2 = time.time()
ip1 = np.array([0])
ip2 = np.array([1])
ipl = np.array([7,8])

dof1 = np.array([],dtype=np.intp)
for p in ip1:
    dof1 = np.concatenate((dof1,np.repeat(np.unique(modeleIGA._IEN[p])-1,2)*2
                           + np.tile(np.array([0,1]),modeleIGA._indCPbyPatch[p].size)))
dof1 = np.intersect1d(dof1,idof)

K1  = Ktot[dof1,:][:,dof1]
F1  = Fb[dof1]
LU1 = sp.linalg.splu(K1)


dof2 = np.array([],dtype=np.intp)
for p in ip2:
    dof2 = np.concatenate((dof2,np.repeat(np.unique(modeleIGA._IEN[p])-1,2)*2
                           + np.tile(np.array([0,1]),modeleIGA._indCPbyPatch[p].size)))

# Inutile car le ss domaine 2 n'a pas de CL de Dirichlet
dof2 = np.intersect1d(dof2,idof)

K2   = Ktot[dof2,:][:,dof2]

sp.save_npz('/home/agagnaire/yeti/temp/K2', K2)

exit()

# Indices des points de controle decrivant les coins
iv2 = manip.get_vertexCPindice(modeleIGA._Nkv,modeleIGA._Jpqr,modeleIGA._dim,num_patch=ip2[0])
print("iv2 : ", iv2)
print("ip2 : ", ip2)
print(modeleIGA._indCPbyPatch[ip2[0]])
print(modeleIGA._indCPbyPatch[ip2[0]][iv2])

# TEST

# FIN TEST

# ddlr= np.array([(modeleIGA._indCPbyPatch[ip2[0]][iv2[0]]-1)*3+0,
#                 (modeleIGA._indCPbyPatch[ip2[0]][iv2[0]]-1)*3+1,
#                 (modeleIGA._indCPbyPatch[ip2[0]][iv2[0]]-1)*3+2,
#                 (modeleIGA._indCPbyPatch[ip2[0]][iv2[1]]-1)*3+0,
#                 (modeleIGA._indCPbyPatch[ip2[0]][iv2[1]]-1)*3+2,
#                 (modeleIGA._indCPbyPatch[ip2[0]][iv2[1]]  )*3+0])

# EN 2D, à tester
ddlr= np.array([(modeleIGA._indCPbyPatch[ip2[0]][iv2[0]]-1)*2+0,
                (modeleIGA._indCPbyPatch[ip2[0]][iv2[0]]-1)*2+1,
                (modeleIGA._indCPbyPatch[ip2[0]][iv2[1]]-1)*2+0])


print(ddlr)


dof2pp = np.setxor1d(dof2,ddlr)
print(dof2pp)
K2pp = Ktot[dof2pp,:][:,dof2pp]
K2pr = Ktot[dof2pp,:][:,ddlr]
print(K2pp.shape)
print(K2pr.shape)

exit()
ddlr_loc = np.zeros(6,dtype=np.intp)
for i in np.arange(0,6):
    ddlr_loc[i] = np.where(dof2 == ddlr[i])[0][0]
dof2pp_loc = np.setxor1d(np.arange(0,dof2.size),ddlr_loc)

LU2 = sp.linalg.splu(K2pp)
R2  = np.zeros((dof2.size,6))
for r in np.arange(0,6):
    R2[dof2pp_loc,r] = -LU2.solve(K2pr[:,r].toarray().flatten())
R2[ddlr_loc,:] = np.identity(6)

F2  = Fb[dof2]


dofl = np.array([],dtype=np.intp)
for p in ipl:
    dofl = np.concatenate((dofl,np.repeat(np.unique(modeleIGA._IEN[p])-1,3)*3
                           + np.tile(np.array([0,1,2]),modeleIGA._indCPbyPatch[p].size)))
dofl = np.intersect1d(dofl,idof)

C1 = Ctot[dofl,:][:,dof1]
C2 = Ctot[dofl,:][:,dof2]


def K2star(f2):
    v2 = np.zeros(dof2.size)
    v2[dof2pp_loc] = LU2.solve(f2[dof2pp_loc])
    return v2

def shurOp_wLU(di):
    ld = di[:-6]
    al = di[-6:]
    y1 = LU1.solve(C1.transpose().dot(ld))
    y2 = K2star(C2.transpose().dot(ld))

    ytot = np.zeros_like(di)
    ytot[:-6] = C1.dot(y1) + C2.dot(y2) -C2.dot(R2.dot(al))
    ytot[-6:] = -(C2.dot(R2).transpose()).dot(ld)
    return ytot
Sd = sp.linalg.LinearOperator((dofl.size+6,dofl.size+6),matvec=shurOp_wLU, dtype='f')

# interface problem
t1 = LU1.solve(F1)
t2 = K2star(F2)
t = np.zeros(dofl.size+6)
t[:-6] = C1.dot(t1) + C2.dot(t2)
t[-6:] = -R2.transpose().dot(F2)
lmbda,info = sp.linalg.cg(Sd,t,maxiter=50000)

# local solutions
u1 = LU1.solve(F1-C1.transpose().dot(lmbda[:-6]))
u2 = K2star(F2-C2.transpose().dot(lmbda[:-6])) + R2.dot(lmbda[-6:])

xDD = np.concatenate((u1,np.concatenate((u2,lmbda[:-6]))))

print('\n Time for DD solving : %.2f s\n\n' % (time.time() - time2))



# Postprocessing
SOL,u = rsol.reconstruction(*modeleIGA.get_inputs4solution(xDD))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'coupling',SOL.transpose(),nb_ref=3*np.array([1,1,1]),
    Flag=np.array([True,True,False])))