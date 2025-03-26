# Python module
import numpy as np
import scipy.sparse as sp
import sys
import time

#IGA module
from yeti_iga.solver import PCPGortho
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.preprocessing.igaparametrization import IGAmanip as manip
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp

modeleIGA = IGAparametrization(filename='twoplatesDD')

ti = time.time()

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

p = 3
r = 3
# domains
nb_deg[:2,:2] = p
nb_ref[:2, 0] = r #r+1
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



if True:
    # resolution with FETI
    print('\n--\nResolution with FETI solver')

    # - dof infos
    ndof = modeleIGA._nb_dof_free
    idof = modeleIGA._ind_dof_free[:ndof]-1

    list_patch = np.where(np.isin(modeleIGA._ELT_TYPE,np.array(['U1','U3','U30'])))[0] + 1
    list_curve = np.where(np.isin(modeleIGA._ELT_TYPE,np.array(['U00'])))[0] + 1
    list_lgrge = np.where(np.isin(modeleIGA._ELT_TYPE,np.array([ 'U4'])))[0] + 1

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
    LU = pseudoLUstep(K2solve,tol=1.e-8)

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
                            M=None, tol=1.e-10, maxiter=300,savetxt=True)
    print(' Resolution with PCPG algo')
    print(' (duration : %.2f s, nb iter : %i).\n\n' % (time.time() - t1,nbiter))


    # - postproc
    alpha = invGG.dot( G.T * (t - Sd(lmbda) ) )
    utot = np.zeros(modeleIGA._nb_dof_tot)
    utot[idof_internal] = LU.solve(f2solve - C2solve.T*lmbda) + R*alpha

    SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(utot[idof]))
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'coupling_pcpg',SOL.transpose(),nb_ref=np.array([3,3,1]),
        Flag=np.array([True,True,True])))


    n_sample = 1000

    x_sample_P1, u_sample_P1, dudx_sample_P1, norm_sample_P1, tan_sample_P1, dudxi_sample_P1 = \
        pp.postproc_curve_2d(**modeleIGA.get_inputs4post_curve_2D(1, 2, n_sample, 'patch_1', SOL.transpose()))
    x_sample_P2, u_sample_P2, dudx_sample_P2, norm_sample_P2, tan_sample_P2, dudxi_sample_P2 = \
        pp.postproc_curve_2d(**modeleIGA.get_inputs4post_curve_2D(2, 1, n_sample, 'patch_2', SOL.transpose()))

    # print(x_sample_P1)

    delta = np.zeros((2, n_sample))
    for i_sample in range(n_sample):
        for i in range(2):
            delta[i, i_sample] = (dudx_sample_P1[i, :, i_sample] @ tan_sample_P1[:, i_sample]) - \
                                 (dudx_sample_P2[i, :, i_sample] @ tan_sample_P2[:, i_sample])


    print(np.shape(u_sample_P1))


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


    import matplotlib.pyplot as plt

    # print(delta[0, :])



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


    print(np.max(delta))
    print(np.max(dudx_sample_P2[0,1]))

    print(np.max(delta/dudx_sample_P1))

    # plt.plot(range(n_sample), delta[0, :], label='comp 1')
    # plt.plot(range(n_sample), delta[1, :], label='comp 2')
    # plt.autoscale()
    # plt.legend()
    # plt.show()

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