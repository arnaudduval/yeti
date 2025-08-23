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


def collocation(ref, ref_l, deg, deg_l, d, u_inter, v_inter):

    modeleIGA = IGAparametrization(filename='twoplatesDDcas1_collocation')

    ti = time.time()

    p = deg #p-1
    r = ref
    diametre = d
    epsi1 = 0.005
    epsi2 = 0.005


    count = 0
    xi_inter = np.zeros((u_inter*v_inter,3))
    xi_master = np.zeros((u_inter*v_inter,3))
    xi_slave = np.zeros((u_inter*v_inter,3))

    for i in range(u_inter):
        for j in range(v_inter):
            xi_inter[count, 0] = epsi1 + ((1-2*epsi1)/(u_inter-1))*i
            xi_inter[count, 1] = epsi2 + ((1-2*epsi2)/(v_inter-1))*j
            count = count + 1
    print(f'{xi_inter = }')
    print(f'{xi_inter.shape[0] = }')

    nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
    nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
    additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}


    # domains
    nb_deg[:2,:2] = p 
    nb_deg[:2,2] = deg_l
    modeleIGA.refine(nb_ref,nb_deg,additional_knots)

    nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
    nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

    additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}
    modeleIGA.refine(nb_ref,nb_deg,additional_knots)

    additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

    nb_ref[:2, :2] = r 
    nb_ref[0, 2] = ref_l - 1
    nb_ref[1, 2] = ref_l 


    modeleIGA.refine(nb_ref,nb_deg,additional_knots)
    # --
    # STATIC STUDY

    # MATRIX Assembly

    ndof = modeleIGA._nb_dof_free
    print(f'{ndof = }')
    idof = modeleIGA._ind_dof_free[:ndof]-1


    t1 = time.time()


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

    print(x_phys_m)

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
    # print(f'{Ktot.shape = }')
    # print(f'{modeleIGA._nnode = }')
    # print(f'{modeleIGA._nb_dof_tot = }')

    print('valeurs singulières')

    idofl = modeleIGA._indCPbyPatch[2].size*2
    print('valeurs singulières')

    singv = svdvals(Ctot[-idofl:,:-idofl].todense())
    print(singv)

    print(np.sum(singv < 1e-14))

    # Monolithic SOLVE
    t2 = time.time()
    K2solve = Ktot[idof,:][:,idof]
    C2solve = Ctot[idof,:][:,idof] * K2solve.max()


    # LU = sp.linalg.splu(K2solve + C2solve)
    # x  = LU.solve(Fb[idof])

    x = sp.linalg.spsolve(K2solve + C2solve, Fb[idof])


    print(('\n Time for monolithique solving : %.2f s\n\n' % (time.time() - t2)))

    # Postprocessing
    # print(modeleIGA.get_inputs4solution(x))
  

    SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        f'collocation_test/coloc_mono_{ref=}_{ref_l=}_{deg=}_{deg_l=}_{d=}_{u_inter=}_{v_inter=}',SOL.transpose(),nb_ref=np.array([3,3,3]),
        Flag=np.array([True,False,False])))

    print(('Total time for Mono analysis : %.2f s' % (time.time() - ti)))

    # print(f'{SOL = }')
   

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

        np.savetxt("/home/agagnaire/yeti/temp/K1K2.txt",  Ktot[:16,:16].todense(), fmt='%.4e')
        np.savetxt("/home/agagnaire/yeti/temp/K22solve.txt", K2solve[4:12,4:12].todense(), fmt='%.4e')
        

        LU = pseudoDense(K2solve,tol=1.e-8)

        moore = (pinv(K2solve.todense()))
        mo = sp.csc_matrix(moore)
        
        moore1 = (pinv(K2solve.todense()))

        null = null_space(K2solve.todense())
        

        np.savetxt("/home/agagnaire/yeti/temp/pinv1.txt",  moore1, fmt='%.4e')
        np.savetxt("/home/agagnaire/yeti/temp/null_space.txt",  null, fmt='%.4e')


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
            f'collocation_test/coloc_pcpg_{ref=}_{ref_l=}_{deg=}_{deg_l=}_{d=}_{u_inter=}_{v_inter=}',SOL.transpose(),nb_ref=np.array([3,3,3]),
            Flag=np.array([True,True,True])))
    return singv

a = np.empty((6,6), dtype=object)
print(a)

with open("/home/agagnaire/yeti/examples/feti/results/Error.txt", "w") as text_file:
                print(f'ref = 3, ref_l = 3, deg = 1, deg_l = 0, d=0.0', file=text_file) 

for i in range(0,6):
    for j in range(0,6):
        try:
            singv = collocation(ref=i, ref_l=j, deg= 1, deg_l= 0, d = 0.1, u_inter=30, v_inter=30)
            if (np.sum(singv < 1e-10)) > 0:
                 a[i,j] = np.sum(singv < 1e-10)
                 with open("/home/agagnaire/yeti/examples/feti/results/Error.txt", "a") as text_file:
                    print(f'singv null = {np.sum(singv < 1e-10)} #colloc_points = {i*j}', file=text_file)
                    print(singv, file=text_file) 
            else:
                 a[i,j] = 'OK'
                
        except:
            a[i,j] = 'error'
            with open("/home/agagnaire/yeti/examples/feti/results/Error.txt", "a") as text_file:
                print(f'error #colloc_points = {i*j}', file=text_file)    
            pass

print(a)
print(a[3,4])
# for i in range(0,5):
#      print(i)  
# a.tofile('foo.csv',sep=',',format='object')
np.savetxt("deg = 1, deg_l = 0, d=0.1, u_inter=30, v_inter=30, petitepsi.csv", a, delimiter=",", fmt='%s')

# with open("/home/agagnaire/yeti/examples/feti/results/Error2.txt", "w") as text_file:
#                 print(f'ref = 2, deg = 1, d = 0.1, u_inter=5, v_inter=5', file=text_file) 

# for i in range(0,4):
#     for j in range(0,3):
#         try:
#             collocation(ref=3, ref_l=i, deg= 1, deg_l= j, d = 0.1, u_inter=5, v_inter=5)
#         except:
#             with open("/home/agagnaire/yeti/examples/feti/results/Error2.txt", "a") as text_file:
#                 print(f'error refl = {i}  degl = {j}', file=text_file)    
#             pass



exit()
