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
This case is described in the following publication :
Hirschler, T., R. Bouclier, A. Duval, T. Elguedj, et J. Morlier.
The Embedded Isogeometric Kirchhoff–Love Shell: From Design to Shape Optimization of Non-Conforming Stiffened Multipatch Structures.
Computer Methods in Applied Mechanics and Engineering 349 (juin 2019): 774‑97.
https://doi.org/10.1016/j.cma.2019.02.042.

The structure is made of 2 coupled Kirchhoff-Love
Result is compared with reference

"""

# Python module
import numpy as np
import scipy.sparse as sp
import sys
import time

#IGA module
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip

EXEMPLE_NB = 5

# Creation of the IGA object
# --------------------------
modeleIGA = IGAparametrization(filename='Tbeam2cplg')

ti = time.time()

# SHAPE Modification
nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

# REFINEMENT for Analysis
nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

# Data read from input file are 
# - two shells of degree 1
# - 4 interface elements of degree 1
# - 2 Lagrange fields of degree 0

p = 2
nb_deg[:2,:2] = p

r = 3

nb_ref[:,0]  = np.array([3+r,2+r,0])
additional_knots["patches"] = np.array([1])
additional_knots["1"] = np.linspace(0,1,10*2**r)[1:-1]
additional_knots["2"] = np.linspace(0,1,4*2**r)[1:-1]
    
nb_ref[0,2:6]= 5+r

nb_deg[0,6:] = np.array([p,p-1])
nb_ref[0,6:] = 3+r

# Final degree for shells is 3
# Final degree for interfaces is 1
# Finale degree for Lagrange fields is 2 (displacements) and 1 (rotations)

modeleIGA.refine(nb_ref,nb_deg,additional_knots)

modeleIGA._NBPINT[ np.where(modeleIGA._ELT_TYPE == 'U00') ] = 6**modeleIGA._dim.min()

# Add boundary conditions : lock displacement for the 2nd row of CPs
icps0 = manip.get_boundCPindice_wEdges(modeleIGA._Nkv,modeleIGA._Jpqr,modeleIGA._dim, 1, 
                                           num_patch=0, offset=1)
icps1 = manip.get_boundCPindice_wEdges(modeleIGA._Nkv,modeleIGA._Jpqr,modeleIGA._dim, 1, 
                                           num_patch=1, offset=1)
for i in np.arange(0,3):
    manip.add_displacementBC(modeleIGA, modeleIGA._indCPbyPatch[0][icps0], i+1, 0.)
    manip.add_displacementBC(modeleIGA, modeleIGA._indCPbyPatch[1][icps1], i+1, 0.)


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
SOL,u = rsol.reconstruction(*modeleIGA.get_inputs4solution(x))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'coupling_mono',SOL.transpose(),nb_ref=np.array([3,3,1]),
    flag=np.array([True,True,False])))

print(('Total time for Mono analysis : %.2f s' % (time.time() - ti)))

print(np.shape(u))
print('shape(sol): ', np.shape(SOL))
print(modeleIGA._mcrd)
print(modeleIGA._nb_cp)


# compute angle 
from postprocessing.postproc import displayangledeformconfig
    
displayangledeformconfig(*modeleIGA.get_inputs4postprocCPLG('test',SOL.transpose(),nb_ref=10))
    
data = np.loadtxt('results/test.txt',delimiter=',')
n  = int(data.shape[0]/2)
print(n)
v1 = data[:n,1:]
v2 = data[n:,1:]
    
deg = np.rad2deg(np.arccos(np.sum(v1*v2,axis=1)))
err = np.abs(deg-90.)/90.
newdata = np.array([data[:n,0],deg,err]).transpose()
np.savetxt('results/rotTbeam%i.txt'%r,newdata,delimiter=',')



exit()

# TODO : check angle between patches at interface

pp.post_cplg(*modeleIGA.get_inputs4postprocCPLG('test',SOL.transpose(),nb_ref=0))

a1 = np.loadtxt('results/test3.res')
a2 = np.loadtxt('results/test4.res')


print(a2[:, 9:])

import matplotlib.pyplot as plt

plt.plot(a1[:,3], a1[:,9], label='a1 - 9')
plt.plot(a1[:,3], a1[:,10], label='a1 - 10')
plt.plot(a2[:,3], a2[:,9], label='a2 - 9')
plt.plot(a2[:,3], a2[:,10], label='a2 - 10')
plt.show()

plt.plot(a1[:,3], np.rad2deg(a1[:,9]), label='a1 - 9')
plt.plot(a1[:,3], np.rad2deg(a1[:,10]), label='a1 - 10')
plt.show()


exit()


if False:
    # energy
    energy = x.dot( K2solve.dot(x) )
    print('energy:',energy)







if False:
    # sensitivity analysis
    if modeleIGA._mcrd == 2:
        SOL3D = np.c_[SOL,np.zeros(SOL.shape[0])]
    else:
        SOL3D = SOL
        
    from optim.gradcompliance import gradcomp_an, gradcplg_an
    dComp_fine = gradcomp_an( 
        *modeleIGA.get_inputs4gradCompliance( SOL3D.T ) )
    dComp_fine+= gradcplg_an( 
        *modeleIGA.get_inputs4gradCoupling( SOL3D.T ) )
    
    modeleIGA.generate_vtk4controlMeshVisu(
        'sensitivityAnalysis',0,sol=dComp_fine )













if False:
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
                            M=None, tol=1.e-8, maxiter=300,savetxt=True)
    print(' Resolution with PCPG algo')
    print((' (duration : %.2f s, nb iter : %i).\n\n' % (time.time() - t1,nbiter)))
    
    
    # - postproc
    alpha = invGG.dot( G.T * (t - Sd(lmbda) ) )
    utot = np.zeros(modeleIGA._nb_dof_tot)
    utot[idof_internal] = LU.solve(f2solve - C2solve.T*lmbda) + R*alpha
    
    SOL,u = rsol.reconstruction(*modeleIGA.get_inputs4solution(utot[idof]))
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'coupling_pcpg',SOL.transpose(),nb_ref=np.array([4,4,1]),
        flag=np.array([True,True,False])))




## compute angle 
#from postprocessing.postproc_cplg import displayangledeformconfig
    
#displayangledeformconfig(*modeleIGA.get_inputs4postprocCPLG('test',SOL.transpose(),nb_ref=10))
    
#data = np.loadtxt('results/test.txt',delimiter=',')
#n  = data.shape[0]/2
#v1 = data[:n,1:]
#v2 = data[n:,1:]
    
#deg = np.rad2deg(np.arccos(np.sum(v1*v2,axis=1)))
#err = np.abs(deg-90.)/90.
#newdata = np.array([data[:n,0],deg,err]).transpose()
#np.savetxt('temp/rotTbeam%i.txt'%r,newdata,delimiter=',')
    
'''


# RESULTION Decomposition de Domaine
time2 = time.time()
ip1 = np.array([1])
ip2 = np.array([2])
ipl = np.array([7,8])

dof1 = np.array([],dtype=np.intp)
for p in ip1:
    dof1 = np.concatenate((dof1,np.repeat(np.unique(modeleIGA._IEN[p])-1,3)*3 
                           + np.tile(np.array([0,1,2]),modeleIGA._indCPbyPatch[p].size)))
dof1 = np.intersect1d(dof1,idof)

K1  = Ktot[dof1,:][:,dof1]
F1  = Fb[dof1]
LU1 = sp.linalg.splu(K1)

dof2 = np.array([],dtype=np.intp)
for p in ip2:
    dof2 = np.concatenate((dof2,np.repeat(np.unique(modeleIGA._IEN[p])-1,3)*3 
                           + np.tile(np.array([0,1,2]),modeleIGA._indCPbyPatch[p].size)))
#dof2 = np.intersect1d(dof2,idof)
K2   = Ktot[dof2,:][:,dof2]

iv2 = manip.get_vertexCPindice(modeleIGA._Nkv,modeleIGA._Jpqr,modeleIGA._dim,num_patch=ip2[0])
ddlr= np.array([(modeleIGA._indCPbyPatch[ip2[0]][iv2[0]]-1)*3+0,
                (modeleIGA._indCPbyPatch[ip2[0]][iv2[0]]-1)*3+1,
                (modeleIGA._indCPbyPatch[ip2[0]][iv2[0]]-1)*3+2,
                (modeleIGA._indCPbyPatch[ip2[0]][iv2[1]]-1)*3+0,
                (modeleIGA._indCPbyPatch[ip2[0]][iv2[1]]-1)*3+2,
                (modeleIGA._indCPbyPatch[ip2[0]][iv2[1]]  )*3+0])
dof2pp = np.setxor1d(dof2,ddlr)
K2pp = Ktot[dof2pp,:][:,dof2pp]
K2pr = Ktot[dof2pp,:][:,ddlr]

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

print '\n Time for DD solving : %.2f s\n\n' % (time.time() - time2)



# Postprocessing
SOL,u = rsol.reconstruction(*modeleIGA.get_inputs4solution(xDD))
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'coupling',SOL.transpose(),nb_ref=3*np.array([1,1,1]),
    flag=np.array([True,True,False])))

'''

