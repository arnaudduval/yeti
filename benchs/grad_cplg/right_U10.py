# Copyright 2021 Arnaud Duval

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

# Python module
import numpy as np
from numpy.lib.arraysetops import intersect1d
import scipy.sparse as sp
import sys
import time

#IGA module
from preprocessing.igaparametrization import IGAparametrization, OPTmodelling
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5
from utils import gausspts
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip
import preprocessing.igaparametrization.IGA_manipulation as igamanip

FILENAME='right_U10'

modeleIGA = IGAparametrization(filename=FILENAME)

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

# Embedded patch
nb_ref[:, 0] = np.array([0, 0, 0])
nb_deg[:, 0] = np.array([0, 0, 0])

# Hull
nb_ref[:, 1] = np.array([0, 0, 0])
nb_deg[:, 1] = np.array([0, 0, 0])


modeleIGA.refine(nb_ref, nb_deg, additional_knots)

if True:
    # STIFFNESS MATRIX
    ndof = modeleIGA._nb_dof_free
    idof = modeleIGA._ind_dof_free[:ndof]-1
    data, row, col, Fb = build_stiffmatrix( *modeleIGA.get_inputs4system_elemStorage() )
    Kside = sp.coo_matrix((data, (row,col)), shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                            dtype='float64').tocsc()
    Ktot = Kside+Kside.transpose()
    del Kside, data, row, col

    print("stiffness matrix assembly done")

    # Monolithic SOLVE
    t2 = time.time()
    K2solve = Ktot[idof,:][:,idof]
    LU = sp.linalg.splu(K2solve)
    x  = LU.solve(Fb[idof])
    t3 = time.time()

    print("Resolution done ", t3-t2, " seconds")

    # Postprocessing
    SOL,u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'couplingU5',SOL.transpose(),nb_ref=np.array([1,1,1]),
        Flag=np.array([True,True,False])))

initcoords = modeleIGA.coords

# control points controlling the hull

ind_cp_patch1_face2 = igamanip.get_boundCPindice_wEdges(
    modeleIGA._Nkv, modeleIGA._Jpqr, modeleIGA._dim, num_bound=2,
    num_patch=0, offset=0)
ind_cp_patch4_face1 = igamanip.get_boundCPindice_wEdges(
    modeleIGA._Nkv, modeleIGA._Jpqr, modeleIGA._dim, num_bound=1,
    num_patch=3, offset=0)

print(ind_cp_patch1_face2)
print(initcoords[:, modeleIGA._indCPbyPatch[0][ind_cp_patch1_face2]-1])
print('---')
print(ind_cp_patch4_face1)
print(initcoords[:, modeleIGA._indCPbyPatch[3][ind_cp_patch4_face1]-1])

print(modeleIGA._indCPbyPatch[0])

ind_cp_p4 = modeleIGA._indCPbyPatch[3][ind_cp_patch4_face1]-1
ind_cp_p1 = modeleIGA._indCPbyPatch[0][ind_cp_patch1_face2]-1

print("CPs patch 2")
print(ind_cp_p4 + 1)
print("CPs patch 1")
print(ind_cp_p1 + 1)

# exit()

NB_VAR = 3 * ind_cp_patch1_face2.size

def shapemodif(coords0, igapara, var):
    igapara.coords[:,:] = coords0[:,:]
    dim = 3

    i = 0
    # OTODO : utiliser un zip
    for i_cp, _ in enumerate(ind_cp_p1):
        for i_dim in range(dim):
            assert (igapara.coords[:, ind_cp_p1[i_cp]] == igapara.coords[:, ind_cp_p1[i_cp]]).all()
            igapara.coords[i_dim, ind_cp_p1[i_cp]] = initcoords[i_dim, ind_cp_p1[i_cp]] + var[i]
            igapara.coords[i_dim, ind_cp_p4[i_cp]] = initcoords[i_dim, ind_cp_p4[i_cp]] + var[i]
            i += 1


# Refinement from optim model to analysis model
nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)

optPB = OPTmodelling(modeleIGA, NB_VAR, shapemodif,
                     nb_degreeElevationByDirection=nb_deg,
                     nb_refinementByDirection=nb_ref)

x0 = np.zeros((NB_VAR))
c0 = optPB.compute_compliance_discrete(x0)
# v0 = optPB.compute_volume(x0, listpatch=[1, 1, 0])
gradC_DF = optPB.compute_gradCompliance_FD(x0)
gradC_AN = optPB.compute_gradCompliance_AN(x0)
# gradV_DF = optPB.compute_gradVolume_DF(x0, listpatch=[1, 1, 0])
# gradV_AN = optPB.compute_gradVolume_AN(x0, listpatch=[1, 1, 0])

for i in range(NB_VAR):
    print(i, gradC_DF[i], gradC_AN[i])

print(f"Compliance : {c0:.02E}")
# print(f"Volume : {v0:.02E}")
print(f"Erreur : {np.linalg.norm(gradC_DF-gradC_AN)/np.linalg.norm(gradC_DF)}")
