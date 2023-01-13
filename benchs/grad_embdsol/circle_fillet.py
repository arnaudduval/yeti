# Copyright 2021 Marie Guerder

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


import sys
import time

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization, OPTmodelling
from preprocessing.igaparametrization import IGAmanip as manip
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# =============================================================================
# Creation of the IGA object
# =============================================================================
# Read data and create IGAparametrization object
modeleIGA = IGAparametrization(filename='circle_fillet')

# Set arguments for model refinement
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
additional_knots = {"patches": np.array([]),
                    "1": np.array([]),
                    "2": np.array([]),
                    "3": np.array([])}

nb_deg[:, 1] = [0, 1, 0]
nb_ref[:, 1] = [2, 3, 2]    # Knot insertion, [xi, eta, zeta]

# Refine model
modeleIGA.refine(nb_ref, nb_deg, additional_knots)

if False:
    # =============================================================================
    # Static study
    # =============================================================================
    # Matrix assembly
    ndof = modeleIGA._nb_dof_free
    idof = modeleIGA._ind_dof_free[:ndof] - 1

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

    # RESOLUTION : direct solver
    t2 = time.time()
    LU = sp.linalg.splu(K2solve)
    x = LU.solve(Fb[idof])

    print(('\n Time for direct solving :    %.2f s\n\n' % (time.time() - t2)))

    # =============================================================================
    # Postprocessing
    # =============================================================================
    # Solution reconstruction
    SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

    # VTU file generation
    nb_ref_visu = np.array([1, 1, 1])       # Refinement for visu: [xi, eta, zeta]
    output = np.array([True, True, False])  # Output type: [disp, stress, VM]
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'circle-fillet',
        SOL.transpose(), nb_ref=3*nb_ref_visu, Flag=output))
    modeleIGA.generate_vtk4controlMeshVisu(
        'circle-fillet_map_cp', 0)
    modeleIGA.generate_vtk4controlMeshVisu(
        'circle-fillet_embded_cp', 1)

# CP of embedded patch not located at u=0 and w=0

indCPembd = modeleIGA._indCPbyPatch[1] - 1

not_u0_w0 = np.where((modeleIGA.coords[0, indCPembd] != 0.)
                     & (modeleIGA.coords[2, indCPembd] != 0.))[0]
ind_movable_cps = indCPembd[not_u0_w0]
for ipt in ind_movable_cps:
    print(modeleIGA.coords[:, ipt])

NB_VAR = ind_movable_cps.size * 2

print(ind_movable_cps)

initcoords = modeleIGA.coords

# Build optim model
def shapemodif(coords0, igapara, var):
    igapara.coords[:, :] = coords0[:, :]
    dim = coords0.shape[0]

    i = 0
    for i_cp in ind_movable_cps:
        for i_dim in [0,2]:
            igapara.coords[i_dim, i_cp] = initcoords[i_dim, i_cp] + var[i]
            i += 1

# Refinement from optim model to analysis model
nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)

optPB = OPTmodelling(modeleIGA, NB_VAR, shapemodif,
                     nb_degreeElevationByDirection=nb_deg,
                     nb_refinementByDirection=nb_ref)

x0 = np.zeros((NB_VAR))

c0 = optPB.compute_compliance_discrete(x0)

gradC_AN = optPB.compute_gradCompliance_AN(x0)
gradC_DF = optPB.compute_gradCompliance_FD(x0)

for i in range(NB_VAR):
    print(i, gradC_AN[i], gradC_DF[i])

error = np.linalg.norm(gradC_DF - gradC_AN) / c0

print(f"Compliance : {c0:.02E}")
print(f"Error : {error:.02E}")