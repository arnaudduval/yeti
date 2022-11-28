# Copyright 2022 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>

# -*- coding: utf-8 -*-

import sys
import numpy as np
import scipy.sparse as sp
import nlopt

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import OPTmodelling
from preprocessing.igaparametrization import IGAmanip as manip
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp


modeleIGA = IGAparametrization(filename='embd_volumetric_beam')

# Set arguments for model refinement
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
additional_knots = {"patches": np.array([]),
                    "1": np.array([]),
                    "2": np.array([]),
                    "3": np.array([])}

nb_ref[:, 1] = [4, 0, 0]
nb_deg[:, 1] = [1, 0, 0]

modeleIGA.refine(nb_ref, nb_deg, additional_knots)

# Shape parametrization
botcps = manip.get_directionCP(modeleIGA, 3, 1, 0)-1
topcps = manip.get_directionCP(modeleIGA, 4, 1, 0)-1
rightcps = manip.get_directionCP(modeleIGA, 5, 1, 0)-1
leftcps = manip.get_directionCP(modeleIGA, 6, 1, 0)-1
edgeR = np.intersect1d(topcps, rightcps)
edgeL = np.intersect1d(topcps, leftcps)

# Defined on embedded entity
hM = 0.1
hm = 0.015


def beamHeight(coords0, igapara, var):
    igapara._COORDS[:, :] = coords0[:, :]
    igapara._COORDS[1, topcps] = (hM - hm)*var[:] + hm
    return None

wM = 0.15/2.
wm = 0.01

def beamWidth(coords0, igapara, var):
    igapara._COORDS[:, :] = coords0[:, :]
    igapara._COORDS[2, edgeR] = wM - ((wM-wm)*var[:edgeR.size] + wm)
    igapara._COORDS[2, edgeL] = wM + ((wM-wm)*var[edgeR.size:] + wm)
    return None

nb_var = topcps.size * 2


def fullShapeParametrization(coords0, igapara, var):
    beamHeight(coords0, igapara, var[:topcps.size])
    coords1 = igapara._COORDS.copy()
    beamWidth(coords1, igapara, var[topcps.size:])
    return None


# Build optim problem
nb_degAN = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_refAN = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_degAN[:, 1] = np.maximum(np.array([1, 1, 1])-nb_deg[:, 1], 0)
nb_refAN[:, 1] = np.maximum(np.array([5, 3, 3])-nb_ref[:, 1], 0)

optPB = OPTmodelling(modeleIGA, nb_var, fullShapeParametrization,
                     nb_degreeElevationByDirection=nb_degAN,
                     nb_refinementByDirection=nb_refAN)

x0 = np.block([4.5/8.5*np.ones(topcps.size), 1.5/6.5*np.ones(topcps.size)])
V0 = optPB.compute_volume(x0, listpatch=[0, 1])
c0 = optPB.compute_compliance_discrete(x0)


def comp(xC,gradC):
    ci = optPB.compute_compliance_discrete(xC)/c0
    if gradC.size>0:
        global i,saveX
        i += 1
        print('\n--')
        print('Iter %3i' % i)
        gradC[:] = optPB.compute_gradCompliance_AN(xC)/c0
        
        # postprocessing
        pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
            'OPT-coarse%0.2d'%i,np.zeros_like(optPB._coarseParametrization._COORDS),
            nb_ref=2*ref_plot,Flag=np.array([False]*3)))
        optPB._coarseParametrization.generate_vtk4controlMeshVisu('OPT-coarse%0.2d'%i,0)
        
        SOL,u = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            'OPT-fine%0.2d'%i,  SOL.transpose(), nb_ref=1*ref_plot, Flag=Output))
    return ci

def vol(xV,gradV):
    if gradV.size>0:
        gradV[:] = optPB.compute_gradVolume_AN(xV)/V0
    return optPB.compute_volume(xV)/V0 - 1.

minimize = nlopt.opt( nlopt.LD_SLSQP, nb_var )

minimize.set_min_objective( comp )
minimize.add_inequality_constraint( vol, 1e-5 )

minimize.set_ftol_rel(1.0e-06)
minimize.set_xtol_rel(1.0e-06)
minimize.set_maxeval(200)

minimize.set_lower_bounds( np.zeros(nb_var) )
minimize.set_upper_bounds( np.ones( nb_var) )

# Run optimization
x = minimize.optimize( x0 )

exit()



initcoords = modeleIGA._COORDS


# Build optim model
def shapemodif(coords0, igapara, var):
    igapara._COORDS[:, :] = coords0[:, :]
    # shape change is made on nodes 10 and 12
    igapara._COORDS[0, 9] = initcoords[0, 9] + 0.5*var[0]
    igapara._COORDS[1, 9] = initcoords[1, 9] + 0.5*var[1]
    igapara._COORDS[2, 9] = initcoords[2, 9] + 0.5*var[2]
    igapara._COORDS[0, 11] = initcoords[0, 11] + 0.5*var[3]
    igapara._COORDS[1, 11] = initcoords[1, 11] + 0.5*var[4]
    igapara._COORDS[2, 11] = initcoords[2, 11] + 0.5*var[5]

    return None


# Refinement from optim model to analysis model
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

nb_var = 6

optPB = OPTmodelling(modeleIGA, nb_var, shapemodif,
                     nb_degreeElevationByDirection=nb_deg,
                     nb_refinementByDirection=nb_ref)

x0 = np.zeros((6))
print(optPB.compute_volume(x0, listpatch=[0, 1]))

print(optPB.compute_compliance_discrete(x0))

print(optPB.compute_gradCompliance_AN(x0))


# Static study
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof] - 1
print(idof)
print(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_free)
data, row, col, Fb = build_stiffmatrix(
                        *modeleIGA.get_inputs4system_elemStorage())
Kside = sp.coo_matrix((data, (row, col)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside + Kside.transpose()
K2solve = Ktot[idof, :][:, idof]
del Kside, data, row, col

x = sp.linalg.spsolve(K2solve, Fb[idof])

# Solution reconstruction
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    'test_vol_beam',
    SOL.transpose(), nb_ref=[3, 3, 3], Flag=[True, False, False]))
modeleIGA.generate_vtk4controlMeshVisu(
    'test_vol_beam_map_cp', 0)
modeleIGA.generate_vtk4controlMeshVisu(
    'test_vol_beam_embded_cp', 1)
