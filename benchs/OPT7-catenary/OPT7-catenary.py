# Copyright 2020 Thibaut Hirschler
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
This cas is described in the following publication : 
Hirschler, T., Bouclier, R., Duval, A. et al. 
A New Lighting on Analytical Discrete Sensitivities in the Context of IsoGeometric Shape Optimization. 
Arch Computat Methods Eng (2020). https://doi.org/10.1007/s11831-020-09458-6

The shape of a Kirchhoff-Love shell catenary is optimized versus its maximal bending moment 
Volume is kept constant
Resuting shape is compared to reference numerical results

"""


# Python module
import numpy as np
import sys
import time

#IGA module
from preprocessing.igaparametrization import IGAparametrization, IGAmanip as manip
import postprocessing.postproc as pp
import reconstructionSOL as rsol


# Selection of .INP and .NB file
# ------------------------------
modeleIGA = IGAparametrization(filename='shellArch')

nb_degDV = np.array([1,0,0])
nb_refDV = np.array([2,0,0])

modeleIGA.refine(nb_refDV,nb_degDV)


# --
# Parametrization

nb_var = int(modeleIGA._Nkv[0,0] - modeleIGA._Jpqr[0,0] - 3)
def altitude(coords0,igapara,var):
    igapara._COORDS[:,:] = coords0[:,:]
    count = 1
    for v in var:
        icps = manip.get_boundCPindice_wEdges(igapara._Nkv,igapara._Jpqr,igapara._dim, 1, 
                                              num_patch=0, offset=count)
        igapara._COORDS[2,igapara._indCPbyPatch[0][icps]-1] += v*10.
        count += 1
    return None


# --
# Build the optimization pb

from preprocessing.igaparametrization import OPTmodelling

nb_degAN = np.maximum(np.array([2,2,0])-nb_degDV,0)
nb_refAN = np.maximum(np.array([5,1,0])-nb_refDV,0)

optPB = OPTmodelling(modeleIGA, nb_var, altitude,
                     nb_degreeElevationByDirection = nb_degAN, 
                     nb_refinementByDirection      = nb_refAN)


# --
# Initialization and Definition of the objective and constraints (using nlopt)

Pnorm = 40
x0 = 0.1*np.ones(nb_var)
phi0 = optPB.compute_stressAggreg(x0,pnorm=Pnorm)[3]
V0 =  7.5

i = 0
ref_plot = np.array([2,2,1])
Output   = np.array([True, True, False])

import nlopt

def momentAggregate(xM,gradM):
    phi = optPB.compute_stressAggreg(xM,pnorm=Pnorm)[3]
    if gradM.size>0:
        global i
        i += 1
        print('\n--')
        print('Iter %3i' % i)
        gradM[:] = optPB.compute_gradStressAggreg_AN(xM,pnorm=Pnorm)[:,3]/phi0
        
        # postprocessing
        pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
            'OPT7-coarse%0.2d'%i,np.zeros_like(optPB._coarseParametrization._COORDS),
            nb_ref=2*ref_plot,Flag=Output))
        optPB._coarseParametrization.generate_vtk4controlMeshVisu('OPT7-coarse%0.2d'%i,0)
        
        SOL,u = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            'OPT7-fine%0.2d'%i,  SOL.transpose(), nb_ref=1*ref_plot, Flag=Output))
    return phi/phi0


def vol(xV,gradV):
    if gradV.size>0:
        gradV[:] = optPB.compute_gradVolume_AN(xV)/V0
    return optPB.compute_volume(xV)/V0 - 1.


minimize = nlopt.opt( nlopt.LD_SLSQP, nb_var )


minimize.set_min_objective( momentAggregate )
minimize.add_equality_constraint( vol, 1e-6 )

minimize.set_ftol_rel(1.0e-06)
minimize.set_xtol_rel(1.0e-06)
minimize.set_maxeval(100)

minimize.set_lower_bounds( 0.*np.ones(nb_var) )
minimize.set_upper_bounds( 1.*np.ones(nb_var) )


# --
# Run optimization
x = minimize.optimize( x0 )

# Verify results
# Numerical reference result
x_ref = np.array([0.29278922, 0.502033, 0.50203325, 0.29278943])

error = np.sqrt(sum((x-x_ref)**2.))

print(error)

if error > 1.e-5:
    sys.exit(-1)
else:
    sys.exit(0)

