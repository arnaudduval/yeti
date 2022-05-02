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

The shape of Kirchhoff-Love shell roof is optimized versus its compliance
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
modeleIGA = IGAparametrization(filename='squareShellPlate')

nb_degDV = np.array([1,1,0])
nb_refDV = np.array([2,2,0])

modeleIGA.refine(nb_refDV,nb_degDV)


# --
# Shape Parametrization

vertex = manip.get_vertexCPindice(modeleIGA._Nkv, modeleIGA._Jpqr, modeleIGA._dim)[:4]
freecp = np.setxor1d(np.arange(0,modeleIGA._nb_cp),vertex)
nb_var = freecp.size
def altitude(coords0,igapara,var):
    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[2,freecp] = coords0[2,freecp] + 5.*var[:]
    return None


# --
# Build the optimization pb

from preprocessing.igaparametrization import OPTmodelling

nb_degAN = np.maximum(np.array([1,1,0])-nb_degDV,0)
nb_refAN = np.maximum(np.array([5,5,0])-nb_refDV,0)

optPB = OPTmodelling(modeleIGA, nb_var, altitude,
                     nb_degreeElevationByDirection = nb_degAN, 
                     nb_refinementByDirection      = nb_refAN)


# --
# Initialization and Definition of the objective and constraints (using nlopt) 

V0 = 1.10*optPB.compute_volume(np.zeros(nb_var))
c0 = optPB.compute_compliance_discrete(np.zeros(nb_var))
x0 = np.ones(nb_var)*0.1

i  = 0
ref_plot = np.array([2,2,1])
Output = np.array([True, False, False])


import nlopt

def comp(xC,gradC):
    ci = optPB.compute_compliance_discrete(xC)/c0
    if gradC.size>0:
        global i
        i += 1
        print('\n--')
        print('Iter %3i' % i)
        gradC[:] = optPB.compute_gradCompliance_AN(xC)/c0

        # postprocessing
        pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
            'OPT2-coarse%0.2d'%i,optPB._coarseParametrization._COORDS-optPB._initialCOORDS,
            nb_ref=2*ref_plot,Flag=Output))
        optPB._coarseParametrization.generate_vtk4controlMeshVisu('OPT2-coarse%0.2d'%i,0)
        
        SOL,u = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            'OPT2-fine%0.2d'%i,  SOL.transpose(), nb_ref=1*ref_plot, Flag=Output))
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

minimize.set_lower_bounds( 0.*np.ones(nb_var))
minimize.set_upper_bounds( 1.*np.ones(nb_var) )


# --
# Run optimization
x = minimize.optimize( x0 )

# Verify results
# Numerical reference result
x_ref = np.array([0.15241722, 0.54611007, 0.54611007, 0.15241722, 0.15241722, 0.25461606,
                  0.52285965, 0.52285965, 0.25461606, 0.15241722, 0.54611007, 0.52285965,
                  0.65837537, 0.65837537, 0.52285965, 0.54611007, 0.54611007, 0.52285965,
                  0.65837537, 0.65837537, 0.52285965, 0.54611007, 0.15241722, 0.25461606,
                  0.52285965, 0.52285965, 0.25461606, 0.15241722, 0.15241722, 0.54611007,
                  0.54611007, 0.15241722])

error = np.sqrt(sum((x-x_ref)**2.))

print(error)

if error > 1.e-5:
    sys.exit(-1)
else:
    sys.exit(0)


