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

The shape of a 2D solid plate with a hole is optimized versus its compliance
Volume is kept constant
Resutling hole shape is compared to a circle

"""

# Python module
import numpy as np
import sys
import time
import math

#IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, IGAmanip as manip
import yeti_iga.postprocessing.postproc as pp
import yeti_iga.reconstructionSOL as rsol


# Creation of the IGA object
# --------------------------
modeleIGA = IGAparametrization(filename='plateWithHole')


# --
# Shape Parametrization

nb_var = 6
def holeshape(coords0,igapara,var):

    moveX = np.arange(0,3, dtype=np.intp)
    moveY = np.arange(1,4, dtype=np.intp)

    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[0,moveX] += var[:moveX.size]
    igapara._COORDS[1,moveY] += var[moveX.size:]
    return None


# --
# Build the optimization pb

from yeti_iga.preprocessing.igaparametrization import OPTmodelling

nb_deg = np.array([1,1,0])
nb_ref = np.array([3,4,0])

optPB = OPTmodelling(modeleIGA, nb_var, holeshape,
                     nb_degreeElevationByDirection = nb_deg,
                     nb_refinementByDirection      = nb_ref)


# --
# Initialization and Definition of the objective and constraints (using nlopt)

x0 = np.zeros(nb_var)
V0 = optPB.compute_volume(x0)
c0 = optPB.compute_compliance_discrete(x0)

i  = 0 # count grad evaluation
ref_plot = np.array([3,3,1])
Output = np.array([True, True, False])


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
            'OPT1-coarse%0.2d'%i,np.zeros_like(optPB._coarseParametrization._COORDS)[:2],
            nb_ref=2*ref_plot,Flag=np.array([False]*3) ))
        optPB._coarseParametrization.generate_vtk4controlMeshVisu('OPT1-coarse%0.2d'%i,0)

        SOL,u = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            'OPT1-fine%0.2d'%i,  SOL.transpose(), nb_ref=1*ref_plot, Flag=Output))
    return ci

def vol(xV,gradV):
    if gradV.size>0:
        gradV[:] = optPB.compute_gradVolume_AN(xV)/V0
    return optPB.compute_volume(xV)/V0 - 1.




minimize = nlopt.opt( nlopt.LD_SLSQP, nb_var )

minimize.set_min_objective( comp )
minimize.add_inequality_constraint( vol, 1e-5 )

minimize.set_ftol_rel(1.0e-08)
minimize.set_xtol_rel(1.0e-08)
minimize.set_maxeval(50)

minimize.set_lower_bounds(-0.50*np.ones(nb_var))
minimize.set_upper_bounds( 0.50*np.ones(nb_var))


# --
# Run optimization
x = minimize.optimize( x0 )


# Verify result
# Since volume should be kept constant, optimum is a hole with radius sqrt(2./pi)
# reference value for design variables
radius = math.sqrt(2.0/math.pi)

x_ref = np.array([1. - radius,
                  0.7 - radius,
                  0.3 - radius*(math.sqrt(2.0) - 1.0),
                  radius*(math.sqrt(2.0) - 1.0) - 0.3,
                  radius - 0.7,
                  radius - 1.])


error =sum((x-x_ref)**2.)

print(error)

if error > 1.e-5:
    sys.exit(-1)
else:
    sys.exit(0)


