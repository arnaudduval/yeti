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

The shape of a solid 2D beam is optimized versus its maximal displacement
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
modeleIGA = IGAparametrization(filename='cantileverBeam')

nb_degDV = np.array([1,0,0])
nb_refDV = np.array([2,0,0])

modeleIGA.refine(nb_refDV,nb_degDV)


# --
# Shape Parametrization

topcps  = manip.get_directionCP(modeleIGA, 4,0,0)-1

nb_var = topcps.size
hM = 10.0
hm =  1.5
def beamHeight(coords0,igapara,var):
    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[1,topcps] = (hM-hm)*var[:] + hm
    return None




# --
# Build the optimization pb

from preprocessing.igaparametrization import OPTmodelling

nb_degAN = np.maximum(np.array([1,1,0])-nb_degDV,0)
nb_refAN = np.maximum(np.array([5,3,0])-nb_refDV,0)

optPB = OPTmodelling(modeleIGA, nb_var, beamHeight,
                     nb_degreeElevationByDirection = nb_degAN, 
                     nb_refinementByDirection      = nb_refAN)


# --
# Initialization and Definition of the objective and constraints (using nlopt) 

x0 = 4.5/8.5*np.ones(topcps.size)
V0 = 210.
xi4disp = np.array([0.75,0.,0.])
n0 = np.linalg.norm( optPB.compute_displacement(x0, xi4disp) )

i  = 0
ref_plot = np.array([2,2,2])
Output = np.array([True, True, False])


import nlopt

def dispNorm(xD,gradD):
    di = optPB.compute_displacement(xD, xi4disp)
    ni = np.linalg.norm( di )
    if gradD.size>0:
        global i
        i += 1
        print('\n--')
        print('Iter %3i' % i)
        gradD[:] = optPB.compute_gradDisplacement_AN(xD, xi4disp).dot(di)/ni/n0
        
        # postprocessing
        pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
            'OPT4-coarse%0.2d'%i,np.zeros((2,optPB._coarseParametrization._nb_cp)),
            nb_ref=2*ref_plot,Flag=np.array([False]*3)))
        optPB._coarseParametrization.generate_vtk4controlMeshVisu('OPT4-coarse%0.2d'%i,0)
        
        SOL,u = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            'OPT4-fine%0.2d'%i,  SOL.transpose(), nb_ref=1*ref_plot, Flag=Output))
    return ni/n0


def vol(xV,gradV):
    if gradV.size>0:
        gradV[:] = optPB.compute_gradVolume_AN(xV)/V0
    return optPB.compute_volume(xV)/V0 - 1.




minimize = nlopt.opt( nlopt.LD_SLSQP, nb_var )

minimize.set_min_objective( dispNorm )
minimize.add_inequality_constraint( vol, 1e-5 )

minimize.set_ftol_rel(1.0e-08)
minimize.set_xtol_rel(1.0e-08)
minimize.set_maxeval(50)

minimize.set_lower_bounds( np.zeros(nb_var) )


# --
# Run optimization
x = minimize.optimize( x0 )

# Verify results
# Numerical reference result
x_ref = np.array([1.16393910e+00, 1.08633243e+00, 9.03526441e-01,
                  5.72507534e-01, 0.00000000e+00, 3.87531924e-17])

error = np.sqrt(sum((x-x_ref)**2.))

print(error)

if error > 1.e-5:
    sys.exit(-1)
else:
    sys.exit(0)

