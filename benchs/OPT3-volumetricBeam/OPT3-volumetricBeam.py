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

The shape of a solid 3D beam is optimized versus its compliance
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
modeleIGA = IGAparametrization(filename='volumetricBeam')

nb_degDV = np.array([1,0,0])
nb_refDV = np.array([4,0,0])

modeleIGA.refine(nb_refDV,nb_degDV)


# --
# Shape Parametrization

botcps  = manip.get_directionCP(modeleIGA, 3,0,0)-1
topcps  = manip.get_directionCP(modeleIGA, 4,0,0)-1
rightcps= manip.get_directionCP(modeleIGA, 5,0,0)-1
leftcps = manip.get_directionCP(modeleIGA, 6,0,0)-1
edgeR   = np.intersect1d(topcps,rightcps)
edgeL   = np.intersect1d(topcps, leftcps)

hM = 10.0
hm =  1.5
def beamHeight(coords0,igapara,var):
    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[1,topcps] = (hM-hm)*var[:] + hm
    return None

wM = 15./2
wm = 1.
def beamWidth(coords0,igapara,var):
    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[2,edgeR] = wM - ( (wM-wm)*var[:edgeR.size] + wm )
    igapara._COORDS[2,edgeL] = wM + ( (wM-wm)*var[edgeR.size:] + wm )
    return None

nb_var = topcps.size * 2
def fullShapeParametrization(coords0,igapara,var):
    beamHeight(coords0,igapara,var[:topcps.size])
    coords1 = igapara._COORDS.copy()
    beamWidth( coords1,igapara,var[topcps.size:])
    return None


# --
# Build the optimization pb

from preprocessing.igaparametrization import OPTmodelling

nb_degAN = np.maximum(np.array([1,1,1])-nb_degDV,0)
nb_refAN = np.maximum(np.array([5,3,3])-nb_refDV,0)

optPB = OPTmodelling(modeleIGA, nb_var, fullShapeParametrization,
                     nb_degreeElevationByDirection = nb_degAN, 
                     nb_refinementByDirection      = nb_refAN)


# --
# Initialization and Definition of the objective and constraints (using nlopt) 

x0 = np.block([4.5/8.5*np.ones(topcps.size), 1.5/6.5*np.ones(topcps.size)])
V0 = optPB.compute_volume(x0) 
c0 = optPB.compute_compliance_discrete(x0)

i  = 0
ref_plot = np.array([2,2,2])
Output = np.array([True, True, False])


import nlopt

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
            'OPT3-coarse%0.2d'%i,np.zeros_like(optPB._coarseParametrization._COORDS),
            nb_ref=2*ref_plot,flag=np.array([False]*3)))
        optPB._coarseParametrization.generate_vtk4controlMeshVisu('OPT3-coarse%0.2d'%i,0)
        
        SOL,u = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            'OPT3-fine%0.2d'%i,  SOL.transpose(), nb_ref=1*ref_plot, flag=Output))
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


# --
# Run optimization
x = minimize.optimize( x0 )


# Verify result
# Reference numerical results
x_ref = np.array([1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                  1.00000000e+00, 1.00000000e+00, 9.51688643e-01, 8.57319065e-01,
                  7.56652137e-01, 6.47480521e-01, 5.37205520e-01, 4.27196601e-01,
                  3.16146907e-01, 2.06193610e-01, 9.83090028e-02, 0.00000000e+00,
                  0.00000000e+00, 0.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                  1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00,
                  9.51688643e-01, 8.57319065e-01, 7.56652137e-01, 6.47480521e-01,
                  5.37205520e-01, 4.27196602e-01, 3.16146907e-01, 2.06193610e-01,
                  9.83090030e-02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                  3.68610530e-01, 2.77434331e-01, 2.23872781e-01, 1.47350801e-01,
                  6.06727437e-02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                  7.93128500e-18, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                  2.67617083e-18, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                  6.15195708e-18, 0.00000000e+00, 3.68610530e-01, 2.77434331e-01,
                  2.23872781e-01, 1.47350801e-01, 6.06727437e-02, 0.00000000e+00,
                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 4.87211811e-17,
                  0.00000000e+00, 0.00000000e+00, 4.60362610e-18, 0.00000000e+00])

error = np.sqrt(sum((x-x_ref)**2.))

print(error)

if error > 1.e-5:
    sys.exit(-1)
else:
    sys.exit(0)



