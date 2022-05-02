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

The shape of a 2D solid tensile specimen is optimized versus its maximal Von Mises stress 
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
modeleIGA = IGAparametrization(filename='tensileSpecimen')

nb_degDV = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_refDV = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_degDV[0,1] = 2
nb_refDV[0,1] = 2

modeleIGA.refine(nb_refDV,nb_degDV)

# --
# Shape Parametrization

topcpsMid = manip.get_directionCP(modeleIGA, 4,1,0) - 1

nb_var = topcpsMid.size - 2
def heightMidPatch(coords0,igapara,var):
    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[1,topcpsMid[1:-1]] += var[:]
    return None


# --
# Build the optimization pb

from preprocessing.igaparametrization import OPTmodelling

nb_degAN = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_refAN = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

nb_degAN[:2,:] = 2
nb_refAN[:2,:] = 4
nb_degAN = np.maximum(nb_degAN-nb_degDV,0)
nb_refAN = np.maximum(nb_refAN-nb_refDV,0)

optPB = OPTmodelling(modeleIGA, nb_var, heightMidPatch,
                     nb_degreeElevationByDirection = nb_degAN, 
                     nb_refinementByDirection      = nb_refAN)


# --
# Initialization and Definition of the objective and constraints (using nlopt) 

Pnorm = 40
x0  = np.zeros(nb_var)
phi0= optPB.compute_vonmisesAggreg(x0,pnorm=Pnorm)

i = 0; iplt = 0
currentmin = phi0+0.1
ref_plot = np.array([3,3,1])
Output = np.array([True, True, False])


import nlopt

def vonmisesAggregate(xVM,gradVM):
    phi = optPB.compute_vonmisesAggreg(xVM,pnorm=Pnorm)
    if gradVM.size>0:
        global i,currentmin,iplt
        i += 1
        print('\n--')
        print('Iter %3i' % i)
        gradVM[:] = optPB.compute_gradVonMisesAggreg_AN(xVM,pnorm=Pnorm)/phi0
        
        # postprocessing
        if phi<currentmin:
            iplt += 1
            currentmin = phi
            pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
                'OPT6-coarse%0.2d'%iplt,np.zeros_like(optPB._coarseParametrization._COORDS)[:2],
                nb_ref=2*ref_plot,Flag=Output))
            optPB._coarseParametrization.generate_vtk4controlMeshVisu('OPT6-coarse%0.2d'%iplt,0)
            
            SOL,u = rsol.reconstruction(
                **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
            pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
                'OPT6-fine%0.2d'%iplt,  SOL.transpose(), nb_ref=1*ref_plot, Flag=Output))
    return phi/phi0



minimize = nlopt.opt( nlopt.LD_SLSQP, nb_var )

minimize.set_min_objective( vonmisesAggregate )

minimize.set_ftol_rel(1.0e-06)
minimize.set_xtol_rel(1.0e-06)
minimize.set_maxeval(200)

minimize.set_lower_bounds(-5*np.ones(nb_var) )
minimize.set_upper_bounds( 5*np.ones(nb_var) )


# --
# Run optimization
x = minimize.optimize( x0 )

# Verify results
# Numerical reference result
x_ref = np.array([-2.36117255, -2.49828454, -1.97837879, -1.11597613, -0.40425204])

error = np.sqrt(sum((x-x_ref)**2.))

print(error)

if error > 1.e-5:
    sys.exit(-1)
else:
    sys.exit(0)
