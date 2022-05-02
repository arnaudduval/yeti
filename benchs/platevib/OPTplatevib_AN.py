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
Shape of a 3D solid plate is optimized versus 2 target eigenfrequencies
Full analytical gradients are used
validation is performed to verify is target frequencies are obtained
"""

# Python module
import numpy as np
import sys
import time

#IGA module
from preprocessing.igaparametrization import IGAparametrization, IGAmanip as manip
import postprocessing.postproc as pp
import reconstructionSOL as rsol

from preprocessing.igaparametrization import OPTmodelling


# Base path for .inp and .NB files
FILENAME = 'plateVolume'

# Creation of the IGA object
modeleIGA = IGAparametrization(filename=FILENAME)

# Refine IGA object along 1st and 2nd directions
nb_deg = np.array([1,1,0])
nb_ref = np.array([2,2,0])
modeleIGA.refine(nb_ref,nb_deg)


# Parametrization
# Get index of CP at top and bottom
botcps = manip.get_directionCP(modeleIGA, 5,0,0) - 1
topcps = manip.get_directionCP(modeleIGA, 6,0,0) - 1

nb_var = botcps.size

# Define thickness map depending on design variables
hM = 0.5
hm = 0.05
def thickness(coords0,igapara,var):
    
    var = (hM-hm)*var + hm
    
    igapara._COORDS[:,:]      = coords0[:,:]
    igapara._COORDS[2,botcps] =-0.5 * ( (hM-hm)*var + hm )
    igapara._COORDS[2,topcps] = 0.5 * ( (hM-hm)*var + hm )
    return None


# Refinement from optimization model to analysis model
nb_deg = np.array([0,0,1])
nb_ref = np.array([2,2,0])

# Define optim problem
# modeleIGA -> initial parametrization
# nb_var    -> number of design variables
# thickness -> function updating coordinates of the coarse (optim) model
# nb_degreeElevationByDirection -> degree elevation for each diretion
# nb_refinementByDirection -> refinement in each direction
optPB = OPTmodelling(modeleIGA, nb_var, thickness,
                     nb_degreeElevationByDirection = nb_deg, 
                     nb_refinementByDirection      = nb_ref)


# OPTIMIZATION

from scipy.optimize import minimize

# Initial values of design variables
x0    = 0.6*np.ones(nb_var)
# Compute initial 1st eigenvalue (m0) and eigenvector (v0)
m0,v0 = optPB.compute_vibrationMode(x0)

# Define target eigenvalue(s)
nb_frqCible = 2
xcible      = 0.2*np.ones(nb_var)
wcible,vcible = optPB.compute_vibrationMode(xcible,nb_frq=nb_frqCible)

# Define a function ans its gradient w.r.t design variables
# Gap with target eigenvalues
def frecIGA(xk):
    valsi,vecti = optPB.compute_vibrationMode(xk,nb_frq=nb_frqCible)
    mi = np.sum((valsi - wcible)**2.)
    return mi / np.sum(wcible**2.)

def gradFrecIGA(xk):
    gradF = np.zeros_like(xk)
    valsi,vecti = optPB.compute_vibrationMode(xk, nb_frq=nb_frqCible)
    gradFreq = optPB.compute_gradVibration_AN(xk, nb_frq=nb_frqCible)
    for i in np.arange(nb_frqCible):
        gradF[:] += 2.*gradFreq[:,i]*(valsi[i] - wcible[i])
    return gradF  / np.sum(wcible**2.)

# Define a function and its gradient w.r.t design variables
# Volume relative to initial volume
V0 = optPB.compute_volume(x0)
def volIGA(xk):
    global V0
    return optPB.compute_volume(xk)/V0
def gradVolIGA(xk):
    global V0
    return optPB.compute_gradVolume_AN(xk)/V0




tab = []
iopt= 0 

# Define callback function to be run at each optimization iteration (iopt)
def saveXk(xk):
    global iopt
    print(('\nIteration%i'%iopt))
    valsi,vecti = optPB.compute_vibrationMode(xk)
    tab.append([xk[0],valsi[0]])
    
    # - Plot Thickness
    ticknessfield = np.zeros_like(optPB._coarseParametrization._COORDS)
    ticknessfield[2,:] = np.tile(  optPB._coarseParametrization._COORDS[2,topcps]
                                   - optPB._coarseParametrization._COORDS[2,botcps], 
                                   (2,1)).flatten('C')
    pp.generatevtu(*optPB._coarseParametrization.get_inputs4postprocVTU(
        'vibOpt%0.2d'%iopt,ticknessfield,
        nb_ref=np.array([4,4,1]),Flag=np.array([True, False, False])))
    np.savetxt('results/cps%0.2d.txt'%iopt,optPB._coarseParametrization._COORDS.T,delimiter=',')

    # - Plot Analysis
    
    SOL,u = rsol.reconstruction(
        **optPB._fineParametrization.get_inputs4solution(vecti))
    pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
        'vibAn%0.2d'%iopt,  SOL.transpose(),
        nb_ref=np.array([2,2,2]), Flag=np.array([True, False, False])))
    
    iopt += 1
    
    return None

# Bounds for design variables : [0 ; 1]
bds = ((0.,1.),)*nb_var

# Run optimization

## 1 - get given 1st eigenvalue
res = minimize(frecIGA,x0,method='SLSQP',jac=gradFrecIGA,bounds=bds,callback=saveXk)

if not res['success'] or res['fun'] > 0.01:
    sys.exit(-1)
else:
    sys.exit(0)




## 2 - minimize volume with a a given 1st eigenvalue
#constraint_frq = ({'type': 'eq', 'fun':frecIGA, 'jac':gradFrecIGA})
#res = minimize(volIGA,x0,method='SLSQP',jac=gradVolIGA,bounds=bds,constraints=constraint_frq,
#               callback=saveXk,tol=1.e-4)


