# Copyright 2022 Arnaud Duval

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

# !/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import time

# yeti modules
from preprocessing.igaparametrization import IGAparametrization, IGAmanip as manip
import postprocessing.postproc as pp
import reconstructionSOL as rsol

from preprocessing.igaparametrization import OPTmodelling

modeleIGA = IGAparametrization(filename='centrif_U1_C0')

# Refinement to create optimisation model
nb_deg = np.zeros((3, modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch),dtype=np.intp)

nb_ref[:, 0] = np.array([0, 2, 0])
nb_deg[:, 0] = np.array([0, 1, 0])
additional_knots = {"patches": np.array([0]),
                    "1": np.array([]), "2": np.array([0.2,0.2]), "3": np.array([])}

#additional_knots = {"patches": np.array([0]),
#                    "1": np.array([]), "2": np.array([]), "3": np.array([])}


# Initial refinement (none)
modeleIGA.refine(nb_ref, nb_deg, additional_knots)

icp = np.where(modeleIGA._COORDS[1, :] > 1.2)[0]
nb_var = np.unique(modeleIGA._COORDS[1, icp]).size

# Define shape change from design variables
# min and max dimensions, assuming that design variables are in [0,1]
mindim = 0.1
maxdim = 3.5

def shapemodif(coords0, igapara, var):
    igapara._COORDS[:,:] = coords0[:,:]
    # shape change is made on points with y coord higher than 3
    icp = np.where(modeleIGA._COORDS[1, :] > 1.2)[0]
    i = 0
    for y in np.unique(modeleIGA._COORDS[1, icp]):
        # WARNING exact real value comparison is unsafe
        jcp = np.where(modeleIGA._COORDS[1, :] == y)[0]

        igapara._COORDS[2, jcp[0]] = - (mindim + var[i]*(maxdim-mindim))/2.
        igapara._COORDS[2, jcp[1]] = - (mindim + var[i]*(maxdim-mindim))/2.
        igapara._COORDS[2, jcp[2]] =   (mindim + var[i]*(maxdim-mindim))/2.
        igapara._COORDS[2, jcp[3]] =   (mindim + var[i]*(maxdim-mindim))/2.

        i += 1
    return None

# Refinement from optim model to analysis model
nb_deg = np.array([1, 0, 1])
nb_ref = np.array([2, 2, 2])# Define optim problem 

optPB = OPTmodelling(modeleIGA, nb_var, shapemodif,
                     nb_degreeElevationByDirection = nb_deg, 
                     nb_refinementByDirection      = nb_ref)

# OPTIMISATION
from scipy.optimize import minimize

# Initial values of design variables
x0 = ((1.5-mindim)/(maxdim-mindim))*np.ones(nb_var)
# Compute initial values
v0 = optPB.compute_volume(x0)
c0 = optPB.compute_compliance_discrete(x0)

# Define functions and gradients
def volIGA(xk):
    global v0
    return (optPB.compute_volume(xk)-v0)/v0
def gradVolIGA(xk):
    global v0
    return optPB.compute_gradVolume_AN(xk)/v0

def compIGA(xk):
    global c0
    return optPB.compute_compliance_discrete(xk)/c0

def gradCompIGA(xk):
    global c0
    return optPB.compute_gradCompliance_AN(xk)/c0
    
iopt = 0

# Define callback function
def saveXk(xk):
    global iopt
    print(('\nIteration%i'%iopt))
    SOL,u = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
    pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            'OPT3-fine%0.2d'%iopt,  SOL.transpose(), nb_ref=3*np.array([1,1,1]),
            Flag=np.array([True, False, False])))
    
    iopt += 1
    
    return None

# Bounds for design variables
bds = ((0., 1.),)*nb_var

constraint = ({'type': 'eq', 'fun':volIGA, 'jac':gradVolIGA})
x0 = ((1.5-mindim)/(maxdim-mindim))*np.ones(nb_var)

res = minimize(compIGA, x0, method='SLSQP',
               jac=gradCompIGA, bounds=bds,
               constraints=constraint, callback=saveXk)

# Verify results
# Numerical reference result
x_ref = np.array([7.16402311e-01, 6.94203856e-01, 3.37998525e-01, 0.0, 0.0])

error = np.sqrt(sum((res['x']-x_ref)**2.))

print(error)

if error > 1.e-5:
    sys.exit(-1)
else:
    sys.exit(0)