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
Test analytical gradients computation in the case of shape optimization of a 3D solid plate
versus 2 target eigenfrequencies
"""

import os
import numpy as np

#IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, IGAmanip as manip
from yeti_iga.preprocessing.igaparametrization import OPTmodelling

N_TARGET_FREQUENCIES = 2

# Creation of the IGA object
script_dir = os.path.dirname(os.path.realpath(__file__))
modeleIGA = IGAparametrization(filename=f'{script_dir}/plateVolume')

# Refine IGA object along 1st and 2nd directions
nb_deg = np.array([1,1,0])
nb_ref = np.array([1,1,0])
modeleIGA.refine(nb_ref,nb_deg)


# Parametrization
# Get index of CP at top and bottom
botcps = manip.get_directionCP(modeleIGA, 5,0,0) - 1
topcps = manip.get_directionCP(modeleIGA, 6,0,0) - 1

nb_var = botcps.size

# Define thickness map depending on design variables
H_MAX = 0.5
H_MIN = 0.05
def thickness(coords0,igapara,var):
    """
    Change thickness depending on design variables values
    """

    var = (H_MAX - H_MIN)*var + H_MIN

    igapara.coords[:,:]      = coords0[:,:]
    igapara.coords[2,botcps] = -0.5 * ( (H_MAX - H_MIN)*var + H_MIN )
    igapara.coords[2,topcps] = 0.5 * ( (H_MAX - H_MIN)*var + H_MIN )


# Refinement from optimization model to analysis model
nb_deg = np.array([0,0,1])
nb_ref = np.array([2,2,0])

# Define optim problem
optPB = OPTmodelling(modeleIGA, nb_var, thickness,
                     nb_degreeElevationByDirection = nb_deg,
                     nb_refinementByDirection      = nb_ref)


# Test values for design variables
xk_test = np.array([0.12009086, 0.31883724, 0.54778358, 0.96904667, 0.10729277,
       0.0791848 , 0.16109069, 0.43826168, 0.10729277, 0.07918475,
       0.16109061, 0.43826165, 0.12009085, 0.31883717, 0.54778359,
       0.96904703])

print(f'{optPB.compute_gradVolume_AN(xk_test) = }')
print(f'{optPB.compute_gradVibration_AN(xk_test, nb_frq=N_TARGET_FREQUENCIES) = }')

# Reference value of volume gradient
ref_gradvol = np.array([14.0625, 28.125 , 28.125 , 14.0625, 28.125 , 56.25  , 56.25  ,
       28.125 , 28.125 , 56.25  , 56.25  , 28.125 , 14.0625, 28.125 ,
       28.125 , 14.0625])

# Reference value of eigenvalues gradient
ref_gradvib = np.array([[ 0.00525155,  0.02234445],
       [ 0.01004653,  0.06093687],
       [ 0.00397719,  0.02923913],
       [-0.00086381, -0.00592409],
       [ 0.00950153,  0.01537288],
       [ 0.01565865,  0.0647121 ],
       [ 0.0036758 ,  0.07163502],
       [-0.00317898,  0.04879141],
       [ 0.00950153,  0.01537288],
       [ 0.01565865,  0.06471209],
       [ 0.0036758 ,  0.07163501],
       [-0.00317898,  0.04879141],
       [ 0.00525155,  0.02234444],
       [ 0.01004652,  0.06093686],
       [ 0.00397719,  0.02923912],
       [-0.00086381, -0.00592409]])

assert np.allclose(optPB.compute_gradVolume_AN(xk_test), ref_gradvol, rtol=1.e-5)
assert np.allclose(optPB.compute_gradVibration_AN(xk_test, nb_frq=N_TARGET_FREQUENCIES),
                   ref_gradvib, rtol=1.e-5)
