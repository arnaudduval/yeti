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

#!/usr/bin/env python
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
nb_deg[:, 0] = np.array([1, 1, 1])
additional_knots = {"patches": np.array([0]),
                    "1": np.array([]), "2": np.array([0.5, 0.5]), "3": np.array([])}


# Initial refinement (none)
modeleIGA.refine(nb_ref, nb_deg, additional_knots)




