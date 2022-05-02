# Copyright 2021 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for 
# more details.
#
# You should have received a copy of the GNU Lesser General Public License along 
# with Yeti. If not, see <https://www.gnu.org/licenses/>

# Python module
import numpy as np

#IGA module
from preprocessing.igaparametrization import IGAparametrization
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip

FILENAME = 'testU5'

modeleIGA = IGAparametrization(filename=FILENAME)

nb_deg = np.zeros((3, modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches": np.array([]), "1": np.array([]), "2": np.array([]), "3": np.array([])}

nb_ref[:, 0] = np.array([3, 3, 3])
nb_deg[:, 0] = np.array([1, 1, 1])

nb_ref[:, 1] = np.array([3, 4, 4])
nb_deg[:, 1] = np.array([1, 1, 1])

nb_ref[:, 2] = np.array([3, 3, 3])
nb_deg[:, 2] = np.array([1, 1, 1])

modeleIGA.refine(nb_ref, nb_deg, additional_knots)

# Output for faces
pp.generate_faces_vtu(**modeleIGA.get_inputs4postproc_faces_vtu(
    FILENAME+'faces', nb_ref=np.array([1, 1, 1])))