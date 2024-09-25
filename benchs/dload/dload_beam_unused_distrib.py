# Copyright 2023 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
A beam subjected to distributed triangular pressure
Two distributions are defined, only one is used
No resolution, just build stiffness matrix and force vector
"""

import sys
import time
import copy

import numpy as np
import scipy.sparse as sp

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.preprocessing.igaparametrization import IGAmanip as manip
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp

iga_model = IGAparametrization(filename='beam-dloadF2F4')

# Model refinement parameters (singla patch and 3D for all)
nb_deg = np.zeros((3, 1), dtype=np.intp)
nb_ref = np.zeros((3, 1), dtype=np.intp)

nb_deg[:, 0] = np.array([1, 1, 1])
nb_ref[:, 0] = np.array([2, 2, 5])

iga_model.refine(nb_ref, nb_deg)


def matrix_assembly(model_iga):
    """
    Assemble stiffness matrix and force vector for a given IGA model
    """
    ndof = model_iga._nb_dof_free
    idof = model_iga._ind_dof_free[:ndof]-1
    data, row, col, Fb = build_stiffmatrix(
                        *model_iga.get_inputs4system_elemStorage())
    Kside = sp.coo_matrix((data, (row, col)),
                          shape=(model_iga._nb_dof_tot,
                          model_iga._nb_dof_tot),
                          dtype='float64').tocsc()
    Ktot = Kside + Kside.transpose()

    del Kside, data, row, col

    return copy.copy(Ktot), copy.copy(Fb), copy.copy(idof)

K, F, idof = matrix_assembly(iga_model)


sys.exit(0)
