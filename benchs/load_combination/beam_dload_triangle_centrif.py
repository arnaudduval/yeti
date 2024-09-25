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
A beam subjected to centrifugal force and distributed triangular pressure
Test to verify if load sum is correct in force vector
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


cases = []

# cases with single loading
cases.append({'filename': 'beam-dloadF2',
              'description': 'a single distributed load on face 2'})
cases.append({'filename': 'beam-dloadF4',
              'description': 'a single distributed load on face 4'})
cases.append({'filename': 'beam-centrif',
              'description': 'centrifugal force only'})

# cases with 2 loadings
cases.append({'filename': 'beam-dloadF2F4',
              'description': '2 distributed loads, face 2 then face 4'})
cases.append({'filename': 'beam-dloadF4F2',
              'description': '2 distributed loads, face 4 then face 2'})
cases.append({'filename': 'beam-centrif-dloadF2',
              'description': 'centrifugal force, then distributed load on \
               face 2'})
cases.append({'filename': 'beam-dloadF2-centrif',
              'description': 'distributed load on face 2, then centrifugal \
               force'})
cases.append({'filename': 'beam-centrif-dloadF4',
              'description': 'centrifugal force, then distributed load on \
               face 4'})
cases.append({'filename': 'beam-dloadF4-centrif',
              'description': 'distributed load on face 4, then centrifugal \
               force'})

# cases with 3 loadings
cases.append({'filename': 'beam-dloadF2-dloadF4-centrif',
              'description': 'distributed load on face 2, then distributed \
               load on face 4, then centrifugal force'})
cases.append({'filename': 'beam-dloadF4-dloadF2-centrif',
              'description': 'distributed load on face 4, then distributed \
               load on face 2, then centrifugal force'})
cases.append({'filename': 'beam-dloadF2-centrif-dloadF4',
              'description': 'distributed load on face 2, then centrifugal \
               force, then distributed load on face 4'})
cases.append({'filename': 'beam-centrif-dloadF2-dloadF4',
              'description': 'centrifugal force, then distributed load on \
               face 2, then distributed load on face 4'})
cases.append({'filename': 'beam-centrif-dloadF4-dloadF2',
              'description': 'centrifugal force, then distributed load on \
               face 4, then distributed load on face 2'})
cases.append({'filename': 'beam-dloadF4-centrif-dloadF2',
              'description': 'distributed load on face 4, then centrifugal \
               force, then distributed load on face 2'})

# Model refinement parameters (single patch and 3D for all)
nb_deg = np.zeros((3, 1), dtype=np.intp)
nb_ref = np.zeros((3, 1), dtype=np.intp)

nb_deg[:, 0] = np.array([0, 0, 1])
nb_ref[:, 0] = np.array([0, 0, 0])


# Compute force vectors
f_vects = {}


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


for case in cases:
    model_iga = IGAparametrization(filename='inputs_beam/'+case['filename'])
    model_iga.refine(nb_ref, nb_deg)
    _, f_vects[case['filename']], _ = matrix_assembly(model_iga)



# Tests on force vectors
# 1 - ensure that F2F4 is the sum of F2 and F4
error = np.linalg.norm(f_vects['beam-dloadF2F4'] -
                       f_vects['beam-dloadF2'] -
                       f_vects['beam-dloadF4'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

# 2 - ensure that F4F2 equal to F2F4
error = np.linalg.norm(f_vects['beam-dloadF2F4'] -
                       f_vects['beam-dloadF4F2'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

# 3 - ensure that F2 centrif is the sum of F2 and centrif
error = np.linalg.norm(f_vects['beam-dloadF2-centrif'] -
                       f_vects['beam-dloadF2'] -
                       f_vects['beam-centrif'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

# 4 - ensure that centrif F2 equal F2 centrif
error = np.linalg.norm(f_vects['beam-dloadF2-centrif'] -
                       f_vects['beam-centrif-dloadF2'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

# 5 - ensure that F4 centrif is the sum of F4 and centrif
error = np.linalg.norm(f_vects['beam-dloadF4-centrif'] -
                       f_vects['beam-dloadF4'] -
                       f_vects['beam-centrif'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

# 6 - ensure that centrif F4 equal F4 centrif
error = np.linalg.norm(f_vects['beam-dloadF4-centrif'] -
                       f_vects['beam-centrif-dloadF4'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

# 7 - ensure that F2 F4 centrif, is the sum of F2, F4 and centrif
error = np.linalg.norm(f_vects['beam-dloadF2-dloadF4-centrif'] -
                       f_vects['beam-dloadF2'] -
                       f_vects['beam-dloadF4'] -
                       f_vects['beam-centrif'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

# 8 - ensure the sum works whatever the order of declaration of the 3 loadings
error = np.linalg.norm(f_vects['beam-dloadF2-dloadF4-centrif'] -
                       f_vects['beam-dloadF4-dloadF2-centrif'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

error = np.linalg.norm(f_vects['beam-dloadF2-dloadF4-centrif'] -
                       f_vects['beam-dloadF2-centrif-dloadF4'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

error = np.linalg.norm(f_vects['beam-dloadF2-dloadF4-centrif'] -
                       f_vects['beam-centrif-dloadF2-dloadF4'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

error = np.linalg.norm(f_vects['beam-dloadF2-dloadF4-centrif'] -
                       f_vects['beam-centrif-dloadF4-dloadF2'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

error = np.linalg.norm(f_vects['beam-dloadF2-dloadF4-centrif'] -
                       f_vects['beam-dloadF4-centrif-dloadF2'])
print(error)
if error > 1.e-12:
    sys.exit(-1)

sys.exit(0)
