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
An full non-linear (immersed and hull) embedded beam
Loading: 2 distributed load (artifical) on faces 1 and 4
+ a centrifugal body force
Test to verify if load sum is correct in force vector
"""

import sys
import random
import copy

import numpy as np
import scipy.sparse as sp

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix

if __name__ == '__main__':
    cases = []

    # Cases with single loading
    cases.append({'filename': 'embeddedbeam-dloadF2',
                  'description': 'a single distributed load on face 2'})
    cases.append({'filename': 'embeddedbeam-dloadF4',
                  'description': 'a single distributed load on face 4'})
    cases.append({'filename': 'embeddedbeam-centrif',
                  'description': 'centrifigal force only'})

    # Cases with 2 loadings
    cases.append({'filename': 'embeddedbeam-dloadF2F4',
                  'description': '2 distributed loads, face 2 then face 4'})
    cases.append({'filename': 'embeddedbeam-dloadF4F2',
                  'description': '2 distributed loads, face 4 then face 2'})
    cases.append({'filename': 'embeddedbeam-centrif-dloadF2',
                  'description': 'centrifugal force, then distributed load on \
                   face 2'})
    cases.append({'filename': 'embeddedbeam-dloadF2-centrif',
                  'description': 'distributed load on face 2, then \
                   centrifugal force'})
    cases.append({'filename': 'embeddedbeam-centrif-dloadF4',
                  'description': 'centrifugal force, then distributed load on \
                   face 4'})
    cases.append({'filename': 'embeddedbeam-dloadF4-centrif',
                  'description': 'distributed load on face 4, then \
                   centrifugal force'})

    # cases with 3 loadings
    cases.append({'filename': 'embeddedbeam-dloadF2-dloadF4-centrif',
                  'description': 'distributed load on face 2, then \
                   distributed load on face 4, then centrifugal force'})
    cases.append({'filename': 'embeddedbeam-dloadF4-dloadF2-centrif',
                  'description': 'distributed load on face 4, then \
                   distributed load on face 2, then centrifugal force'})
    cases.append({'filename': 'embeddedbeam-dloadF2-centrif-dloadF4',
                  'description': 'distributed load on face 2, then \
                   centrifugal force, then distributed load on face 4'})
    cases.append({'filename': 'embeddedbeam-centrif-dloadF2-dloadF4',
                  'description': 'centrifugal force, then distributed load on \
                   face 2, then distributed load on face 4'})
    cases.append({'filename': 'embeddedbeam-centrif-dloadF4-dloadF2',
                  'description': 'embeddedcentrifugal force, then distributed \
                   load on face 4, then distributed load on face 2'})
    cases.append({'filename': 'embeddedbeam-dloadF4-centrif-dloadF2',
                  'description': 'distributed load on face 4, then \
                   centrifugal force, then distributed load on face 2'})

    # Model refinement parameters (2 patches and 3D for all)
    nb_deg = np.zeros((3, 2), dtype=np.intp)
    nb_ref = np.zeros((3, 2), dtype=np.intp)

    # Hull
    nb_deg[:, 0] = np.array([0, 0, 0], dtype=int)
    nb_ref[:, 0] = np.array([2, 1, 1], dtype=int)

    # Embedded entity
    nb_deg[:, 1] = np.array([0, 0, 0], dtype=int)
    nb_ref[:, 1] = np.array([1, 1, 2], dtype=int)

    # Force vectors
    f_vects = {}

    def matrix_assembly(model):
        """
        Assemble stiffness matrix and force vector for a given IGA model
        """
        ndof = model.nb_dof_free
        idof = model.ind_dof_free[:ndof]-1
        data, row, col, fvect = build_stiffmatrix(
                        *model.get_inputs4system_elemStorage())
        kside = sp.coo_matrix((data, (row, col)),
                              shape=(iga_model.nb_dof_tot, model.nb_dof_tot),
                              dtype='float64').tocsc()
        ktot = kside + kside.transpose()

        del kside, data, row, col

        return copy.copy(ktot), copy.copy(fvect), copy.copy(idof)

    # Movable control points coordinates to introduce non linearity
    # Array : index of cp, index of coordinate
    movable = np.array([[1, 1], [3, 2], [4, 1], [4, 2], [5, 2], [7, 1],
                        [9, 0], [10, 1], [10, 0], [11, 0], [12, 2], [12, 0],
                        [13, 1], [13, 2], [13, 0], [14, 2], [14, 0], [15, 0],
                        [16, 1], [16, 0], [17, 0],
                        [19, 1], [21, 2], [22, 1], [22, 2], [23, 2], [25, 1],
                        # End of hull
                        [28, 0], [30, 1], [31, 0], [31, 1], [32, 1], [34, 0],
                        [36, 2], [37, 0], [37, 2], [38, 2], [39, 1], [39, 2],
                        [40, 0], [40, 1], [40, 2], [41, 1], [41, 2], [42, 2],
                        [43, 0], [43, 2], [44, 2],
                        [46, 0], [48, 1], [49, 0], [49, 1], [50, 1], [52, 0]])

    for case in cases:
        iga_model = IGAparametrization(filename='inputs_embeddedbeam/' +
                                       case['filename'])

        # Move control points
        random.seed(3)
        # Hull size
        hull_size = np.array([20., 3., 1.])
        # Move hull CPs
        for move in movable[:movable.shape[0]//2]:
            # CP at index move[0] is moved in move[1] direction
            iga_model.coords[move[1], move[0]] = \
                hull_size[move[1]]*(0.1 + 0.8*random.random())

        # Move embedded entity CPs
        for move in movable[movable.shape[0]//2:]:
            # CP at index move[0] is moved in move[1] direction
            iga_model.coords[move[1], move[0]] = 0.1 + 0.8*random.random()

        # Refine model
        iga_model.refine(nb_ref, nb_deg)

        _, f_vects[case['filename']], _ = matrix_assembly(iga_model)

    # Tests on force vectors
    # 1 - ensure that F2F4 is the sum of F2 and F4
    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2F4'] -
                           f_vects['embeddedbeam-dloadF2'] -
                           f_vects['embeddedbeam-dloadF4'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    # 2 - ensure that F4F2 equal to F2F4
    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2F4'] -
                           f_vects['embeddedbeam-dloadF4F2'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    # 3 - ensure that F2 centrif is the sum of F2 and centrif
    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2-centrif'] -
                           f_vects['embeddedbeam-dloadF2'] -
                           f_vects['embeddedbeam-centrif'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    # 4 - ensure that centrif F2 equal F2 centrif
    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2-centrif'] -
                           f_vects['embeddedbeam-centrif-dloadF2'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    # 5 - ensure that F4 centrif is the sum of F4 and centrif
    error = np.linalg.norm(f_vects['embeddedbeam-dloadF4-centrif'] -
                           f_vects['embeddedbeam-dloadF4'] -
                           f_vects['embeddedbeam-centrif'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    # 6 - ensure that centrif F4 equal F4 centrif
    error = np.linalg.norm(f_vects['embeddedbeam-dloadF4-centrif'] -
                           f_vects['embeddedbeam-centrif-dloadF4'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    # 7 - ensure that F2 F4 centrif, is the sum of F2, F4 and centrif
    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2-dloadF4-centrif'] -
                           f_vects['embeddedbeam-dloadF2'] -
                           f_vects['embeddedbeam-dloadF4'] -
                           f_vects['embeddedbeam-centrif'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    # 8 - ensure the sum works whatever the order of declaration
    # of the 3 loadings
    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2-dloadF4-centrif'] -
                           f_vects['embeddedbeam-dloadF4-dloadF2-centrif'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2-dloadF4-centrif'] -
                           f_vects['embeddedbeam-dloadF2-centrif-dloadF4'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2-dloadF4-centrif'] -
                           f_vects['embeddedbeam-centrif-dloadF2-dloadF4'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2-dloadF4-centrif'] -
                           f_vects['embeddedbeam-centrif-dloadF4-dloadF2'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    error = np.linalg.norm(f_vects['embeddedbeam-dloadF2-dloadF4-centrif'] -
                           f_vects['embeddedbeam-dloadF4-centrif-dloadF2'])
    print(error)
    if error > 1.e-12:
        sys.exit(-1)

    sys.exit(0)
