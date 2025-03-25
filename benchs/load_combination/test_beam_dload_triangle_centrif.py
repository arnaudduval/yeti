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

import os
import copy

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix


def test_beam_dload_triangle_centrif():
    """
    Compare returned linear systems
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
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
                  'description': 'distributed load on face 2, then \
                   centrifugal force'})
    cases.append({'filename': 'beam-centrif-dloadF4',
                  'description': 'centrifugal force, then distributed load on \
                  face 4'})
    cases.append({'filename': 'beam-dloadF4-centrif',
                  'description': 'distributed load on face 4, then \
                  centrifugal force'})

    # cases with 3 loadings
    cases.append({'filename': 'beam-dloadF2-dloadF4-centrif',
                  'description': 'distributed load on face 2, then \
                  distributed load on face 4, then centrifugal force'})
    cases.append({'filename': 'beam-dloadF4-dloadF2-centrif',
                  'description': 'distributed load on face 4, then \
                  distributed load on face 2, then centrifugal force'})
    cases.append({'filename': 'beam-dloadF2-centrif-dloadF4',
                  'description': 'distributed load on face 2, then \
                  centrifugal force, then distributed load on face 4'})
    cases.append({'filename': 'beam-centrif-dloadF2-dloadF4',
                  'description': 'centrifugal force, then distributed load on \
                  face 2, then distributed load on face 4'})
    cases.append({'filename': 'beam-centrif-dloadF4-dloadF2',
                  'description': 'centrifugal force, then distributed load on \
                  face 4, then distributed load on face 2'})
    cases.append({'filename': 'beam-dloadF4-centrif-dloadF2',
                  'description': 'distributed load on face 4, then \
                  centrifugal force, then distributed load on face 2'})

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
        ndof = model_iga.nb_dof_free
        idof = model_iga.ind_dof_free[:ndof]-1
        data, row, col, rhs = build_stiffmatrix(
                            *model_iga.get_inputs4system_elemStorage())
        stiff_side = sp.coo_matrix((data, (row, col)),
                                   shape=(model_iga.nb_dof_tot,
                                   model_iga.nb_dof_tot),
                                   dtype='float64').tocsc()
        stiff_tot = stiff_side + stiff_side.transpose()

        del stiff_side, data, row, col

        return copy.copy(stiff_tot), copy.copy(rhs), copy.copy(idof)

    for case in cases:
        model_iga = IGAparametrization(
            filename=f'{script_dir}/inputs_beam/{case["filename"]}')
        model_iga.refine(nb_ref, nb_deg)
        _, f_vects[case['filename']], _ = matrix_assembly(model_iga)

    # Tests on force vectors
    # 1 - ensure that F2F4 is the sum of F2 and F4
    assert np.allclose(f_vects['beam-dloadF2F4'],
                       f_vects['beam-dloadF2'] + f_vects['beam-dloadF4'],
                       rtol=1.e-12)

    # 2 - ensure that F4F2 equal to F2F4
    assert np.allclose(f_vects['beam-dloadF2F4'], f_vects['beam-dloadF4F2'],
                       rtol=1.e-12)

    # 3 - ensure that F2 centrif is the sum of F2 and centrif
    assert np.allclose(f_vects['beam-dloadF2-centrif'],
                       f_vects['beam-dloadF2'] + f_vects['beam-centrif'],
                       rtol=1.e-12)

    # 4 - ensure that centrif F2 equal F2 centrif
    assert np.allclose(f_vects['beam-dloadF2-centrif'],
                       f_vects['beam-centrif-dloadF2'],
                       rtol=1.e-12)

    # 5 - ensure that F4 centrif is the sum of F4 and centrif
    assert np.allclose(f_vects['beam-dloadF4-centrif'],
                       f_vects['beam-dloadF4'] + f_vects['beam-centrif'],
                       rtol=1.e-12)

    # 6 - ensure that centrif F4 equal F4 centrif
    assert np.allclose(f_vects['beam-dloadF4-centrif'],
                       f_vects['beam-centrif-dloadF4'],
                       rtol=1.e-12)

    # 7 - ensure that F2 F4 centrif, is the sum of F2, F4 and centrif
    assert np.allclose(f_vects['beam-dloadF2-dloadF4-centrif'],
                       f_vects['beam-dloadF2'] + f_vects['beam-dloadF4'] +
                       f_vects['beam-centrif'],
                       rtol=1.e-12)

    # 8 - ensure the sum works whatever the order of declaration of the 3
    # loadings
    assert np.allclose(f_vects['beam-dloadF2-dloadF4-centrif'],
                       f_vects['beam-dloadF4-dloadF2-centrif'],
                       rtol=1.e-12)
    assert np.allclose(f_vects['beam-dloadF2-dloadF4-centrif'],
                       f_vects['beam-dloadF2-centrif-dloadF4'],
                       rtol=1.e-12)
    assert np.allclose(f_vects['beam-dloadF2-dloadF4-centrif'],
                       f_vects['beam-centrif-dloadF2-dloadF4'],
                       rtol=1.e-12)
    assert np.allclose(f_vects['beam-dloadF2-dloadF4-centrif'],
                       f_vects['beam-centrif-dloadF4-dloadF2'],
                       rtol=1.e-12)
    assert np.allclose(f_vects['beam-dloadF2-dloadF4-centrif'],
                       f_vects['beam-dloadF4-centrif-dloadF2'],
                       rtol=1.e-12)


if __name__ == '__main__':
    test_beam_dload_triangle_centrif()
