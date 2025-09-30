# Copyright 2025 Arnaud Duval

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
Build stiffness matrix and force vector for a beam with distributed load
Matrix is build with sequential routine and with OpenMP routine.
Obtained stiffness matrix and force vector are compared
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module
# pylint: disable=c-extension-no-member
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
from yeti_iga.stiffmtrx_elemstorage_omp import sys_linmat_lindef_static_omp \
    as build_stiffmatrix_omp


def test_build_stoffmatrix_omp():
    """
    Test built matrix and vector
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/beam-dist')

    # Refine modele
    iga_model.refine(np.array([2, 2, 2]), np.array([2, 2, 2]))

    # Matrix assembly
    data, row, col, rhs_seq = build_stiffmatrix(
                            *iga_model.get_inputs4system_elemStorage())

    stiff_side_seq = sp.coo_matrix((data, (row, col)),
                                   shape=(iga_model.nb_dof_tot,
                                          iga_model.nb_dof_tot),
                                   dtype='float64').tocsc()

    data, row, col, rhs_omp = build_stiffmatrix_omp(
        **iga_model.get_inputs4system_elemStorage_OMP(num_threads=2))

    stiff_side_omp = sp.coo_matrix((data, (row, col)),
                                   shape=(iga_model.nb_dof_tot,
                                          iga_model.nb_dof_tot),
                                   dtype='float64').tocsc()

    assert sp.linalg.norm(stiff_side_seq - stiff_side_omp)  \
        / sp.linalg.norm(stiff_side_seq) < 1.e-15
    assert np.linalg.norm(rhs_seq - rhs_omp)/np.linalg.norm(rhs_seq) < 1.e-15


if __name__ == '__main__':
    test_build_stoffmatrix_omp()
