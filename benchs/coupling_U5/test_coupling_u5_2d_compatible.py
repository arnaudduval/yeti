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
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
In this test case, 2 compatible domains are coupled with U5 interface element
(mortar coupling with projection)
Each domain is degree 1 and made of a single element, so 4 conditions can be
expressed (2 coincident control points x 2 dof)
U5 element is defined with degree 0 and elevated to degree 1 in order to have
2 functions for the expression of 2 independant conditions
A test is made to verify that the correct coupling condition are expressed.
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5

TOL = 1.E-10


def test_coupling_u5_2d_compatible(tmp_path):
    """
    Build system matrices and verify couling conditions
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/testU5_2D_compatible')

    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

    # Interface patch is elevated to degree 1
    nb_deg[:, 2] = np.array([1, 0, 0])

    iga_model.refine(nb_ref, nb_deg)

    # Stiffness matrix
    data, row, col, _ = \
        build_stiffmatrix(*iga_model.get_inputs4system_elemStorage())
    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()

    # Coupling matrix
    data, row, col = cplg_matrixU5(*iga_model.get_inputs4cplgmatrixU5(
        output_path=tmp_path))
    cplg_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()

    mat_solve = stiff_side + stiff_side.transpose() + \
        cplg_side + cplg_side.transpose()

    # Coupling conditions for 1st DOF of CPs 2->5 and 4->7 are expressed in
    # lines with indices 16 and 18 of mat_solve (DOFS 3->9 and 7->13)
    sum_lines = mat_solve[16, :] + mat_solve[18, :]

    assert abs(sum_lines[0, 2] + sum_lines[0, 8]) < TOL
    assert abs(sum_lines[0, 6] + sum_lines[0, 12]) < TOL

    # Coupling conditions for 2nd DOF of CPs 2->5 and 4->7 are expressed in
    # lines with indices 17 and 18 of mat_solve (DOFS 4->10 and 8->14)
    sum_lines = mat_solve[17, :] + mat_solve[19, :]

    assert abs(sum_lines[0, 3] + sum_lines[0, 9]) < TOL
    assert abs(sum_lines[0, 7] + sum_lines[0, 13]) < TOL


if __name__ == '__main__':
    test_coupling_u5_2d_compatible('results')
