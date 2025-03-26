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
In this test, 2 incompatible meshes are coupled with U5 inreface element
Solution at a given point is compared with a reference value computed with
Abaqus.
"""

import os

import numpy as np
from numpy import intersect1d
import scipy.sparse as sp

# pylint: disable=c-extension-no-member
# pylint: disable=no-name-in-module

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5
import yeti_iga.reconstructionSOL as rsol


# Reference solution computed with Abaqus
REF = -1.136E-2
TOL = 0.03


def test_coupling_u5(tmp_path):
    """
    Build system matrices, solve the problem and compare with reference
    solution.
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/testU5')

    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

    nb_ref[:, 0] = np.array([2, 2, 2])
    nb_deg[:, 0] = np.array([1, 1, 1])

    nb_ref[:, 1] = np.array([2, 3, 3])
    nb_deg[:, 1] = np.array([1, 1, 1])

    nb_ref[:, 2] = np.array([2, 2, 0])
    nb_deg[:, 2] = np.array([1, 1, 0])

    iga_model.refine(nb_ref, nb_deg)

    idof = iga_model.ind_dof_free[:iga_model.nb_dof_free]-1

    # Stiffness matrix
    data, row, col, rhs = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())
    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()

    # Coupling matrix
    data, row, col = cplg_matrixU5(*iga_model.get_inputs4cplgmatrixU5(
        output_path=tmp_path
    ))
    cplg_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()

    # Full monolithic matrix
    mat_solve = (stiff_side+stiff_side.transpose())[idof, :][:, idof] + \
        (cplg_side + cplg_side.transpose())[idof, :][:, idof] \
        * stiff_side.max()

    x = sp.linalg.spsolve(mat_solve, rhs[idof])

    # Postprocessing
    sol, _ = rsol.reconstruction(**iga_model.get_inputs4solution(x))

    # Verify result
    assert np.allclose(
        sol[search_point(iga_model, np.array([6., 0., 2.]))[0], 2],
        REF,
        rtol=TOL)


def search_point(iga_model, xtest):
    """
    search index of DOF for given point coordinates
    """

    search = np.arange(np.shape(iga_model.coords)[1])
    for i in range(3):
        search = intersect1d(search,
                             np.where(iga_model.coords[i, :] == xtest[i]))
    return search


if __name__ == '__main__':
    test_coupling_u5('results')
