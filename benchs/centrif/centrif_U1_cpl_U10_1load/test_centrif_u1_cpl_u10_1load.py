# Copyright 2020 Arnaud Duval

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
Compute deflexion due to distributed pressure field on a beam
Compare with analytic result
1 U1 patch coupled with 1 U10 patch, 1 single load defined on both patchs
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=c-extension-no-member
# pylint: disable=no-name-in-module

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5
import yeti_iga.reconstructionSOL as rsol

REF_PT_COORDS = np.array([12.0, 1.5, 0.0])
REF_ABAQUS_SOLUTION = np.array([1.60785E-6, -1.94110E-8, 3.15129E-7])


def test_centrif_u1_cpl_u10_1load(tmp_path):
    """
    Build system matrices, solce the problem and compare with reference
    solution.
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/centrif_U1_cpl_U10_1load')

    # Initialise refinement arrays
    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    additional_knots = {
        "patches": np.array([]),
        "1": np.array([]),
        "2": np.array([]),
        "3": np.array([])
    }

    # U1 patch
    nb_ref[:, 0] = np.array([2, 2, 2])
    nb_deg[:, 0] = np.array([1, 1, 1])
    # Hull
    nb_ref[:, 1] = np.array([0, 0, 0])
    nb_deg[:, 1] = np.array([0, 0, 0])
    # U10 patch
    nb_ref[:, 2] = np.array([2, 2, 2])
    nb_deg[:, 2] = np.array([1, 1, 1])
    # coupling
    nb_ref[:, 3] = np.array([2, 2, 0])
    nb_deg[:, 3] = np.array([1, 1, 0])

    # Refine
    iga_model.refine(nb_ref, nb_deg, additional_knots)
    idof = iga_model.ind_dof_free[:iga_model.nb_dof_free] - 1

    # Build stiffness matrix
    data, row, col, rhs = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())
    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64'
    ).tocsc()

    # Build coupling matrix
    data, row, col = cplg_matrixU5(*iga_model.get_inputs4cplgmatrixU5(
        output_path=tmp_path
    ))
    cplg_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()

    x = sp.linalg.spsolve(
        (stiff_side + stiff_side.transpose())[idof, :][:, idof] +
        (cplg_side + cplg_side.transpose())[idof, :][:, idof] *
        stiff_side.max(),
        rhs[idof])
    sol, _ = rsol.reconstruction(**iga_model.get_inputs4solution(x))

    # Get solution at interpolating control points and compare with Abaqus
    # reference result
    # Tolerance = max 2% on quasi radial component
    found = False
    for idx in range(np.shape(iga_model.coords)[1]):
        if np.all(iga_model.coords[:, idx] == REF_PT_COORDS):
            found = True
            break

    assert found
    assert np.abs(
        (sol[idx, :]-REF_ABAQUS_SOLUTION)/REF_ABAQUS_SOLUTION
    )[0] < 0.02


if __name__ == '__main__':
    test_centrif_u1_cpl_u10_1load('results')
