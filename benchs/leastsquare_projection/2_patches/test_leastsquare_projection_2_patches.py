# Copyright 2021 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
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
A single domain with 2 patches
Test generation of a vtu file with with computation of values at nodes points
using least square projection of values at Gauss points.
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=c-extension-no-member
# pylint: disable=no-name-in-module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp


def test_leastsquare_projection_2_patches(tmp_path):
    """
    Compute solution and generate VTU file
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/plateWithRoundHole2P')

    # odel refinement
    iga_model.refine(np.array([[3, 3, 0], [3, 3, 0]]).T,
                     np.array([[0, 0, 0], [0, 0, 0]]).T)

    # Matrix assembly
    ndof = iga_model.nb_dof_free
    idof = iga_model.ind_dof_free[:ndof]-1

    data, row, col, rhs = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())
    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    stiff_tot = stiff_side + stiff_side.transpose()

    # Resolution
    x = sp.linalg.spsolve(stiff_tot[idof, :][:, idof], rhs[idof])

    # Standard postprocessing
    sol, _ = rsol.reconstruction(**iga_model.get_inputs4solution(x))
    # Stress and strain post processing is deactivated due to a singular point
    # leading to NaN values
    pp.generatevtu(*iga_model.get_inputs4postprocVTU(
        'linear_analysis', sol.transpose(), nb_ref=np.zeros(3),
        Flag=np.array([True, False, False]),
        output_path=tmp_path))

    # Postprocessing with least square projection
    # Compute gram matrix
    data, row, col = pp.build_cgrammatrix(*iga_model.get_inputs4grammat())
    # RHS vector
    rhs = pp.compute_svars_solid_rhs(
        *iga_model.get_inputs4svarsrhs(sol.transpose()))

    gram_side = sp.coo_matrix((data, (row, col)),
                              shape=(max(col)+1, max(col)+1),
                              dtype='float64').tocsc()
    gram_tot = gram_side + gram_side.transpose()

    projected_svars = sp.linalg.spsolve(gram_tot, rhs.transpose())

    pp.generate_proj_vtu(*iga_model.get_inputs4proj_vtu('least_square',
                         sol.transpose(),
                         projected_svars.transpose(),
                         nb_ref=np.array([3, 3, 0]),
                         output_path=tmp_path))


if __name__ == '__main__':
    test_leastsquare_projection_2_patches('results')
