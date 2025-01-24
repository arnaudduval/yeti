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
3 coupled domain with U4 coupling
Test VTU output the coupling conditions on the interfaces
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=c-extension-no-member
# pylint: disable=no-name-in-module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp


def refine_model(iga_model):
    """
    Refine patches of an IGA model
    """

    # Model refinement
    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

    # domain 1 : chamfer
    nb_deg[:, 0] = [0, 0, 0]
    nb_ref[:, 0] = [2, 3, 2]
    # domain 2 : tube
    nb_deg[:, 1] = [1, 0, 1]
    nb_ref[:, 1] = [1, 3, 2]
    # domain 3 : base
    nb_deg[:, 2] = [1, 0, 1]
    nb_ref[:, 2] = [3, 3, 2]

    additional_knots = {"patches": np.array([1]),
                        "1": np.array([]),
                        "2": np.array([]),
                        "3": np.array([0.25])
                        }

    # Interfaces
    nb_deg[:2, (2, 3, 5, 6, 7, 8)] = 0
    nb_ref[:2, (3, 4)] = ([1, 3], [1, 3])
    nb_ref[:2, (5, 6)] = ([1, 3], [1, 3])
    nb_ref[:2, (7, 8)] = ([1, 3], [1, 3])

    # Lagrange multiplier
    nb_deg[:2, (9, 10, 11)] = 1
    nb_ref[:2, 9] = [1, 3]
    nb_ref[:2, 10] = [1, 3]
    nb_ref[:2, 11] = [1, 3]

    iga_model.refine(nb_ref, nb_deg, additional_knots)

    return iga_model


def test_coupling_postprocessing(tmp_path):
    """
    Compute solution and generate VTU files
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/chamfer')

    iga_model = refine_model(iga_model=iga_model)

    # Specify number of intergration points
    iga_model.num_integration_points[np.where(iga_model.elt_type == 'U00')] = \
        6**iga_model.dim.min()

    # Matrix Assembly
    idof = iga_model.ind_dof_free[:iga_model.nb_dof_free] - 1

    data, row, col, rhs = \
        build_stiffmatrix(*iga_model.get_inputs4system_elemStorage())
    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    stiff_tot = stiff_side + stiff_side.T

    data, row, col = cplg_matrix(*iga_model.get_inputs4cplgmatrix())
    cplg_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    cplg_tot = cplg_side + cplg_side.T

    # Monolithic resolution
    x = sp.linalg.spsolve(stiff_tot[idof, :][:, idof] +
                          cplg_tot[idof, :][:, idof] * stiff_tot.max(),
                          rhs[idof])

    # Global postprocessing
    sol, _ = rsol.reconstruction(**iga_model.get_inputs4solution(x))
    pp.generatevtu(*iga_model.get_inputs4postprocVTU(
        'coupling_mono', sol.T, nb_ref=np.array([3, 3, 3]),
        Flag=np.array([True, False, False]),
        output_path=tmp_path))

    # Postprocessing at interfaces
    pp.generate_coupling_vtu(*iga_model.get_inputs4postproc_cplg_vtu(
        'interfaces_10', 10, sol.T, nb_ref=np.array([3, 3]),
        output_path=tmp_path))
    return
    pp.generate_coupling_vtu(*iga_model.get_inputs4postproc_cplg_vtu(
        'interfaces_11', 11, sol.T, nb_ref=np.array([3, 3]),
        output_path=tmp_path))
    pp.generate_coupling_vtu(*iga_model.get_inputs4postproc_cplg_vtu(
        'interfaces_12', 12, sol.T, nb_ref=np.array([3, 3]),
        output_path=tmp_path))


if __name__ == '__main__':
    test_coupling_postprocessing('results')
