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
Compute deflection due to 2 triangular pressure fields applied on 2 different
faces of a beam
Compare with the combination of analytic results
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module
# pylint: disable=c-extension-no-member
import yeti_iga.reconstructionSOL as rsol
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix


def test_beam_dlod_combination():
    """
    Tets computed solution
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/inputs_beam/beam-dloadF2F4')

    # Refine modele
    iga_model.refine(np.array([2, 2, 5]), np.array([1, 1, 1]))

    # Matrix assembly
    idof = iga_model.ind_dof_free[:iga_model.nb_dof_free]-1

    data, row, col, rhs = build_stiffmatrix(
                            *iga_model.get_inputs4system_elemStorage())

    stiff_side = sp.coo_matrix((data, (row, col)),
                               shape=(iga_model.nb_dof_tot,
                                      iga_model.nb_dof_tot),
                               dtype='float64').tocsc()
    stiff_tot = stiff_side + stiff_side.transpose()

    x = sp.linalg.spsolve(stiff_tot[idof, :][:, idof], rhs[idof])

    # Solution reconstruction
    sol, _ = rsol.reconstruction(**iga_model.get_inputs4solution(x))

    # Search (interpolating) control point to verify solution
    # Not optimal, we should evaluate the solution at given parameters
    pt_coords = np.array([20., 0., 0.])

    found = False
    for idx in range(np.shape(iga_model.coords)[1]):
        if np.all(iga_model.coords[:, idx] == pt_coords):
            found = True
            break
    assert found

    # Analytic solution (x value is not relevant)
    ref_disp = np.array([np.nan, -20.9592, -300.0611])

    assert np.allclose(sol[idx, 1:], ref_disp[1:], rtol=2.e-2)


if __name__ == '__main__':
    test_beam_dlod_combination()
