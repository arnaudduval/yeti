# Copyright 2021 Marie Guerder
# Copyright 2024 Arnaud Duval

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
A non linear cube embedded in a linear hull
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix


def test_non_linear_embedded_solid():
    """
    Compute solution and compare to a reference solution
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/test_embd_cube')

    # Set arguments for model refinement
    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    additional_knots = {"patches": np.array([]),
                        "1": np.array([]),
                        "2": np.array([]),
                        "3": np.array([])}

    nb_ref[:, 1] = [1, 1, 1]    # Knot insertion, [xi, eta, zeta]

    # Refine model
    iga_model.refine(nb_ref, nb_deg, additional_knots)

    # Matrix assembly
    ndof = iga_model.nb_dof_free
    idof = iga_model.ind_dof_free[:ndof] - 1

    data, row, col, rhs = build_stiffmatrix(
                            *iga_model.get_inputs4system_elemStorage())
    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    stiff_tot = stiff_side + stiff_side.T

    # Problem resolution with direct solver
    x = sp.linalg.spsolve(stiff_tot[idof, :][:, idof], rhs[idof])

    ref_sol = np.loadtxt(f'{script_dir}/ref_sol.txt')

    assert np.allclose(x, ref_sol, rtol=1.e-6)


if __name__ == '__main__':
    test_non_linear_embedded_solid()
