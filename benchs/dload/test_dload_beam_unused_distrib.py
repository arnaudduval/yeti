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
A beam subjected to distributed triangular pressure
Two distributions are defined, only one is used
No resolution, just build stiffness matrix and force vector
"""

import os
import copy

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix


def test_dload_beam_unused_distrib():
    """
    Test if stiffness matrix and force vector are build (no value check)
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/beam-dloadF2F4')

    # Refine modele
    iga_model.refine(np.array([2, 2, 5]), np.array([1, 1, 1]))

    def matrix_assembly(model):
        """
        Assemble stiffness matrix and force vector for a given IGA model
        """
        ndof = model.nb_dof_free
        idof = model.ind_dof_free[:ndof]-1
        data, row, col, rhs = build_stiffmatrix(
            *model.get_inputs4system_elemStorage())
        stiff_side = sp.coo_matrix((data, (row, col)),
                                  shape=(model.nb_dof_tot,
                                  model.nb_dof_tot),
                                  dtype='float64').tocsc()
        stiff_tot = stiff_side + stiff_side.transpose()

        del stiff_side, data, row, col

        return copy.copy(stiff_tot), copy.copy(rhs), copy.copy(idof)

    _, __, ___ = matrix_assembly(iga_model)


if __name__ == '__main__':
    test_dload_beam_unused_distrib()
