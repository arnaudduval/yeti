# Copyright 2025 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>

# !/usr/bin/env python
# -*- coding: utf-8 -*-


"""
Compute deflexion due to volume force
Compare with reference Abaqus result
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module

# pylint: disable=c-extension-no-member
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp

REF_SOLUTION = [1.49192, 8.68649]


def test_volume_force_2d():
    """
    Compute solution and compare with reference value.
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/volume_force_2D')

    # Refine modele
    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

    nb_deg[0, 0] = 1
    nb_deg[1, 0] = 1

    nb_ref[0, 0] = 4
    nb_ref[1, 0] = 3

    iga_model.refine(nb_ref, nb_deg)

    # Matrix assembly
    idof = iga_model.ind_dof_free[:iga_model.nb_dof_free]-1

    data, row, col, rhs = build_stiffmatrix(
                            *iga_model.get_inputs4system_elemStorage())

    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()

    x = sp.linalg.spsolve(
        (stiff_side + stiff_side.transpose())[idof, :][:, idof],
        rhs[idof])

    # Solution reconstruction
    sol, _ = rsol.reconstruction(**iga_model.get_inputs4solution(x))

    disp = pp.evaldisp(*iga_model.get_inputs4evaldisp(
        sol.transpose(), np.array([1., 0., 0.]), numpatch=1))[:2]

    error = np.linalg.norm(disp-REF_SOLUTION)/np.linalg.norm(REF_SOLUTION)

    assert error < 2.e-3


if __name__ == '__main__':
    test_volume_force_2d()
