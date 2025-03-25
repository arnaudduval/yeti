# Copyright 2023 Arnaud Duval

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
A single shell patch clamped at corner
Shell is divided in 9 elements (1 element read from input files and division
made in this script) and only central
element is loaded with pressure.
Normal displacement at center is compared with reference result computed in
Abaqus.
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


# Reference normal displacement at center computed with Abaqus
# (S8R elements, mesh size 0.1)
REF_ABA = -2.140E-4


def test_shell_load_1elt():
    """
    Test computed solution
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/1_element_degree_3')

    # Model refinement
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    additional_knots = {'patches': np.array([0]),
                        '1': np.array([0.33, 0.67]),
                        '2': np.array([0.33, 0.67]),
                        '3': np.array([])}

    iga_model.refine(nb_ref, np.array([0, 0, 0]),
                     additional_knots=additional_knots)

    # Add loading on central element
    iga_model._indDLoad = np.array([[5]])
    iga_model._JDLType = np.array([60])
    iga_model._ADLMAG = np.array([100.])
    iga_model._load_target_nbelem = np.array([1])
    iga_model._nb_load = 1
    iga_model._indDLoad_flat = np.array([], dtype=np.intp)
    for load in iga_model._indDLoad:
        iga_model._indDLoad_flat = np.hstack((iga_model._indDLoad_flat, load))

    nb_ref[:, 0] = np.array([2, 2, 0])
    iga_model.refine(nb_ref, np.array([0, 0, 0]))

    # Matrix assembly
    idof = iga_model.ind_dof_free[:iga_model.nb_dof_free]-1

    data, row, col, rhs = build_stiffmatrix(
                            *iga_model.get_inputs4system_elemStorage())
    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()
    stiff_tot = stiff_side + stiff_side.transpose()

    # Resolution
    x = sp.linalg.spsolve(stiff_tot[idof, :][:, idof], rhs[idof])

    # Solution reconstruction
    sol, _ = rsol.reconstruction(**iga_model.get_inputs4solution(x))

    disp_center = pp.evaldisp(*iga_model.get_inputs4evaldisp(
        sol.transpose(), np.array([0.5, 0.5, 0.]), numpatch=1))

    print('yeti solution : ', disp_center[2])
    print('reference solution : ', REF_ABA)
    print('error : ', np.abs((REF_ABA-disp_center[2])/REF_ABA))

    assert np.allclose(disp_center[2], REF_ABA, rtol=4.e-2)


if __name__ == '__main__':
    test_shell_load_1elt()
