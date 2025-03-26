# Copyright 2020 Arnaud Duval

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
Compute eigenfrequencies of a plate modeled with 3D solid elements
Results are compared with numerical reference
"""

import os
import numpy as np
import scipy.sparse as sp

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
from yeti_iga.massmtrx import build_cmassmatrix


def test_platevib():
    """
    Compute eigenfrequencies of a plate modeled with 3D solid elements
    Results are compared with numerical reference
    """

    # Creation of the IGA object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/plateVolume')

    # Refine model
    iga_model.refine(nb_refinementByDirection=np.array([[2, 2, 1]]).T,
                     nb_degreeElevationByDirection=np.array([[1, 1, 1]]).T)

    # Build stiffness matrix
    data, row, col, _ = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())
    k_side = sp.coo_matrix((data, (row, col)),
                           shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
                           dtype='float64').tocsc()
    ndof = iga_model.nb_dof_free
    idof = iga_model.ind_dof_free[:ndof]-1
    k_solve = (k_side + k_side.transpose())[idof, :][:, idof]
    del k_side

    # Build mass matrix
    data, row, col = build_cmassmatrix(*iga_model.get_inputs4massmat())
    m_side = sp.coo_matrix((data, (row, col)),
                           shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
                           dtype='float64').tocsc()
    m_solve = (m_side + m_side.transpose())[idof, :][:, idof]
    del m_side, data, row, col

    # Compute eigenvalues and eigenfrequencies
    nb_frq = 5
    vals, _ = sp.linalg.eigsh(k_solve, k=nb_frq, M=m_solve, sigma=0.)
    frq = np.sqrt(vals[:])/2./np.pi

    # Reference eigenfrequencies
    ref_frq = np.array([0.2548, 0.6217, 1.9546, 2.3870, 2.8585])

    assert np.allclose(frq, ref_frq, rtol=1.e-3)


if __name__ == '__main__':
    test_platevib()
