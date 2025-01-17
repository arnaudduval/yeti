# Copyright 2025 Arnaud Duval

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
Comparaison of U4 and U5 coupling methods in 2D case
Compare computed coupling matrices
"""

import os

import numpy as np
import scipy.sparse as sp

# pylint: disable=no-name-in-module

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.coupling.cplgmatrix import cplg_matrix, cplg_matrixu5

# Models refinement parameters
P = 0
R = 0


def test_compare_u4_u5_2d(tmp_path):
    """
    Build coupling matrices and compare values.
    """

    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model_u4 = IGAparametrization(filename=f'{script_dir}/twoplates_U4')
    iga_model_u5 = IGAparametrization(filename=f'{script_dir}/twoplates_U5')



    # Refine U4 model
    nb_deg = np.zeros((3,iga_model_u4.nb_patch),dtype=np.intp)
    nb_ref = np.zeros((3,iga_model_u4.nb_patch),dtype=np.intp)


    # domains
    nb_deg[:2,:2] = P
    nb_ref[:2, 0] = R
    nb_ref[:2, 1] = R
    # curves
    nb_ref[0,(2,3)] = R
    # lgrge
    nb_deg[0,4] = P
    nb_ref[0,4] = R

    iga_model_u4.refine(nb_ref,nb_deg)
    iga_model_u4.num_integration_points[ np.where(iga_model_u4.elt_type == 'U00') ] = 3


    # refine U5 model

    nb_deg = np.zeros((3,iga_model_u5.nb_patch),dtype=np.intp)
    nb_ref = np.zeros((3,iga_model_u5.nb_patch),dtype=np.intp)

    # domains
    nb_deg[:2,:2] = P
    nb_ref[:2, 0] = R
    nb_ref[:2, 1] = R
    # lgrge
    nb_ref[:,2] = np.array([R,0,0])
    nb_deg[:,2] = np.array([R,0,0])

    iga_model_u5.refine(nb_ref,nb_deg)

    cplg_side = {}

    # Build coupling matrix for U4 model
    data,row,col = cplg_matrix( *iga_model_u4.get_inputs4cplgmatrix() )
    cplg_side['U4'] = sp.coo_matrix((data,(row,col)),
                            shape=(iga_model_u4.nb_dof_tot, iga_model_u4.nb_dof_tot),
                            dtype='float64').tocsc()

    # Build coupling matrix for U5 model
    data, row, col = cplg_matrixu5(*iga_model_u5.get_inputs4cplgmatrixU5(
        integrationOrder=3,
        output_path=tmp_path))

    cplg_side['U5'] = sp.coo_matrix((data,(row,col)),
                            shape=(iga_model_u5.nb_dof_tot, iga_model_u5.nb_dof_tot),
                            dtype='float64').tocsc()

    idof = {}

    # Slice to get components used in linear solver
    idof['U4'] = iga_model_u4.ind_dof_free[:iga_model_u4.nb_dof_free]-1
    idof['U5'] = iga_model_u5.ind_dof_free[:iga_model_u5.nb_dof_free]-1

    error = np.linalg.norm((cplg_side['U4'][idof['U4'],:][:,idof['U4']] - \
                            cplg_side['U5'][idof['U5'],:][:,idof['U5']]).todense()) / \
        np.linalg.norm(cplg_side['U4'][idof['U4'],:][:,idof['U4']].todense())

    print(f"{error = }")

    assert error < 1.E-15


if __name__ == '__main__':
    test_compare_u4_u5_2d('results')
