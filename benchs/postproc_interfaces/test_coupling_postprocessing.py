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


from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp


def test_coupling_postprocessing(tmp_path):
    """
    Compute solution and generate VTU files
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/chamfer')

    # Model refinement
    nb_deg = np.zeros((3, iga_model._nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model._nb_patch), dtype=np.intp)

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

    # Specify number of intergration points
    iga_model._NBPINT[np.where(iga_model._ELT_TYPE == 'U00')] = \
        6**iga_model._dim.min()

    # Matrix Assembly

    ndof = iga_model._nb_dof_free
    idof = iga_model._ind_dof_free[:ndof]-1

    data, row, col, Fb = \
        build_stiffmatrix(*iga_model.get_inputs4system_elemStorage())
    Kside = sp.coo_matrix((data, (row, col)),
                        shape=(iga_model._nb_dof_tot, iga_model._nb_dof_tot),
                        dtype='float64').tocsc()
    Ktot = Kside + Kside.transpose()
    del Kside, data, row, col

    Cdata, Crow, Ccol = cplg_matrix(*iga_model.get_inputs4cplgmatrix())
    Cside = sp.coo_matrix((Cdata, (Crow, Ccol)),
                        shape=(iga_model._nb_dof_tot, iga_model._nb_dof_tot),
                        dtype='float64').tocsc()
    Ctot = Cside + Cside.transpose()
    del Cdata, Crow, Ccol, Cside


    # Monolithic resolution
    K2solve = Ktot[idof, :][:, idof]
    C2solve = Ctot[idof, :][:, idof] * K2solve.max()
    x = sp.linalg.spsolve(K2solve + C2solve,Fb[idof])

    # Global postprocessing
    sol, u = rsol.reconstruction(**iga_model.get_inputs4solution(x))
    pp.generatevtu(*iga_model.get_inputs4postprocVTU(
        'coupling_mono', sol.transpose(), nb_ref=np.array([3, 3, 3]),
        Flag=np.array([True, False, False]),
        output_path=tmp_path))

    # Postprocessing at interfaces
    pp.generate_coupling_vtu(*iga_model.get_inputs4postproc_cplg_vtu(
        'interfaces_10', 10, sol.transpose(), nb_ref=np.array([3, 3]),
        output_path=tmp_path))
    pp.generate_coupling_vtu(*iga_model.get_inputs4postproc_cplg_vtu(
        'interfaces_11', 11, sol.transpose(), nb_ref=np.array([3, 3]),
        output_path=tmp_path))
    pp.generate_coupling_vtu(*iga_model.get_inputs4postproc_cplg_vtu(
        'interfaces_12', 12, sol.transpose(), nb_ref=np.array([3, 3]),
        output_path=tmp_path))


if __name__ == '__main__':
    test_coupling_postprocessing('results')