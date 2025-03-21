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
Generate VTU output using Bezier en assembly of 2 patchs with Mortar
coupling
No test is made on the content of output files.
WARNING : THIS CASE CREATES A SINGULAR MATRIX
"""

import os

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp


def test_vtu_output_2patchs_coupling():
    """
    Solve analysis and write VTU result file
    """
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/inputs/BentPipe_2patchs_coupling_U5')

    # Model refinement
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)


    nb_ref[:, 0] = np.array([0, 0, 0])
    nb_ref[:, 1] = np.array([0, 0, 0])

    nb_deg[:, 0] = np.array([0, 1, 0])
    nb_deg[:, 1] = np.array([0, 1, 0])

    # Coupling patches
    nb_ref[:, 2] = np.array([0, 0, 0])
    nb_ref[:, 3] = np.array([0, 0, 0])

    nb_deg[:, 2] = np.array([1, 1, 0])
    nb_deg[:, 3] = np.array([1, 1, 0])

    print(iga_model.nb_patch)
    iga_model.refine(nb_ref, nb_deg)

    # Matrix assembly
    ndof = iga_model._nb_dof_free
    idof = iga_model._ind_dof_free[:ndof]-1

    print("Build stiffmatrix")

    data, row, col, Fb = build_stiffmatrix(
                            *iga_model.get_inputs4system_elemStorage())

    Kside = sp.coo_matrix((data, (row, col)),
                        shape=(iga_model._nb_dof_tot, iga_model._nb_dof_tot),
                        dtype='float64').tocsc()
    Ktot = Kside + Kside.transpose()
    del Kside, data, row, col


    print("Build coupling matrix")
    Cdata, Crow, Ccol = cplg_matrixU5( *iga_model.get_inputs4cplgmatrixU5() )
    Cside = sp.coo_matrix((Cdata,(Crow,Ccol)), shape=(iga_model._nb_dof_tot,iga_model._nb_dof_tot),
                        dtype='float64').tocsc()
    Ctot  = Cside + Cside.transpose()
    del Cdata,Crow,Ccol,Cside


    # Resolution
    K2solve = Ktot[idof,:][:,idof]
    C2solve = Ctot[idof,:][:,idof] * K2solve.max()
    plt.spy(K2solve + C2solve, markersize=0.1)
    # plt.show()

    x = sp.linalg.spsolve(K2solve + C2solve, Fb[idof], use_umfpack=True)

    # Solution reconstruction
    SOL, u = rsol.reconstruction(**iga_model.get_inputs4solution(x))

    pp.generatevtu(*iga_model.get_inputs4postprocVTU(
        'couplingU5',SOL.transpose(),nb_ref=np.array([3,3,3]),
        Flag=np.array([True,True,False])))

    exit()

    # VTK output using Bezier elements
    for i_patch in range(iga_model._nb_patch):
        pp.generate_vtu_bezier(**iga_model.get_inputs4postproc_bezier(
            i_patch+1,
            f'BentPipeBezier_2patchs_P{i_patch+1}',
            SOL.transpose(),
            ))

if __name__ == '__main__':
    test_vtu_output_2patchs_coupling()