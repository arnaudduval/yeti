# Copyright 2021 Arnaud Duval

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
Coupled domains woth monolithic resolution
Test generation of a vtu file with with computation of values at nodes points
using least square projection of values at Gauss points.
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


def refine_model(iga_model, p, r):
    """
    Refine patch of an IGA model with 2 domains couples

    Parameters
    ----------
    iga_model: IGAparametrization object
        Model to refine
    p : int
        Reference degree elevation
    r : int
        Reference nomber of subdivisions

    Returns
    -------
    model: IGAparametrization object
        refined model
    """

    # Model refinement
    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

    # domain 1 (initially at degree 2)
    nb_deg[:, 0] = [p - 2, p - 2, p - 2]
    nb_ref[:, 0] = [r, r, 0]
    # domain 2 (initially at degree 2)
    nb_deg[:, 1] = [p - 2, p - 2, p - 2]
    nb_ref[:, 1] = [r+1, r+1, 0]
    additional_knots = {"patches": np.array([1]),
                        "1": np.array([]),
                        "2": np.array([0.3]),
                        "3": np.array([])}
    # interface
    nb_ref[:1, 2] = [r+2]
    nb_ref[:1, 3] = [r+2]
    # Lagrange field (initially at degree 0)
    nb_deg[:1, 4] = p - 1
    nb_ref[:1, 4] = r

    iga_model.refine(nb_ref, nb_deg, additional_knots)

    return iga_model


def test_leastsquare_projection_2_patchescplg(tmp_path):
    """
    Compute solution and generate VTU file
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/plateWithRoundHoleCPLG')

    iga_model = refine_model(iga_model, p=2, r=2)

    iga_model.num_integration_points[np.where(iga_model.elt_type == 'U00')] = \
        6**iga_model.dim.min()

    # Matrix assembly
    idof = iga_model.ind_dof_free[:iga_model.nb_dof_free]-1

    # Stiffness matrix
    data, row, col, rhs = build_stiffmatrix(
        *iga_model.get_inputs4system_elemStorage())
    stiff_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()

    # Coupling matrix
    data, row, col = cplg_matrix(*iga_model.get_inputs4cplgmatrix())
    cplg_side = sp.coo_matrix(
        (data, (row, col)),
        shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
        dtype='float64').tocsc()

    # Monolithic solving
    x = sp.linalg.spsolve(
        (stiff_side + stiff_side.transpose())[idof, :][:, idof] +
        (cplg_side + cplg_side.transpose())[idof, :][:, idof] *
        stiff_side.max(),
        rhs[idof])

    # Standard Postprocessing
    sol, _ = rsol.reconstruction(**iga_model.get_inputs4solution(x))
    pp.generatevtu(*iga_model.get_inputs4postprocVTU(
            'coupling_mono', sol.transpose(), nb_ref=np.array([3, 3, 3]),
            Flag=np.array([True, False, False]),
            output_path=tmp_path))

    # Least square projection post processing
    # Gram matrix
    data, row, col = pp.build_cgrammatrix(*iga_model.get_inputs4grammat())
    # RHS vector
    rhs = pp.compute_svars_solid_rhs(
        *iga_model.get_inputs4svarsrhs(sol.transpose()))

    # get max node index for U1 elements
    # WARNING : Assumption : U1 elements are given first
    max_ind = 0
    for i in np.where(iga_model.elt_type == 'U1')[0]:
        max_ind = max(max_ind, np.max(iga_model.ien[i]))

    # Add identity for Gram matrix coefficients not related to U1 elements
    data = np.concatenate((data, 0.5*np.ones(iga_model.nb_cp-max_ind)),
                          axis=None)
    row = np.concatenate((row, np.arange(max_ind, iga_model.nb_cp)), axis=None)
    col = np.concatenate((col, np.arange(max_ind, iga_model.nb_cp)), axis=None)

    gram_side = sp.coo_matrix((data, (row, col)),
                              shape=(iga_model.nb_cp, iga_model.nb_cp),
                              dtype='float64').tocsc()

    projected_svars = sp.linalg.spsolve(gram_side + gram_side.transpose(),
                                        rhs.transpose())

    pp.generate_proj_vtu(
        *iga_model.get_inputs4proj_vtu('coupling_mono_least_square',
                                       sol.transpose(),
                                       projected_svars.transpose(),
                                       nb_ref=3*np.array([1, 1, 0]),
                                       output_path=tmp_path))


if __name__ == '__main__':
    test_leastsquare_projection_2_patchescplg('results')
