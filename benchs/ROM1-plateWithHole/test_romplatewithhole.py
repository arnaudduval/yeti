# Copyright 2022 Thibaut Hirschler

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
Use a parameterized model.
Create a simple Reduced Order Model.
"""

import os

import numpy as np
import scipy.sparse as sp
from scipy.stats import qmc
from scipy.linalg import svd as scipySVD, solve as scipy_solve
from scipy import interpolate

# pylint: disable=no-name-in-module

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static \
    as build_stiffmatrix


def test_romplatewithhole():
    """
    Build a ROM and compare with reference IGA solution for several samples.
    """

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(
        filename=f'{script_dir}/parametericPlateWithHole')

    # Refine model
    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_deg[0, :] = 0
    nb_deg[1, :] = 1
    nb_ref[:2, :] = 1
    iga_model.refine(nb_ref, nb_deg)

    # Get design parameters ('<DESIGNVAR>' in .inp file)
    designvar_name = list(iga_model.design_parameters.keys())[0]
    designvar_value = iga_model.design_parameters[designvar_name]

    # a. Initial geometry

    # b. Shape modified geometry
    # The shape modification is done in two steps:
    # - step1: assign the new value of the design variable.
    # - step2: call `iga_model.shapeupdate()` to update the control points.
    designvar_value = 0.25
    iga_model.design_parameters[designvar_name] = designvar_value   # step1
    iga_model.shapeupdate()     # step2

    # Let us try to build a simple Reduced Order Model by using this
    # parameterized model.

    # Preliminaries: define a function to automatically build the FE systems.
    ndof = iga_model.nb_dof_free
    idof = iga_model.ind_dof_free[:ndof]-1

    def buildsystem(designvar):
        '''
        Build the finite element operators for a given geometric configuration.

        Parameters
        ----------
        designvar : float
            The value of the design variable, i.e. the radius of the hole.
            Should be between 0.05 and 0.95.

        Returns
        -------
        mat_k : sparse matrix
            Stiffness matrix.
        vec_f : array of floats
            Load vector.
        '''
        if designvar < 0.05 or designvar > 0.95:
            raise ValueError(
                'Input parameter should be a float in interval [0.05, 0.95].'
            )

        iga_model.design_parameters[designvar_name] = designvar
        iga_model.shapeupdate()

        data, row, col, rhs = build_stiffmatrix(
            *iga_model.get_inputs4system_elemStorage())
        k_side = sp.coo_matrix(
            (data, (row, col)),
            shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
            dtype='float64').tocsc()
        k_tot = k_side + k_side.transpose()
        mat_k = k_tot[idof, :][:, idof]
        mat_k.sort_indices()
        vec_f = rhs[idof]
        return mat_k, vec_f

    def refsolution(designvar):
        '''Perform the IGA for a given geometrical configuration.

        Parameters
        ----------
        designvar : float
            The value of the design variable, i.e. the radius of the hole.
            Should be between 0.05 and 0.95.

        Returns
        -------
        vec_u : array of floats
            Solution vector.
        '''
        mat_k, vec_f = buildsystem(designvar)
        vec_u = sp.linalg.spsolve(mat_k, vec_f)
        return vec_u

    # Step1: build snapshots of K and F.

    rmin = 0.05
    rmax = 0.95
    sample = np.ravel((rmax-rmin)*qmc.LatinHypercube(1).random(n=100)+rmin)
    sample = np.append(sample, [0.05, 0.95])

    snapshots_k = []
    snapshots_f = []
    for i in range(sample.size):
        mat_k, vec_f = buildsystem(sample[i])
        snapshots_k.append([mat_k.data.copy()])
        snapshots_f.append([vec_f.copy()])
    snapshots_k = np.block(snapshots_k).T
    snapshots_f = np.block(snapshots_f).T

    # Step2: build a reduced basis for K and F.
    # We do it by using a truncated SVD.
    rtol = 1e-6     # tolerance for truncating the basis.

    u, s, v_h = scipySVD(snapshots_k)
    n_basis_k = np.count_nonzero(s/s.max() > rtol)
    basis_k = u[:, :n_basis_k]*s[:n_basis_k]
    coef_k = [interpolate.interp1d(sample, v_h[i, :], 'cubic')
              for i in range(n_basis_k)]

    u, s, v_h = scipySVD(snapshots_f)
    n_basis_f = np.count_nonzero(s/s.max() > rtol)
    basis_f = u[:, :n_basis_f]*s[:n_basis_f]
    coef_f = [interpolate.interp1d(sample, v_h[i, :], 'cubic')
              for i in range(n_basis_f)]

    # Check approximation
    # Define a function to automatically build the FE systems as given by
    # the ROM.
    indices = mat_k.indices.copy()
    indptr = mat_k.indptr.copy()

    def buildapproxsystem(designvar):
        '''Build the finite element operators as approximated by the ROM for
        a given geometrical configuration.

        Parameters
        ----------
        designvar : float
            The value of the design variable, i.e. the radius of the hole.
            Should be between 0.05 and 0.95.

        Returns
        -------
        mat_k : sparse matrix
            Stiffness matrix.
        vec_f : array of floats
            Load vector.
        '''
        if designvar < 0.05 or designvar > 0.95:
            raise ValueError(
                'Input parameter should be a float in interval [0.05, 0.95].'
            )

        iga_model.design_parameters[designvar_name] = designvar
        iga_model.shapeupdate()

        alpha_k = np.array([c(designvar) for c in coef_k])
        data_k = basis_k.dot(alpha_k)
        mat_k = sp.csc_matrix((data_k, indices, indptr))

        alpha_f = np.array([c(designvar) for c in coef_f])
        vec_f = basis_f.dot(alpha_f)
        return mat_k, vec_f

    designvar_value = 0.45

    # Step3: build a reduced basis for U.

    rmin = 0.05
    rmax = 0.95
    sample = np.ravel((rmax-rmin)*qmc.LatinHypercube(1).random(n=100)+rmin)
    sample = np.append(sample, [0.05, 0.95])

    snapshots_u = []
    for i in range(sample.size):
        mat_k, vec_f = buildapproxsystem(sample[i])
        vec_u = sp.linalg.spsolve(mat_k, vec_f)
        snapshots_u.append([vec_u.copy()])
    snapshots_u = np.block(snapshots_u).T
    print('  Done.')

    # Step4: build a reduced basis for U.
    # We do it by using a truncated SVD.
    rtolsol = 1e-4  # tolerance for truncating the basis.

    u, s, v_h = scipySVD(snapshots_u)
    n_basis_u = np.count_nonzero(s/s.max() >= rtolsol)
    basis_u = u[:, :n_basis_u]

    def romsolution(designvar):
        """
        Get the solution using the Reduced Order Model.

        Parameters
        ----------
        designvar : float
            The value of the design variable, i.e. the radius of the hole.
            Should be between 0.05 and 0.95.

        Returns
        -------
        vec_u : array of floats
            Solution vector.
        """

        if designvar < 0.05 or designvar > 0.95:
            raise ValueError(
                'Input parameter should be a float in interval [0.05, 0.95].'
            )

        mat_k, vec_f = buildapproxsystem(designvar)

        mat_uku = basis_u.T.dot(mat_k.dot(basis_u))
        vec_uf = basis_u.T.dot(vec_f)
        alpha_u = scipy_solve(mat_uku, vec_uf, assume_a='pos')

        vec_u = basis_u.dot(alpha_u)
        return vec_u

    # Check approximation
    sample = np.linspace(rmin, rmax, 15)

    for designvar_value in sample:
        u_iga = refsolution(designvar_value)
        u_rom = romsolution(designvar_value)
        assert np.linalg.norm(u_iga - u_rom)/np.linalg.norm(u_iga) < 1.e-3


if __name__ == '__main__':
    test_romplatewithhole()
