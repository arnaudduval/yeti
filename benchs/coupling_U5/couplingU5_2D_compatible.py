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


# In this test case, 2 compatible domains are coupled with U5 interface element
# (mortar coupling with projection)
# Each domain is degree 1 and made of a single element, so 4 conditions can be
# expressed (2 coincident control points x 2 dof)
# U5 element is defined with degree 0 and elevated to degree 1 in order to have
# 2 functions for the expression of 2 independant conditions
# A test is made to verify that the correct couling condition are expressed.


# Python module
import sys
import numpy as np
import scipy.sparse as sp

# IGA module
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5

FILENAME = 'testU5_2D_compatible'

modeleIGA = IGAparametrization(filename=FILENAME)

nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
additional_knots = {"patches": np.array([]),
                    "1": np.array([]),
                    "2": np.array([]),
                    "3": np.array([])}

# Interface patch is elevated to degree 1
nb_deg[:, 2] = np.array([1, 0, 0])

modeleIGA.refine(nb_ref, nb_deg, additional_knots)


# STIFFNESS MATRIX
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof]-1

data, row, col, Fb = \
    build_stiffmatrix(*modeleIGA.get_inputs4system_elemStorage())
Kside = sp.coo_matrix((data, (row, col)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside+Kside.transpose()
del Kside, data, row, col

print("stiffness matrix assembly done")

# COUPLING MATRIX
Cdata, Crow, Ccol = cplg_matrixU5(*modeleIGA.get_inputs4cplgmatrixU5())

Cside = sp.coo_matrix((Cdata, (Crow, Ccol)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ctot = Cside + Cside.transpose()

print("Coupling matrix assembly done")

Csolve = Ktot + Ctot

# Coupling conditions for 1st DOF of CPs 2->5 and 4->7 are expressed in lines
# with indices 16 and 18 of Csolve (DOFS 3->9 and 7->13)
sum_lines = Csolve[16, :] + Csolve[18, :]
TOL = 1.E-6
if abs(sum_lines[0, 2] + sum_lines[0, 8]) > TOL:
    sys.exit(-1)
if abs(sum_lines[0, 6] + sum_lines[0, 12]) > TOL:
    sys.exit(-1)

# Coupling conditions for 2nd DOF of CPs 2->5 and 4->7 are expressed in lines
# with indices 17 and 18 of Csolve (DOFS 4->10 and 8->14)
sum_lines = Csolve[17, :] + Csolve[19, :]
if abs(sum_lines[0, 3] + sum_lines[0, 9]) > TOL:
    sys.exit(-1)
if abs(sum_lines[0, 7] + sum_lines[0, 13]) > TOL:
    sys.exit(-1)
