# Copyright 2022 Arnaud Duval

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

# -*- coding: utf-8 -*-

"""
A test to validate compliance gradient computation for embedded solid element.
"""

import numpy as np
import scipy.sparse as sp
from termcolor import colored

from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import OPTmodelling
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp




if __name__ == "__main__":

    # modeleIGA = IGAparametrization(filename='embd_cube_lin_press')
    # modeleIGA = IGAparametrization(filename='unit_square')
    # modeleIGA = IGAparametrization(filename='embd_cube_nonli_hull_press')
    # modeleIGA = IGAparametrization(filename='embd_cube_nonli_hull_unit_embd')
    # modeleIGA = IGAparametrization(filename='embd_cube_nonli_embd_press')
    modeleIGA = IGAparametrization(filename='test')

    # Set arguments for model refinement
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    additional_knots = {"patches": np.array([]),
                        "1": np.array([]),
                        "2": np.array([]),
                        "3": np.array([])}

    # nb_deg[:, 0] = [1, 1, 1]
    # nb_ref[:, 0] = [1, 1, 1]
    
    # nb_deg[:, 1] = [1, 1, 1]
    # nb_ref[:, 1] = [1, 1, 1]
    modeleIGA.refine(nb_ref, nb_deg, additional_knots)

    # # Static study
    # ndof = modeleIGA._nb_dof_free
    # idof = modeleIGA._ind_dof_free[:ndof] - 1
    # print(idof)
    # print(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_free)
    # data, row, col, Fb = build_stiffmatrix(
    #                         *modeleIGA.get_inputs4system_elemStorage())
    # Kside = sp.coo_matrix((data, (row, col)),
    #                     shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
    #                     dtype='float64').tocsc()
    # Ktot = Kside + Kside.transpose()
    # K2solve = Ktot[idof, :][:, idof]
    # del Kside, data, row, col

    # x = sp.linalg.spsolve(K2solve, Fb[idof])

    # # Solution reconstruction
    # SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

    # pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    #     'test_embd_cube',
    #     SOL.transpose(), nb_ref=[3, 3, 3], Flag=[True, False, False]))
    # modeleIGA.generate_vtk4controlMeshVisu(
    #     'test_embd_cube_map_cp', 0)
    # modeleIGA.generate_vtk4controlMeshVisu(
    #     'test_embd_cube_embded_cp', 1)

    # # Fin static study

    # exit()



    initcoords = modeleIGA.coords

    # Build optim model
    def shapemodif(coords0, igapara, var):
        """
        A shape change that apply on all coordinates of all control points
        """
        igapara.coords[:, :] = coords0[:, :]
        dim = coords0.shape[0]
        n_cp = coords0.shape[1]
        i = 0

        for i_cp in range(n_cp):
            for i_dim in range(dim):
                igapara.coords[i_dim, i_cp] = initcoords[i_dim, i_cp] + var[i]
                i += 1

    # Refinement from optim model to analysis model
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)

    nb_var = modeleIGA.coords.size

    optPB = OPTmodelling(modeleIGA, nb_var, shapemodif,
                         nb_degreeElevationByDirection=nb_deg,
                         nb_refinementByDirection=nb_ref)

    x0 = np.zeros((nb_var))

    c0 = optPB.compute_compliance_discrete(x0)
    gradC_DF = optPB.compute_gradCompliance_FD(x0)
    gradC_AN = optPB.compute_gradCompliance_AN(x0)

    #exit()

    error = np.linalg.norm(gradC_DF - gradC_AN) / c0
    print(f"Compliance : {c0:.02E}")
    print(f"Error : {error:.02E}")

    for i_var in range(nb_var):
        success = (abs(gradC_DF[i_var]/c0) < 1.e-5) or          \
                    (abs((gradC_AN[i_var]-gradC_DF[i_var])/gradC_DF[i_var]) < 1.e-3)
        if success:
            print(i_var, int(i_var/3)+1, gradC_AN[i_var], gradC_DF[i_var])
        else:
            print(colored((i_var, int(i_var/3)+1, gradC_AN[i_var], gradC_DF[i_var]),'red'))

    error = np.linalg.norm(gradC_DF - gradC_AN)/c0

    print(error)              

    # print(modeleIGA._IEN)
    # print(modeleIGA._ind_dof_bloq)
    # assert error < 1.e-6
