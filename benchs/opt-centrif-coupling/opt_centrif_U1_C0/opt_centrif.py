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

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Verify optimisation of a solid structure subjectd to centrifugal body force
Results (design variables values) are compared with reference numerical results.
"""

import numpy as np

# yeti modules
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, OPTmodelling
import yeti_iga.postprocessing.postproc as pp
import yeti_iga.reconstructionSOL as rsol

if __name__ == "__main__":

    modeleIGA = IGAparametrization(filename='centrif_U1_C0')

    # Refinement to create optimisation model
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=np.intp)

    nb_ref[:, 0] = np.array([0, 2, 0])
    nb_deg[:, 0] = np.array([0, 1, 0])
    additional_knots = {"patches": np.array([0]),
                        "1": np.array([]), "2": np.array([0.2,0.2]), "3": np.array([])}

    #additional_knots = {"patches": np.array([0]),
    #                    "1": np.array([]), "2": np.array([]), "3": np.array([])}


    # Initial refinement (none)
    modeleIGA.refine(nb_ref, nb_deg, additional_knots)

    icp = np.where(modeleIGA.coords[1, :] > 1.2)[0]
    nb_var = np.unique(modeleIGA.coords[1, icp]).size

    # Define shape change from design variables
    # min and max dimensions, assuming that design variables are in [0,1]
    MINDIM = 0.1
    MAXDIM = 3.5

    def shapemodif(coords0, igapara, var):
        """"
        Shape parametrization
        Design variables drive coordinates of control point located at y > 1.2
        """
        igapara.coords[:, :] = coords0[:, :]
        # shape change is made on points with y coord higher than 3
        i = 0
        for y_loc in np.unique(modeleIGA.coords[1, icp]):
            # WARNING exact real value comparison is unsafe
            jcp = np.where(modeleIGA.coords[1, :] == y_loc)[0]

            igapara.coords[2, jcp[0]] = - (MINDIM + var[i]*(MAXDIM-MINDIM))/2.
            igapara.coords[2, jcp[1]] = - (MINDIM + var[i]*(MAXDIM-MINDIM))/2.
            igapara.coords[2, jcp[2]] = (MINDIM + var[i]*(MAXDIM-MINDIM))/2.
            igapara.coords[2, jcp[3]] = (MINDIM + var[i]*(MAXDIM-MINDIM))/2.

            i += 1

    # Refinement from optim model to analysis model
    nb_deg = np.array([1, 0, 1])
    nb_ref = np.array([2, 2, 2])

    # Define optim problem
    optPB = OPTmodelling(modeleIGA, nb_var, shapemodif,
                         nb_degreeElevationByDirection=nb_deg,
                         nb_refinementByDirection=nb_ref)

    # OPTIMISATION
    from scipy.optimize import minimize

    # Initial values of design variables
    x0 = ((1.5-MINDIM)/(MAXDIM-MINDIM))*np.ones(nb_var)
    # Compute initial values
    V0 = optPB.compute_volume(x0)
    C0 = optPB.compute_compliance_discrete(x0)

    # Define functions and gradients
    def vol_iga(x_k):
        """
        Compute relative volume variation of the structure
        """
        return (optPB.compute_volume(x_k)-V0)/V0

    def grad_vol_iga(x_k):
        """
        Compute gradient of the relative volume variation with respect to the design variables
        """
        return optPB.compute_gradVolume_AN(x_k)/V0

    def comp_iga(x_k):
        """
        Compute relative compliance of the structure
        """
        return optPB.compute_compliance_discrete(x_k)/C0

    def grad_comp_iga(x_k):
        """
        Compute gradient of the relative compliance with respect to the design variables
        """
        return optPB.compute_gradCompliance_AN(x_k)/C0

    iopt = 0

    def save_x_k(x_k):
        """
        Callback function saving results at each optimization iteration
        """
        global iopt
        print((f'\nIteration {iopt:03}'))
        sol, _ = rsol.reconstruction(
                **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
                        f'OPT3-fine{iopt:02}',  sol.transpose(),
                        nb_ref=3*np.array([1, 1, 1]),
                        Flag=np.array([True, False, False])))
        iopt += 1

    # Bounds for design variables
    bds = ((0., 1.),)*nb_var

    constraint = ({'type': 'eq', 'fun': vol_iga, 'jac': grad_vol_iga})
    x0 = ((1.5-MINDIM)/(MAXDIM-MINDIM))*np.ones(nb_var)

    res = minimize(comp_iga, x0, method='SLSQP',
                   jac=grad_comp_iga, bounds=bds,
                   constraints=constraint, callback=save_x_k)

    # Verify results
    # Numerical reference result
    x_ref = np.array([1.0, 0.843, 0.0, 0.0, 0.0])

    error = np.linalg.norm(res['x']-x_ref)

    print(error)

    assert error < 1.e2
