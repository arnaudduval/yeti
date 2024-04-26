# Copyright 2023 Arnaud Duval

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
Optimize thickness of a quarter cylinder subjected to internal hydrostatic
pressure
"""

import numpy as np
from scipy.optimize import minimize

# yeti modules
from preprocessing.igaparametrization import IGAparametrization, OPTmodelling
from preprocessing.igaparametrization import IGAmanip as manip
import postprocessing.postproc as pp
import reconstructionSOL as rsol

if __name__ == "__main__":

    ALLOWED_VOLUME_VAR = 0.025
    ALLOWED_COMP_VAR = 0.05
    INIT_THICKNESS = 10.        # Thickness in the initial model
    MIN_THICKNESS = 2.
    MAX_THICKNESS = 25.

    modeleIGA = IGAparametrization(filename='inputFiles/tube')

    # Refinement to create optimization model
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=int)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=int)
    nb_deg[:, 0] = np.array([0, 0, 1])
    nb_ref[:, 0] = np.array([0, 0, 4])

    modeleIGA.refine(nb_ref, nb_deg)

    # Get number of control points along w direction
    # (assume 2 CP in tube width)
    nb_var = int(len(manip.get_directionCP(modeleIGA, 3))/2)

    # Get indices of control points on external diameter
    ext_cps = manip.get_directionCP(modeleIGA, 2)-1
    n_pt_per_height = len(ext_cps) // nb_var

    assert n_pt_per_height == 3

    # Get indices of control points on intrernal diameter
    int_cps = manip.get_directionCP(modeleIGA, 1)-1

    def shapemodif(coords0, igapara, var):
        """
        Shape parametrization
        Design variables drive radial coordinates of control points for the
        surface defining the external diameter
        """

        igapara.coords[:, :] = coords0[:, :]
        for i, ctrlpt in enumerate(ext_cps):
            new_thickness = MIN_THICKNESS + var[i // n_pt_per_height]*(MAX_THICKNESS-MIN_THICKNESS)
            if i % 3 == 0:
                igapara.coords[0, ctrlpt] = coords0[0, ctrlpt] + new_thickness - INIT_THICKNESS
            elif i % 3 == 1:
                igapara.coords[:2, ctrlpt] = coords0[:2, ctrlpt] + new_thickness - INIT_THICKNESS
            elif i % 3 == 2:
                igapara.coords[1, ctrlpt] = coords0[1, ctrlpt] + new_thickness - INIT_THICKNESS

    def shapemodif_der(coords0, _, __, i_var):
        """
        Derivative of shape parametrization with respect to design varables
        """
        d_shape = np.zeros_like(coords0)
        for i, ctrlpt in enumerate(ext_cps):
            if i // n_pt_per_height == i_var:
                if i % 3 == 0:
                    d_shape[0, ctrlpt] = MAX_THICKNESS-MIN_THICKNESS
                elif i % 3 == 1:
                    d_shape[:2, ctrlpt] = MAX_THICKNESS-MIN_THICKNESS
                elif i % 3 == 2:
                    d_shape[1, ctrlpt] = MAX_THICKNESS-MIN_THICKNESS

        return d_shape

    # Refinement from optim model to analysis model
    nb_deg = np.zeros((3, modeleIGA.nb_patch), dtype=int)
    nb_ref = np.zeros((3, modeleIGA.nb_patch), dtype=int)
    nb_deg[:, 0] = np.array([1, 0, 0])
    nb_ref[:, 0] = np.array([2, 2, 0])

    # Define optim problem
    opt_pb = OPTmodelling(modeleIGA, nb_var, shapemodif,
                          nb_degreeElevationByDirection=nb_deg,
                          nb_refinementByDirection=nb_ref,
                          fct_dervShapeParam=shapemodif_der)

    # Initial value of design variables
    x0 = ((INIT_THICKNESS-MIN_THICKNESS)/(MAX_THICKNESS-MIN_THICKNESS))*np.ones(nb_var)

    # Initial volume and compliance
    v0 = opt_pb.compute_volume(x0)
    c0 = opt_pb.compute_compliance_discrete(x0)

    print('Initial volume : ', v0)
    print('Initial compliance : ', c0)

    # Define functions and gradients
    def var_vol_iga(x_k):
        """
        Compute relative volume variation of the structure compared to an
        allowed value.
        Positive if volume variation is lower than allowed
        """
        return ALLOWED_VOLUME_VAR - abs(opt_pb.compute_volume(x_k) - v0) / v0

    def grad_var_vol_iga(x_k):
        """
        Compute gradient of var_vol_iga
        """
        v_k = opt_pb.compute_volume(x_k)
        return -np.sign(v_k - v0)*opt_pb.compute_gradVolume_AN(x_k)/v0

    def vol_iga(x_k):
        """
        Compute relative volume of the structure
        """
        return opt_pb.compute_volume(x_k)/v0

    def grad_vol_iga(x_k):
        """
        Compute gradient of vol_iga
        """
        return opt_pb.compute_gradVolume_AN(x_k)/v0


    def var_comp_iga(x_k):
        """
        Compute relative compliance variation of the structure compared to an
        allowed value
        Positive if compliance variation is lower than allowed
        """
        return ALLOWED_COMP_VAR - abs(opt_pb.compute_compliance_discrete(x_k) - c0) / c0

    def grad_var_comp_iga(x_k):
        """
        Compute gradient of var_comp_iga
        """
        c_k = opt_pb.compute_compliance_discrete(x_k)
        return -np.sign(c_k - c0)*opt_pb.compute_gradCompliance_AN(x_k)/c0

    def comp_iga(x_k):
        """
        Compute relative complianjce of the structure
        """
        return opt_pb.compute_compliance_discrete(x_k)/c0

    def grad_comp_iga(x_k):
        """
        Compute gradient of comp_iga
        """
        return opt_pb.compute_gradCompliance_AN(x_k)/c0

    iopt = 0

    def save_xk(_):
        """
        Callback function saving results at each optimization iteration
        """
        global iopt
        print((f'\nIteration {iopt:03}'))
        pp.generatevtu(*opt_pb.coarse_parametrization.get_inputs4postprocVTU(
            f'opt-coarse{iopt:02}', np.zeros_like(opt_pb.coarse_parametrization.coords),
            nb_ref=3*np.array([1, 1, 1]),Flag=np.array([False]*3)))

        opt_pb.coarse_parametrization.generate_vtk4controlMeshVisu(f'opt-coarse{iopt:02}',0)

        sol, _ = rsol.reconstruction(
                **opt_pb.fine_parametrization.get_inputs4solution(opt_pb.save_sol_fine))
        pp.generatevtu(*opt_pb.fine_parametrization.get_inputs4postprocVTU(
                        f'opt-fine{iopt:02}',  sol.transpose(),
                        nb_ref=3*np.array([1, 1, 1]),
                        Flag=np.array([True, False, False])))
        iopt += 1

    # Bounds for design variables
    bounds = ((0., 1.),) * nb_var

    save_xk(x0)

    # OPTIM COMPLIANCE A VARIATION DE VOLUME DONNEE
    constraints = ({'type': 'ineq', 'fun': var_vol_iga, 'jac': grad_var_vol_iga})
    res = minimize(comp_iga, x0, method='SLSQP',
                   jac=grad_comp_iga, bounds=bounds,
                   constraints=constraints, callback=save_xk)

    # OPTIM VOLUME A COMPLIANCE CONSTANTE
    # constraints = ({'type': 'eq', 'fun': comp_iga, 'jac': grad_comp_iga})
    # res = minimize(vol_iga, x0, method='SLSQP',
    #                jac=grad_vol_iga, bounds=bounds,
    #                constraints=constraints, callback=save_xk)

    # OPTIM VOLUME A VARIATIONDE COMPLIANCE DONNEE
    # constraints = ({'type': 'ineq', 'fun': var_comp_iga, 'jac': grad_var_comp_iga})
    # res = minimize(vol_iga, x0, method='SLSQP',
    #                jac=grad_vol_iga, bounds=bounds,
    #                constraints=constraints, callback=save_xk)

    print('Optimization succes : ', res['success'])
    print(res['message'])
    print('Result : ', res['x'])
    print('Objective function value : ', res['fun'])
    print('# of evaluations of objective function : ', res['nfev'])
    print('# of evaluations of jacobian : ', res['njev'])

    v_f = opt_pb.compute_volume(res['x'])
    c_f = opt_pb.compute_compliance_discrete(res['x'])
    print('Volume: ', v0, '->', v_f, (100.*(v_f - v0)/v0), ' %')
    print('Compliance: ', c0, '->',  c_f, (100.*(c_f - c0)/c0), ' %')
