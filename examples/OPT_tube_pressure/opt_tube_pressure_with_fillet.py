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
pressure.
A fillet is present at clamping location
"""

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib import gridspec

# yeti modules
from yeti_iga.preprocessing.igaparametrization import IGAparametrization, OPTmodelling
from yeti_iga.preprocessing.igaparametrization import IGAmanip as manip
import yeti_iga.postprocessing.postproc as pp
import yeti_iga.reconstructionSOL as rsol

if __name__ == "__main__":

    CREATE_PLOTS = True         # Save compliance and volume evolution for plot
    ALLOWED_VOLUME_VAR = 0.025
    ALLOWED_COMP_VAR = 0.05
    INIT_THICKNESS = 10.        # Thickness in the initial model read from file
    MIN_THICKNESS = 2.
    MAX_THICKNESS = 25.

    modeleIGA = IGAparametrization(filename='inputFiles/tube_fillet')

    ### Refinement to create optimization model
    nb_deg_opt = np.zeros((3, modeleIGA.nb_patch), dtype=int)
    nb_ref_opt = np.zeros((3, modeleIGA.nb_patch), dtype=int)

    # Tube and hull (they share directions v and w)
    nb_deg_opt[:, 0] = np.array([0, 0, 1])
    nb_deg_opt[:, 1] = np.array([0, 0, 1])
    nb_ref_opt[:, 0] = np.array([0, 0, 3])
    nb_ref_opt[:, 1] = np.array([0, 0, 3])


    ### Refinement of analysis model (defined form initial parametrization)
    nb_deg_an = np.zeros((3, modeleIGA.nb_patch), dtype=int)
    nb_ref_an = np.zeros((3, modeleIGA.nb_patch), dtype=int)

    # Tube and Hull (they share directions v and w)
    nb_ref_an[:, 0] = np.array([2, 3, 4])
    nb_ref_an[:, 1] = np.array([0, 3, 4])
    nb_deg_an[:, 0] = np.array([1, 0, 1])
    nb_deg_an[:, 1] = np.array([0, 0, 1])
    # Fillet
    nb_ref_an[:, 2] = np.array([2, 3, 2])
    nb_deg_an[:, 2] = np.array([0, 1, 0])
    # Lagrange
    nb_ref_an[:, 3] = np.array([3, 2, 0])
    nb_deg_an[:, 3] = np.array([1, 1, 0])

    # Create optimization model
    modeleIGA.refine(nb_ref_opt, nb_deg_opt)

    # Get number of control points along w direction on the tube
    # (assume 2 CP in tube width)
    nb_var = int(len(manip.get_directionCP(modeleIGA, 3, 0))/2) - 1

    # Get indices of control points on external diameter of the tube
    ext_cps = manip.get_directionCP(modeleIGA, 2, 0)[3:]-1
    n_pt_per_height = len(ext_cps) // nb_var

    assert n_pt_per_height == 3

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

    # Define optimization problem
    opt_pb = OPTmodelling(modeleIGA, nb_var, shapemodif,
                          nb_degreeElevationByDirection=nb_deg_an-nb_deg_opt,
                          nb_refinementByDirection=nb_ref_an-nb_ref_opt,
                          fct_dervShapeParam=shapemodif_der)

    # Initial value of design variables
    x0 = ((INIT_THICKNESS-MIN_THICKNESS)/(MAX_THICKNESS-MIN_THICKNESS))*np.ones(nb_var)

    # Flags for patch taken into account for volume computation
    listpatch_vol = [1, 0, 1, 0]

    # Initial volume and compliance
    v0 = opt_pb.compute_volume(x0, listpatch=listpatch_vol)
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
        return ALLOWED_VOLUME_VAR \
            - abs(opt_pb.compute_volume(x_k, listpatch=listpatch_vol) - v0) / v0

    def grad_var_vol_iga(x_k):
        """
        Compute gradient of var_vol_iga
        """
        v_k = opt_pb.compute_volume(x_k, listpatch=listpatch_vol)
        return -np.sign(v_k - v0)*opt_pb.compute_gradVolume_AN(x_k, listpatch=listpatch_vol)/v0

    def vol_iga(x_k):
        """
        Compute relative volume of the structure
        """
        return opt_pb.compute_volume(x_k, listpatch=listpatch_vol)/v0

    def grad_vol_iga(x_k):
        """
        Compute gradient of vol_iga
        """
        return opt_pb.compute_gradVolume_AN(x_k, listpatch=listpatch_vol)/v0


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
        Compute relative compliance of the structure
        """
        return opt_pb.compute_compliance_discrete(x_k)/c0

    def grad_comp_iga(x_k):
        """
        Compute gradient of comp_iga
        """
        return opt_pb.compute_gradCompliance_AN(x_k)/c0

    iopt = 0
    comp_values = []
    vol_values = []

    def save_xk(x_k):
        """
        Callback function saving results at each optimization iteration
        """
        global iopt
        print((f'\nIteration {iopt:03}'))
        pp.generatevtu(*opt_pb.coarse_parametrization.get_inputs4postprocVTU(
            f'opt_coarse_w_fillet{iopt:02}', np.zeros_like(opt_pb.coarse_parametrization.coords),
            nb_ref=6*np.array([1, 1, 1]),Flag=np.array([False]*3)))

        opt_pb.coarse_parametrization.generate_vtk4controlMeshVisu(
            f'opt_ctrl_coarse_w_fillet{iopt:02}',0)

        sol, _ = rsol.reconstruction(
                **opt_pb.fine_parametrization.get_inputs4solution(opt_pb.save_sol_fine))
        pp.generatevtu(*opt_pb.fine_parametrization.get_inputs4postprocVTU(
                        f'opt_fine_w_fillet{iopt:02}',  sol.transpose(),
                        nb_ref=3*np.array([1, 1, 1]),
                        Flag=np.array([True, False, False])))

        # Postprocessing of hull
        if 'U0' in opt_pb.coarse_parametrization._ELT_TYPE:
            # Set refinement and variables to process
            output = np.array([False, False, False])
            nb_ref_visu = np.array([6, 6, 6])
            # Loop on U0 patches
            ii = 1
            for idx in np.where(opt_pb.coarse_parametrization._ELT_TYPE == 'U0')[0]:
                filename = \
                    'couplingU5_{}_hullobj{}'.format(ii, iopt)
                # # Generate control mesh (.vtk. file)
                # modeleIGA.generate_vtk4controlMeshVisu(filename, idx)
                # Initialise
                idx = np.array([idx])
                inputs = [filename, output, nb_ref_visu]
                # Get geometrical settings
                geoSet = opt_pb.coarse_parametrization.get_geometricSettings_somePatch(idx)
                # Get mechanical settings
                mechSet = opt_pb.coarse_parametrization.get_mechanicalSettings_somePatch(idx)
                inputs.append(np.zeros_like(np.array(mechSet[1][0])))
                # Build inputs list to pass to function
                inputs.append(mechSet[1][0])  # COORDS
                inputs.append(np.concatenate(mechSet[2]).ravel())  # IEN
                inputs.append(geoSet[6])  # elementsByPatch
                inputs.append(geoSet[1])  # Nkv
                inputs.append(np.concatenate(geoSet[2][0]).ravel())  # Ukv
                inputs.append(geoSet[4])  # Nijk
                inputs.append(np.concatenate(geoSet[5]).ravel())  # weigth
                inputs.append(geoSet[3])  # Jpqr
                inputs.append('U1')  # ELT_TYPE
                inputs.append(np.array([[0.0], [0.0]]))  # MATERIAL PROPERTIES
                inputs.append('/THREED')  # TENSOR
                inputs.append(mechSet[4][0][0].tolist())  # PROPS
                inputs.append(mechSet[4][1])  # JPROPS
                inputs.append(geoSet[8])  # nnode
                inputs.append(geoSet[7])  # nb_patch
                inputs.append(geoSet[9])  # nb_elem
                inputs.append(mechSet[1][1])  # nb_cp
                inputs.append(mechSet[0][3])  # MCRD
                # Generate .vtu file
                pp.generatevtu(*inputs)
                ii += 1

        if CREATE_PLOTS:
            comp_values.append(opt_pb.compute_compliance_discrete(x_k))
            vol_values.append(opt_pb.compute_volume(x_k, listpatch=listpatch_vol))

        iopt += 1

    # Bounds for design variables
    bounds = ((0., 1.),) * nb_var

    save_xk(x0)

    # OPTIM COMPLIANCE A VARIATION DE VOLUME DONNEE
    constraints = ({'type': 'ineq', 'fun': var_vol_iga, 'jac': grad_var_vol_iga})
    res = minimize(comp_iga, x0, method='SLSQP',
                   jac=grad_comp_iga, bounds=bounds,
                   constraints=constraints, callback=save_xk)

    # OPTIM VOLUME A VARIATION DE COMPLIANCE DONNEE
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

    v_f = opt_pb.compute_volume(res['x'], listpatch=listpatch_vol)
    c_f = opt_pb.compute_compliance_discrete(res['x'])
    print('Volume: ', v0, '->', v_f, (100.*(v_f - v0)/v0), ' %')
    print('Compliance: ', c0, '->',  c_f, (100.*(c_f - c0)/c0), ' %')

    if CREATE_PLOTS:
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])

        ax0 = plt.subplot(gs[0])
        ax0.plot(range(len(comp_values)), comp_values)
        # ax0.set_xlim([0, nopt])
        ax0.set_ylim(0, np.amax(comp_values))
        ax0.set_ylabel('Compliance')

        ax1 = plt.subplot(gs[1], sharex = ax0)
        ax1.plot(range(len(comp_values)), vol_values)
        ax1.set_ylim(0, np.amax(vol_values))
        ax1.set_ylabel('Volume')
        ax1.set_xlabel('Iterations')

        plt.savefig('results/graphs_evol.png')
