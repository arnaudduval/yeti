# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test optimization of a shell structure
Optimization can be run with scipy or nlopt
With scipy, the result is not what is expected
"""

import numpy as np
from scipy.optimize import minimize
import nlopt

from yeti_iga import IgaModel, Patch, ElasticMaterial, \
    BoundaryCondition, DistributedLoad, IgaOptimization, \
    Refinement

class OptimSaver:
    def __init__(self, optim, file_suffix=''):
        """
        Parameters
        ----------
        optim : IgaOptimization
            IgaOptimization objetc on which this saver apply
        file_suffix : string
            Suffix added to result file name
        """
        self.iopt = 0
        self.optim = optim
        self.file_suffix = '_' + file_suffix

    def save_x_k(self, x):
        """
        Save analysis results for a model modified by a set of design variables x
        This function is called at each iteration of optimization

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Value of design variables
        """
        tmp_path = '.'
        print(f'Iteration {self.iopt}')
        self.optim.write_analysis_solution_vtu(x,
                                          f"{tmp_path}/result{self.file_suffix}_{self.iopt:02}.vtu")
        self.iopt += 1


def test_optim_shell():
    # Build model
    model = IgaModel("3D shell")
    material = ElasticMaterial(young_modulus=210.e9, poisson_ratio=0.3)
    patch = Patch(element_type='U3',
                  degrees=np.array([1, 1]),
                  knot_vectors=[np.array([0., 0., 1., 1.]),
                                np.array([0., 0., 1., 1.])],
                  control_points=np.array([[0., 0., 0.],
                                           [10., 0., 0.],
                                           [0., 10., 0.],
                                           [10., 10., 0.]]),
                  weights=np.array([1., 1., 1., 1.]),
                  connectivity=np.array([[3, 2, 1, 0]]),
                  spans=np.array([1, 1]),
                  material=material,
                  properties=np.array([0.1])
                  )
    model.add_patch(patch)
    for icp in range(4):
        bc = BoundaryCondition(cp_index=np.array([icp]),
                               dof=np.array([0, 1, 2]),
                               value=0.
                               )
        model.add_boundary_condition(0, bc)

    dload = DistributedLoad(el_index=np.array([0]),
                            dl_type=66,
                            magnitude=-500.)
    model.add_distributed_load(0, dload)

    model.refine_patch(ipatch=0,
                       nb_degree_elevation=np.array([1, 1]),
                       nb_subdivision=np.array([2, 2])
                       )

    # Get indices of control points located at corners
    corner_cps = model.get_corners_cp_indices()
    # Deduce indices of control points NOT located at corers
    free_cps = np.setxor1d(model.cp_indices, corner_cps)
    # # Elevate those control points along z axis
    # model.cp_coordinates[free_cps, 2] += 0.5

    # Design variables : free CPs can move along z axis
    nb_var = free_cps.size


    def altitude(coords_0, iga_model, x):
        # TODO use directly an IgaModel object
        """
        Shape modification function
        Set z coordinate of CP (except corners)

        Parameters
        ----------
        coords0 : np.array
            Initial CP coordinates
        iga_model : IgaModel
            IGA model on which shape update applies
        x : np.array(dtype=float)
            Design variables
        """

        model.cp_coordinates[:, :] = coords_0[:, :]
        model.cp_coordinates[free_cps, 2] = coords_0[free_cps, 2] + 5.*x[:]




    # Define refinement from design model to analysis model
    refinement = Refinement(model.nb_patch)
    refinement.set_refinement(0, np.array([0, 0]), np.array([3, 3]))

    optim = IgaOptimization(model, nb_var, altitude, refinement)
    # Handle results saving
    solver = 'scipy'
    solver = 'nlopt'
    saver = OptimSaver(optim, solver)

    # Compute reference compliance and volume
    compliance_0 = optim.compliance(np.zeros(nb_var))
    volume_0 = 1.1 * optim.volume(np.zeros(nb_var))



    if solver == 'nlopt':
        def rel_compliance(x, grad):
            c = optim.compliance(x)/compliance_0
            if grad.size > 0:
                saver.save_x_k(x)
                grad[:] = optim.grad_compliance_analytic(x)/compliance_0


            return c

        def rel_volume(x, grad):
            if grad.size > 0:
                grad[:] = optim.grad_volume_analytic(x)/volume_0
            return optim.volume(x)/volume_0 - 1.


        x_0 = np.ones(nb_var)*0.1

        minimizer = nlopt.opt(nlopt.LD_SLSQP, nb_var)

        minimizer.set_min_objective(rel_compliance)
        minimizer.add_inequality_constraint(rel_volume, 1e-5 )

        minimizer.set_ftol_rel(1.0e-06)
        minimizer.set_xtol_rel(1.0e-06)
        minimizer.set_maxeval(200)

        minimizer.set_lower_bounds( 0.*np.ones(nb_var))
        minimizer.set_upper_bounds( 1.*np.ones(nb_var) )

        x = minimizer.optimize(x_0)

    if solver == 'scipy':
        def rel_compliance(x):
            """
            Compute relative compliance of the model for a given set of design
            variables x
            """
            return optim.compliance(x)/compliance_0

        def rel_volume(x):
            """
            Compute relative volume of the model for a given set of design
            variables x
            """
            return optim.volume(x)/volume_0 - 1.

        def grad_rel_volume(x):
            """
            Compute gradient of the relative volume with respect to a set of
            design variables x
            """
            return optim.grad_volume_analytic(x)/volume_0

        def grad_rel_compliance(x):
            """
            Compute gradient of the remative compliance with respect to a set of
            design variables x
            """
            return optim.grad_compliance_analytic(x)/compliance_0

        # Initial value of design variables for optimization
        x_0 = np.ones(nb_var) * 0.1

        # Set bounds for design variable [0., 1.]
        bounds = ((0., 1.),) * nb_var

        saver.save_x_k(x_0)

        # Set up optimization : minimize compliance and keep volmume < 1.1 x initial volume
        constraints = ({'type': 'ineq', 'fun': rel_volume, 'jac': grad_rel_volume})
        result = minimize(rel_compliance,
                        x_0,
                        method='SLSQP',
                        jac=grad_rel_compliance,
                        bounds=bounds,
                        constraints=constraints,
                        callback=saver.save_x_k)








if __name__ == '__main__':
    test_optim_shell()