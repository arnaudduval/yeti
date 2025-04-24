# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test optimization of a solid structure structure
Optimization can be run with scipy or nlopt
"""

import numpy as np
import nlopt

from yeti_iga import IgaModel, Patch, ElasticMaterial, \
    BoundaryCondition, DistributedLoad, IgaOptimization, \
    Refinement


class OptimSaver:
    def __init__(self, optim, output_path='.', file_suffix=''):
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
        if file_suffix != '':
            self.file_suffix = '_' + file_suffix
        else:
            self.file_suffix = ''

        self.output_path = output_path

    def save_x_k(self, x):
        """
        Save analysis results for a model modified by a set of design variables x
        This function is called at each iteration of optimization

        Parameters
        ----------
        x : numpy.array(dtype=float)
            Value of design variables, not used here
        """

        print(f'Iteration {self.iopt}')
        self.optim.write_analysis_solution_vtu(f"{self.output_path}/result{self.file_suffix}_{self.iopt:02}.vtu")

        self.optim.write_design_model_control_mesh(f"{self.output_path}/control_net{self.file_suffix}_{self.iopt:02}.vtk")
        self.iopt += 1


def test_optim_solid(tmp_path):
    # Build model
    model = IgaModel('3D solid')

    material = ElasticMaterial(young_modulus=210.e9, poisson_ratio=0.3)

    patch = Patch(element_type='U1',
                  degrees=np.array([1, 1, 1]),
                  knot_vectors=[np.array([0., 0., 1., 1.]),
                                np.array([0., 0., 1., 1.]),
                                np.array([0., 0., 1., 1.])],
                  control_points=np.array([[0., 0., 0.],
                                           [0., 3., 0.],
                                           [0., 0., 1.],
                                           [0., 3., 1.],
                                           [20., 0., 0.],
                                           [20., 3., 0.],
                                           [20., 0., 1.],
                                           [20., 3., 1.]]),
                  weights=np.array([1., 1., 1., 1., 1., 1., 1., 1.]),
                  connectivity=np.array([[7, 6, 5, 4, 3, 2, 1, 0]]),
                  spans=np.array([[1, 1, 1]]),
                  material=material
                  )

    model.add_patch(patch)

    bc = BoundaryCondition(cp_index=np.array([0, 1, 2, 3]),
                            dof=np.array([0, 1, 2]),
                            value=0.
                            )
    model.add_boundary_condition(0, bc)

    dload = DistributedLoad(el_index=np.array([0]),
                            dl_type=40,
                            magnitude=1000.)
    model.add_distributed_load(0, dload)

    model.refine_patch(ipatch=0,
                       nb_subdivision=np.array([0, 0, 2]),
                       nb_degree_elevation=np.array([0, 0, 1]))

    # Get indices of movables CPs
    low_z_cps = model.get_face_cp_indices(3)
    low_y_cps = model.get_face_cp_indices(1)
    movable_cps = np.intersect1d(low_z_cps, low_y_cps)
    symmetry_cps = np.setdiff1d(low_z_cps, low_y_cps)

    nb_var = movable_cps.size*2
    print(nb_var)

    def beam_shape(coords_0, iga_model, x):
        """
        Shape modification function
        Set y and z coordinates of CP on the low z face

               Parameters
        ----------
        coords0 : np.array
            Initial CP coordinates
        iga_model : IgaModel
            IGA model on which shape update applies
        x : np.array(dtype=float)
            Design variables
        """
        iga_model.cp_coordinates[:, :] = coords_0[:, :]
        # Update y coordinate
        iga_model.cp_coordinates[movable_cps, 1] = coords_0[movable_cps, 1] -1.5 + 3.0*x[:nb_var//2]
        iga_model.cp_coordinates[symmetry_cps, 1] = coords_0[symmetry_cps, 1] +1.5 - 3.0*x[:nb_var//2]
        # Update z coordinate
        iga_model.cp_coordinates[movable_cps, 2] = coords_0[movable_cps, 2] -1.0 + 2.0*x[nb_var//2:]
        iga_model.cp_coordinates[symmetry_cps, 2] = coords_0[symmetry_cps, 2] -1.0 + 2.0*x[nb_var//2:]

    def d_beam_shape_dx(coords_0, iga_model, x, ivar):
        """
        Derivative of the shape modification function with respect to the
        design variables

        Parameters
        ----------
        coords0 : np.array
            Initial CP coordinates
        iga_model : IgaModel
            IGA model on which shape update applies
        x : np.array(dtype=float)
            Design variables
        ivar : int
            Index of design variable to derive

        Returns
        -------
        derivative : np.array(dtype=float)
            Derivative of the shape modification fucntipon with respect to
            design variable ivar
        """

        deriv = np.zeros_like(coords_0)

        if ivar < nb_var // 2:
            # Update y coordinate
            deriv[movable_cps[ivar], 1] = 3.0
            deriv[symmetry_cps[ivar], 1] = -3.0
        else:
            # Update z coordinate
            deriv[movable_cps[ivar - nb_var//2], 2] = 2.0
            deriv[symmetry_cps[ivar - nb_var//2], 2] = 2.0

        return deriv

    # Define refinement from design model to analysis model
    refinement = Refinement(model.nb_patch)
    refinement.set_refinement(0, np.array([1, 1, 0]), np.array([2, 1, 0]))

    optim = IgaOptimization(model, nb_var, beam_shape, refinement, d_beam_shape_dx)
    # Handle results saving
    saver = OptimSaver(optim, tmp_path)

    # Compute reference compliance and volume
    x_0 = 0.5*np.ones(nb_var)
    compliance_0 = optim.compliance(x_0)
    volume_0 = optim.volume(x_0)
    point = np.array([0.5, 1.0, 1.0])
    displacement_0 = optim.displacement(x_0, point)
    print(f'{compliance_0 = }')
    print(f'{volume_0 = }')
    print(f'{displacement_0 = }')


    def rel_displacement(x, grad):
        d = optim.displacement(x, point)[2]/displacement_0[2]
        if grad.size > 0:
            # saver.save_x_k(x)
            grad[:] = optim.grad_displacement_analytic(x, point)[:, 2]/displacement_0[2]
        return d - 1.

    def rel_volume(x, grad):
        v = optim.volume(x)/volume_0
        if grad.size > 0:
            saver.save_x_k(x)
            grad[:] = optim.grad_volume_analytic(x)/volume_0
        return v

    minimizer = nlopt.opt(nlopt.LD_SLSQP, nb_var)

    minimizer.set_min_objective(rel_volume)
    minimizer.add_inequality_constraint(rel_displacement, 1e-5 )

    minimizer.set_ftol_rel(1.0e-06)
    minimizer.set_xtol_rel(1.0e-06)
    minimizer.set_maxeval(200)

    minimizer.set_lower_bounds( 0.*np.ones(nb_var))
    minimizer.set_upper_bounds( 0.9*np.ones(nb_var) )

    x = minimizer.optimize(x_0)

    opt_val = minimizer.last_optimum_value()
    result = minimizer.last_optimize_result()
    print(f'{x = }')
    print(f'optimization result : {result}')
    print(f'optimal value : {opt_val}')




if __name__ == '__main__':
    test_optim_solid('.')