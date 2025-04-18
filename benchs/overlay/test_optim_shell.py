# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test optimization of a shell structure
"""

import numpy as np

from yeti_iga import IgaModel, Patch, ElasticMaterial, \
    BoundaryCondition, DistributedLoad, IgaOptimization, \
    Refinement


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
        bc = BoundaryCondition(cp_index=np.array([0]),
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
    # Elevate those control points along z axis
    model.cp_coordinates[free_cps, 2] += 0.5

    # Design variables : free CPs can move along z axis
    nb_var = free_cps.size

    def altitude(coords_0, iga_param, x):
        # TODO use directly an IgaModel object
        """
        Shape modification function
        Set z coordinate of CP (except corners)

        Parameters
        ----------
        coords0 : np.array
            Initial CP coordinates
        iga_param : IgaParametrization
            Legacy IgaParametrization object
            TODO Use IgaModel instead
        x : np.array(dtype=float)
            Design variables
        """

        iga_param.coords[:, :] = coords_0[:, :]
        # TODO if using IgaModel object, coords are translated
        iga_param.coords[2, free_cps] = coords_0[2, free_cps] + 5.*x[:]




    # Define refinement from design model to analysis model
    refinement = Refinement(model.nb_patch)
    refinement.set_refinement(0, np.array([1, 1]), np.array([2, 2]))

    optim = IgaOptimization(model, nb_var, altitude, refinement)



if __name__ == '__main__':
    test_optim_shell()