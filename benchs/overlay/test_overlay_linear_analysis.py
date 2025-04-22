# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test definition of an IGA model for linear elasticity analysis
"""

import numpy as np
import scipy.sparse as sp

from yeti_iga import IgaModel, Patch, ElasticMaterial, \
    BoundaryCondition, DistributedLoad


def test_linear_analysis_single_patch_3d_solid(tmp_path):
    """
    Build a 3D model, refine it with degree elevation, knot insertion and
    subdivision.
    Compute linear analysis solution
    Write VTU result file
    """

    model = IgaModel("3D solid")

    material = ElasticMaterial(young_modulus=210000., poisson_ratio=0.3)

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

    bc1 = BoundaryCondition(cp_index=np.array([0, 1, 2, 3]),
                            dof=np.array([0, 1, 2]),
                            value=0.
                            )
    model.add_boundary_condition(0, bc1)

    dload = DistributedLoad(el_index=np.array([0]),
                            dl_type=40,
                            magnitude=1000.)
    model.add_distributed_load(0, dload)

    model.refine_patch(ipatch=0,
                       nb_subdivision=np.array([1, 1, 1]),
                       nb_degree_elevation=np.array([1, 1, 1]),
                       additional_knots=[np.array([0.1,]),
                                         np.array([0.25,]),
                                         np.array([0.85])])

    stiff, rhs = model.build_stiffness_matrix()

    x = sp.linalg.spsolve(
        stiff[model.idof_free, :][:, model.idof_free],
        rhs[model.idof_free]
        )

    model.write_solution_vtu(x,
                             f"{tmp_path}/result.vtu",
                             data_flag=np.array([True, True, True]))


def test_linear_analysis_single_patch_3d_shell(tmp_path):
    """
    Build a 3D shel model
    Make a first refinement
    Modify geometry
    Make a second refinement
    Compute linear analysis solution
    Write VTU result file
    """
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
    # Deduce indices of control points NOT located at corners
    free_cps = np.setxor1d(model.cp_indices, corner_cps)
    # Elevate those control points along z axis
    model.cp_coordinates[free_cps, 2] += 0.5

    # Make a new refinement for analysis
    model.refine_patch(ipatch=0,
                       nb_subdivision=np.array([1, 1]))

    stiff, rhs = model.build_stiffness_matrix()

    x = sp.linalg.spsolve(
        stiff[model.idof_free, :][:, model.idof_free],
        rhs[model.idof_free]
        )

    model.write_control_mesh_vtu(f'{tmp_path}/control_net.vtu')
    model.write_solution_vtu(x,
                             f"{tmp_path}/result.vtu")


if __name__ == '__main__':
    # test_linear_analysis_single_patch_3d_solid('.')
    test_linear_analysis_single_patch_3d_shell('.')
