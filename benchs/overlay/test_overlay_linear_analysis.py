# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test definition of an IGA model for linear elasticity analysis
"""

import numpy as np
import scipy.sparse as sp

from yeti_iga import IgaModel, Patch, ElasticMaterial, \
    BoundaryCondition, DistributedLoad

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
                       dof=0,
                       value=0.
                       )
bc2 = BoundaryCondition(cp_index=np.array([0, 1, 2, 3]),
                       dof=1,
                       value=0.
                       )
bc3 = BoundaryCondition(cp_index=np.array([0, 1, 2, 3]),
                       dof=2,
                       value=0.
                       )

model.add_boundary_condition(0, bc1)
model.add_boundary_condition(0, bc2)
model.add_boundary_condition(0, bc3)

dload = DistributedLoad(el_index=np.array([0]),
                        dl_type=40,
                        magnitude=1000.)

model.add_distributed_load(0, dload)


model.refine_patch(ipatch=0,
                   nb_refinements=np.array([1, 1, 1]),
                   nb_degree_elevation=np.array([1, 1, 1]),
                   additional_knots=[np.array([0.1, 0.2]),
                                     np.array([0.23, 0.24]),
                                     np.array([0.85, 0.95])])

# print(vars(model.iga_param))


stiff, rhs = model.build_stiffness_matrix()

x = sp.linalg.spsolve(
    stiff[model.idof_free, :][:, model.idof_free],
    rhs[model.idof_free]
    )
