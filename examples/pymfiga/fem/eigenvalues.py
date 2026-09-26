from yeti_iga.pymfiga.common.material import ThermalMaterial, LinearElasticity
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.model import MechanicalModel, ThermalModel
from yeti_iga.pymfiga.common.physics import EigenProblem
import numpy as np

# Eigenvalues in a thermal model
material = ThermalMaterial()
material.add_conductivity(inpt=1, is_uniform=True, ndim=2)
material.add_capacity(inpt=1, is_uniform=True)

mesh = MeshConverter(
    MeshpyGenerator(filename="square", max_volume=0.01).mesh_geometry()
)

boundary = BoundaryCondition(mesh, dofs=(DOF.T,))
boundary.add_constraint({DOF.T: True}, constraint_type="dirichlet")

model = ThermalModel(mesh, material, FiniteElementQuadrature(), boundary)
eigenvalues, eigenvectors = EigenProblem().solve(model, k=5)

true_values = [20.1926543, 52.16350324, 52.44877112, 87.02293159, 109.4249949]
np.testing.assert_allclose(eigenvalues, true_values)

################################################################################

# Eigenvalues in a mechanical model
material = LinearElasticity({"elastic_modulus": 1e3, "poisson_ratio": 0.25})
material.add_density(10, is_uniform=True)


def condition_disp_x(pts):
    return np.isclose(pts[0], 0.0)


def condition_disp_y(pts):
    return np.isclose(pts[1], 0.0)


boundary = BoundaryCondition(mesh, dofs=(DOF.UX, DOF.UY))
boundary.add_constraint(
    {DOF.UX: condition_disp_x, DOF.UY: condition_disp_y},
    constraint_type="dirichlet",
)

model = MechanicalModel(mesh, material, FiniteElementQuadrature(), boundary)
eigenvalues, eigenvectors = EigenProblem().solve(model, k=5)

true_values = [199.02221027, 267.90712612, 363.91169082, 769.54722346, 1194.00440299]
np.testing.assert_allclose(eigenvalues, true_values)
