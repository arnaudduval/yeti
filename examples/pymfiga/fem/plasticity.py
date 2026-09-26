from yeti_iga.pymfiga.common.material import J2General
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.model import MechanicalModel
from yeti_iga.pymfiga.common.physics import IncrementalElastoPlasticity
import numpy as np

# Material definition
material = J2General(
    {
        "elastic_modulus": 1e3,
        "poisson_ratio": 0.3,
        "elastic_limit": 5,
        "iso_hardening": {"name": "linear", "Eiso": 5e2},
    }
)

# Geometry definition
mesh = MeshConverter(
    MeshpyGenerator(filename="plate_hole", max_volume=0.001).mesh_geometry()
)


# Boundary conditions
def condition_disp_x(pts):
    return np.isclose(pts[0], 0.0)


def condition_disp_y(pts):
    return np.isclose(pts[1], 0.0)


boundary = BoundaryCondition(mesh, dofs=(DOF.UX, DOF.UY))
boundary.add_constraint(
    {DOF.UX: condition_disp_x, DOF.UY: condition_disp_y},
    constraint_type="dirichlet",
)

# Model definition
model = MechanicalModel(
    mesh,
    material,
    FiniteElementQuadrature(interior_order=1, boundary_order=1),
    boundary,
)


# Computation of external force
def condition_force(pts):
    return np.isclose(pts[1], 1.0)


def surface_force(pts):
    force = 10 * np.ones(len(pts))
    return force


timespan = np.linspace(0, 1, 21)
ref_force = model.assemble_surface_force({DOF.UY: (condition_force, surface_force)})
external_force = np.outer(timespan, ref_force)
displacement = np.zeros_like(external_force)

# Solve
IncrementalElastoPlasticity().solve(model, displacement, external_force)
