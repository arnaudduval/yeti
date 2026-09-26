from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.model import MechanicalModel
from yeti_iga.pymfiga.common.io import FeaPostprocessing
import numpy as np

TX, RINT = 1, 0.2

# Material definition
material = LinearElasticity({"elastic_modulus": 1e3, "poisson_ratio": 0.25})

# Geometry definition
mesh = MeshConverter(
    MeshpyGenerator(filename="plate_hole", max_volume=0.0005).mesh_geometry()
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
    mesh, material, FiniteElementQuadrature(interior_order=2), boundary
)


# Computation of external force
def stress_tensor(x, y):
    r_square = x**2 + y**2
    r_fourth = r_square**2
    r = np.sqrt(r_square)
    cos_theta, sin_theta = x / r, y / r
    rotation = np.zeros((2, 2, *np.shape(x)))
    rotation[0, 0, ...] = cos_theta
    rotation[1, 1, ...] = cos_theta
    rotation[0, 1, ...] = -sin_theta
    rotation[1, 0, ...] = sin_theta

    cos_2theta = cos_theta**2 - sin_theta**2
    sin_2theta = 2 * cos_theta * sin_theta
    polar_tensor = np.zeros_like(rotation)
    polar_tensor[0, 0, ...] = (
        TX / 2 * (1 - RINT**2 / r_square)
        + TX / 2 * (1 - 4 * RINT**2 / r_square + 3 * RINT**4 / r_fourth) * cos_2theta
    )
    polar_tensor[1, 1, ...] = (
        TX / 2 * (1 + RINT**2 / r_square)
        - TX / 2 * (1 + 3 * RINT**2 / r_fourth) * cos_2theta
    )
    polar_tensor[0, 1, ...] = (
        -TX / 2 * (1 + 2 * RINT**2 / r_square - 3 * RINT**4 / r_fourth) * sin_2theta
    )
    polar_tensor[1, 0, ...] = np.copy(polar_tensor[0, 1, ...])
    return np.einsum("il...,lm...,jm...->ij...", rotation, polar_tensor, rotation)


def surface_force_single(position):
    x = position[0]
    return np.ones_like(x)


def surface_force_all(position):
    x = position[0]
    f = np.zeros((2, *np.shape(x)))
    f[0] = 1.0
    return f


def surface_force_1(position):
    x = position[0]
    y = position[1]
    stress = stress_tensor(x, y)
    normal = np.zeros((2, *np.shape(x)))
    normal[0, ...] = 1.0
    return np.einsum("ij...,j...->i...", stress, normal)


def surface_force_2(position):
    x = position[0]
    y = position[1]
    stress = stress_tensor(x, y)
    normal = np.zeros((2, *np.shape(x)))
    normal[1, ...] = 1.0
    return np.einsum("ij...,j...->i...", stress, normal)


def condition_force_1(pts):
    return np.isclose(pts[0], 1.0)


def condition_force_2(pts):
    return np.isclose(pts[1], 1.0)


force_sing = model.assemble_volumetric_force({DOF.UX: surface_force_single})
force_all = model.assemble_volumetric_force({DOF.ALL: surface_force_all})
err = np.linalg.norm(force_all - force_sing) / np.linalg.norm(force_sing)
np.testing.assert_almost_equal(err, 0.0)

force_sing = model.assemble_surface_force(
    {DOF.UX: (condition_force_1, surface_force_single)}
)
force_all = model.assemble_surface_force(
    {DOF.ALL: (condition_force_1, surface_force_all)}
)
err = np.linalg.norm(force_all - force_sing) / np.linalg.norm(force_sing)
np.testing.assert_almost_equal(err, 0.0)

force = model.assemble_surface_force({DOF.ALL: (condition_force_1, surface_force_1)})
force += model.assemble_surface_force({DOF.ALL: (condition_force_2, surface_force_2)})

displacement = model.solve_linearized_system(force)

# Postprocessing
strain_at_quadpts_cart = model.interpolate_strain(
    displacement=displacement, convert_to_3d=True
)
stress_at_quadpts_cart = material.eval_elastic_stress(strain_at_quadpts_cart)
vm_stress_at_quadpts_cart = material.eval_von_mises_stress(stress_at_quadpts_cart)

stress_at_cell_cart = np.mean(stress_at_quadpts_cart, axis=-1)
vm_stress_at_cell_cart = np.mean(vm_stress_at_quadpts_cart, axis=-1)

vm_stress_max = np.max(vm_stress_at_cell_cart)
np.testing.assert_allclose(vm_stress_max, 2.858481744362727)

FeaPostprocessing.export_meshpy(
    mesh,
    filename="fem_elasticity_cartesian",
    cell_data={
        "dual_Sxx": stress_at_cell_cart[0][0],
        "dual_Syy": stress_at_cell_cart[1][1],
        "dual_Sxy": stress_at_cell_cart[0][1],
        "dual_VM": vm_stress_at_cell_cart,
    },
)
