from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.model import MechanicalModel
from yeti_iga.pymfiga.common.io import FeaPostprocessing
from yeti_iga.pymfiga.common.numerics.coordinates.coordinate_transformation import (
    Cartesian2Polar,
    Cartesian2Cylindrical,
)
import numpy as np

TX, RINT, REXT = 1.0, 0.2, 1.0

# Material definition
material = LinearElasticity({"elastic_modulus": 1e3, "poisson_ratio": 0.25})

# Geometry definition
mesh = MeshConverter(
    MeshpyGenerator(filename="quarter_annulus", max_volume=0.0005).mesh_geometry()
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
    mesh, material, FiniteElementQuadrature(interior_order=1), boundary
)


# Computation of external force
def stress_tensor_polar(r, th):
    r_square, r_third, r_fourth = r**2, r**3, r**4
    cos_theta, sin_theta = np.cos(th), np.sin(th)
    cos_2theta = cos_theta**2 - sin_theta**2
    sin_2theta = 2 * cos_theta * sin_theta
    polar_tensor = np.zeros((2, 2, *np.shape(r)))
    polar_tensor[0, 0, ...] = (1 - RINT**2 / r_square) + (
        1 - 4 * RINT**2 / r_square + 3 * RINT**4 / r_fourth
    ) * cos_2theta
    polar_tensor[0, 1, ...] = (
        -(1 + 2 * RINT**2 / r - 3 * RINT**4 / r_third) * sin_2theta
    )
    polar_tensor[1, 1, ...] = (1 + RINT**2) - (
        r_square + 3 * RINT**2 / r_square
    ) * cos_2theta

    polar_tensor[1, 0, ...] = np.copy(polar_tensor[0, 1, ...])
    return TX * polar_tensor / 2


def surface_force(position):
    cart2pol = Cartesian2Polar()
    x = position[0]
    y = position[1]
    r = np.sqrt(x**2 + y**2)
    th = np.arctan2(y, x)
    stress_polar = stress_tensor_polar(r, th)
    normal_polar = np.zeros_like(position)
    normal_polar[0, ...] = 1.0
    force_polar = np.einsum("ij...,j...->i...", stress_polar, normal_polar)
    force_carte = np.zeros_like(force_polar)
    for i, _ in enumerate(x):
        force_carte[:, i] = cart2pol.transform_vector(
            vector=force_polar[:, i],
            from_base="polar",
            point=[r[i], th[i]],
            point_base="polar",
            comp_type_in="upper",
        )
    return force_carte


def condition_force(pts):
    radius = np.sqrt(pts[0] ** 2 + pts[1] ** 2)
    return np.isclose(radius, 1.0, rtol=1e-2)


force = model.assemble_surface_force({DOF.ALL: (condition_force, surface_force)})

displacement = model.solve_linearized_system(force)

#  Postprocessing
strain_at_quadpts_cart = model.interpolate_strain(
    displacement=displacement, convert_to_3d=True
)
stress_at_quadpts_cart = material.eval_elastic_stress(strain_at_quadpts_cart)
vm_stress_at_quadpts_cart = material.eval_von_mises_stress(stress_at_quadpts_cart)

stress_at_cell_cart = np.mean(stress_at_quadpts_cart, axis=-1)
vm_stress_at_cell_cart = np.mean(vm_stress_at_quadpts_cart, axis=-1)
vm_stress_max = np.max(vm_stress_at_cell_cart)
np.testing.assert_allclose(vm_stress_max, 2.7345531726459744)


# ========================================
# Change of system of reference
# ========================================
def normalization(tensor, position):
    radius = np.linalg.norm(position)
    tensor[0, 1, ...] /= radius
    tensor[1, 0, ...] /= radius
    tensor[1, 1, ...] /= radius**2
    return tensor


def compute_quadrature_points(model: MechanicalModel):
    order = model.part.mesh_order
    output = model.prepare_lagrange_integral_on_element(order)
    funbasis, _, elements, points = output
    quadrature_point = [np.array([]) for _ in range(len(elements))]
    for el, idx_nodes_element in enumerate(elements):
        elem_coords = np.array([points[idx_nd] for idx_nd in idx_nodes_element])
        evalpoints = model.evaluate_field(funbasis, elem_coords)
        quadrature_point[el] = evalpoints
    out = np.insert(np.asarray(quadrature_point), 2, 0, axis=1)
    out = np.moveaxis(out, 0, 1)
    return out


cart2cyl = Cartesian2Cylindrical()
quadpts_at_elements = compute_quadrature_points(model)
stress_at_quadpts_cyl = np.zeros_like(stress_at_quadpts_cart)

for el in range(quadpts_at_elements.shape[1]):
    for qp in range(quadpts_at_elements.shape[2]):
        stress_cyl = cart2cyl.transform_tensor(
            tensor=stress_at_quadpts_cart[:, :, el, qp],
            from_base="cartesian",
            point=quadpts_at_elements[:, el, qp],
            point_base="cartesian",
            comp_type_in=["lower"] * 2,
            tensor_rank=2,
        )
        stress_cyl = normalization(stress_cyl, quadpts_at_elements[:, el, qp])
        stress_at_quadpts_cyl[:, :, el, qp] = stress_cyl

vm_stress_at_quadpts_cyl = material.eval_von_mises_stress(stress_at_quadpts_cyl)

stress_at_cell_cyl = np.mean(stress_at_quadpts_cyl, axis=-1)
vm_stress_at_cell_cyl = np.mean(vm_stress_at_quadpts_cyl, axis=-1)

FeaPostprocessing.export_meshpy(
    mesh,
    filename="fem_elasticity_polar",
    cell_data={
        "dual_Srr": stress_at_cell_cyl[0][0],
        "dual_Stt": stress_at_cell_cyl[1][1],
        "dual_Srt": stress_at_cell_cyl[0][1],
        "dual_VM": vm_stress_at_cell_cyl,
    },
)

error = np.linalg.norm(
    vm_stress_at_quadpts_cyl - vm_stress_at_quadpts_cart
) / np.linalg.norm(vm_stress_at_quadpts_cart)

np.testing.assert_almost_equal(error, 0.0)
print(f"Error in transformation {error}")
