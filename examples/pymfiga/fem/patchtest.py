from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.model import MechanicalModel
from yeti_iga.pymfiga.common.io import FeaPostprocessing
import numpy as np

# Material definition
material = LinearElasticity({"elastic_modulus": 1e2, "poisson_ratio": 0.3})

# Geometry definition
mesh = MeshConverter(
    MeshpyGenerator(filename="square", max_volume=0.01).mesh_geometry()
)
mesh.convert_linear_to_quadratic()


# Boundary condition
def condition_disp_x(pts):
    return np.isclose(pts[0], 0.0) or np.isclose(pts[0], 1.0)


def condition_disp_y(pts):
    return np.isclose(pts[1], 0.0)


boundary = BoundaryCondition(
    mesh, dofs=(DOF.UX, DOF.UY), dofs_funspace=(("lagrange", 2), ("lagrange", 2))
)
boundary.add_constraint(
    {DOF.UX: condition_disp_x, DOF.UY: condition_disp_y},
    constraint_type="dirichlet",
)

# Model definition
model = MechanicalModel(
    mesh,
    material,
    FiniteElementQuadrature(interior_order=2, boundary_order=2),
    boundary,
    lagrange_order="quadratic",
)


# Computation of external force
def condition_force(pts):
    return np.isclose(pts[1], 1.0)


def surface_force(pts):
    return np.ones(len(pts[0]))


force = model.assemble_surface_force({DOF.UY: (condition_force, surface_force)})
displacement = model.solve_linearized_system(force)
np.testing.assert_almost_equal(displacement.max(), 0.00742857142857145)

# Postprocessing
strain_at_quadpts = model.interpolate_strain(displacement, convert_to_3d=True)
stress_at_quadpts = material.eval_elastic_stress(strain_at_quadpts)
vm_stress_at_quadpts = material.eval_von_mises_stress(stress_at_quadpts)
displacement = np.reshape(displacement, (2, -1))

strain_at_cell = np.mean(strain_at_quadpts, axis=-1)
stress_at_cell = np.mean(stress_at_quadpts, axis=-1)
vm_stress_at_cell = np.mean(vm_stress_at_quadpts, axis=-1)
np.testing.assert_allclose(vm_stress_at_cell.max(), 0.5714285714285812)
np.testing.assert_allclose(vm_stress_at_cell.min(), 0.5714285714285812)


FeaPostprocessing.export_meshpy(
    mesh,
    filename="fem_patchtest",
    point_data={
        "primal_disp_x": displacement[0],
        "primal_disp_y": displacement[1],
        "primal_VM": model.smoothing_dual_field_at_cell(vm_stress_at_cell),
    },
    mesh_order="quadratic",
    cell_data={
        "dual_exx": strain_at_cell[0, 0, :],
        "dual_eyy": strain_at_cell[1, 1, :],
        "dual_exy": strain_at_cell[0, 1, :],
        "dual_sxx": strain_at_cell[0, 0, :],
        "dual_syy": strain_at_cell[1, 1, :],
        "dual_VM": vm_stress_at_cell,
    },
)
