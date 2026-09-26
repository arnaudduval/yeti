from yeti_iga.pymfiga.common.material import J2General
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.model import MechanicalModel
from yeti_iga.pymfiga.common.physics.fem_contact import ContactProblem
from yeti_iga.pymfiga.common.io import FeaPostprocessing
import numpy as np

# Material definition
material = J2General({"elastic_modulus": 1e3, "poisson_ratio": 0.3})

# Time span definition
timespan = np.linspace(0.0, 1.0, 2)

# Define solid s1: master
#########################
# Geometry
mesh_s1 = MeshConverter(
    MeshpyGenerator(filename="square_modified", max_volume=0.005).mesh_geometry()
)


# Boundary conditions
def condition_disp_x(pts):
    return np.isclose(pts[0], 0.0)


def condition_disp_y(pts):
    return np.isclose(pts[1], 0.0)


def condition_s1_contact(pts):
    return np.isclose(pts[1], 1.0)


boundary_s1 = BoundaryCondition(mesh_s1, dofs=(DOF.UX, DOF.UY))
boundary_s1.add_constraint(
    {DOF.UX: condition_disp_x, DOF.UY: condition_disp_y},
    constraint_type="dirichlet",
)
boundary_s1.add_constraint([condition_s1_contact], constraint_type="contact")

# Problem definition
model_s1 = MechanicalModel(
    mesh_s1, material, FiniteElementQuadrature(boundary_order=2), boundary_s1
)

# Define solid s2: slave
#########################
# Geometry
mesh_s2 = MeshConverter(
    MeshpyGenerator(filename="semicircle", max_volume=0.0001).mesh_geometry()
)


# Boundary conditions
def condition_s2_contact(pts):
    return pts[1] < 1.1


boundary_s2 = BoundaryCondition(mesh_s2, dofs=(DOF.UX, DOF.UY))
boundary_s2.add_constraint({DOF.UX: condition_disp_x}, constraint_type="dirichlet")
boundary_s2.add_constraint([condition_s2_contact], constraint_type="contact")

# Problem definition
model_s2 = MechanicalModel(
    mesh_s2, material, FiniteElementQuadrature(boundary_order=2), boundary_s2
)


# Solve coupled model
# Compute external forces
def condition_force(pts):
    return np.isclose(pts[1], 1.2)


def surface_force(pts):
    force = -10 * np.ones(len(pts))
    return force


external_forces_s1 = np.zeros((len(timespan), model_s1.get_size_of_arrays()))
displacements_s1 = np.zeros_like(external_forces_s1)

external_forces_s2 = np.outer(
    timespan,
    model_s2.assemble_surface_force({DOF.UY: (condition_force, surface_force)}),
)
displacements_s2 = np.zeros_like(external_forces_s2)

# TODO: study why slow convergence (penalty ?)
ContactProblem(
    model_s1, model_s2, maxiters_nonlinear=20, tolerance_nonlinear=1e-3
).solve(
    displacements_s1,
    displacements_s2,
    external_forces_s1,
    external_forces_s2,
)


# Post processing
def postprocessing(model: MechanicalModel, displacement, filename="output"):
    strain = model.interpolate_strain(displacement=displacement, convert_to_3d=True)
    stress = material.eval_elastic_stress(strain)
    vonmises = material.eval_von_mises_stress(stress)
    displacement = np.reshape(displacement, (2, -1))
    FeaPostprocessing.export_meshpy(
        model.part,
        filename=filename,
        warping=displacement.T,
        point_data={
            "primal_disp_x": displacement[0],
            "primal_disp_y": displacement[1],
            "primal_VM": model.smoothing_dual_field_at_cell(vonmises[..., 0]),
        },
        cell_data={
            "dual_exx": strain[0, 0, :, 0],
            "dual_eyy": strain[1, 1, :, 0],
            "dual_txy": strain[0, 1, :, 0],
            "dual_VM": vonmises[..., 0],
        },
    )
    return


postprocessing(model_s1, displacements_s1[-1], filename="fem_contact_master")
postprocessing(model_s2, displacements_s2[-1], filename="fem_contact_slave")
