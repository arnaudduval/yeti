from yeti_iga.pymfiga.common.material import J2General
from yeti_iga.pymfiga.common.physics import StaticElastoPlasticity
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel
import numpy as np


def surface_force(args):
    position = args["position"]
    prop = np.zeros_like(position)
    prop[1] = 1.0
    return prop


# Create model
degree, nbel, length = 2, 32, 1.0
geometry = GeomdlGenerator(
    filename="square", geo_args={"degree": degree, "nbel": nbel}
).export_geometry()
patch = SinglePatch(geometry, quadclass="gs")
patch.generate()

# Add material
material = J2General({"elastic_modulus": 1e2, "poisson_ratio": 0.3})

# Set Dirichlet boundaries
boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX, DOF.UY))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.BOTH,
            "dofs": (DOF.UX,),
        },
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MIN,
            "dofs": (DOF.UY,),
        },
    ],
    constraint_type="dirichlet",
)

# Elasticity model
model = MechanicalModel(material, patch, boundary)
external_force = model.assemble_surface_force(
    {
        DOF.ALL: (
            {"direction": ParametricDirection.ETA, "face": BoundarySide.MAX},
            surface_force,
        )
    }
)

# Solve problem as linear
contraint_nodes = model.get_free_and_constraint_nodes()[-1]
linear_solver = LinearSolver(
    tolerance=1e-8,
    maxiters=100,
    cleandod=contraint_nodes,
)

output = model.solve_linearized_system(
    external_force,
    linear_solver_backend=linear_solver,
)

# Solve problem as nonlinear
displacement = np.zeros_like(external_force)
StaticElastoPlasticity().solve(model, displacement, external_force)

error = np.linalg.norm(output - displacement) / np.linalg.norm(displacement)
np.testing.assert_almost_equal(error, 2.566345e-13)
np.testing.assert_almost_equal(displacement.max(), 0.007428571428570144)

# Post processing
strain_3d = model.interpolate_strain(displacement, convert_to_3d=True)

stress_3d = material.eval_elastic_stress(strain_3d)
vonmises = material.eval_von_mises_stress(stress_3d)
error = abs(vonmises.max() - vonmises.min())
np.testing.assert_almost_equal(vonmises.min(), 0.5714285714282623)
np.testing.assert_almost_equal(error, 0.0)
print("Von misses max:%.4e, min:%.4e" % (vonmises.max(), vonmises.min()))
print("Difference: %.4e" % (error))
