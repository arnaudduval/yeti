from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel
from yeti_iga.pymfiga.common.physics import StaticElastoPlasticity
from yeti_iga.pymfiga.common.io import IgaPostprocessing
from time import time
import numpy as np

# Set global variables
DEGREE, NBEL = 4, 16
TRACTION, RINT, REXT = 1.0, 1.0, 2.0
YOUNG, POISSON = 1e3, 0.25


def surface_force(args: dict):
    position = args["position"]
    x = position[0]
    y = position[1]

    r_square = x**2 + y**2
    b = RINT**2 / r_square
    b2 = b**2

    r = np.sqrt(r_square)
    cos_theta = x / r
    sin_theta = y / r

    cos_3theta = 4 * cos_theta**3 - 3 * cos_theta
    sin_3theta = 3 * sin_theta - 4 * sin_theta**3

    force = np.zeros_like(position)
    force[0] = (
        TRACTION
        / 2
        * (2 * cos_theta - b * (2 * cos_theta + 3 * cos_3theta) + 3 * b2 * cos_3theta)
    )
    force[1] = TRACTION / 2 * 3 * sin_3theta * (b2 - b)
    return force


# Create geometry
geo_args = {
    "degree": DEGREE,
    "nbel": NBEL,
    "parameters": {"Rin": RINT, "Rex": REXT},
}
geometry = GeomdlGenerator(
    filename="quarter_annulus", geo_args=geo_args
).export_geometry()
patch = SinglePatch(geometry, quadclass="gs")
patch.generate()

# Set material
material = LinearElasticity({"elastic_modulus": YOUNG, "poisson_ratio": POISSON})

# Set Dirichlet boundaries
boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX, DOF.UY))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MAX,
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

# Define model
model = MechanicalModel(material, patch, boundary)

external_force = model.assemble_surface_force(
    {
        DOF.ALL: (
            {"direction": ParametricDirection.XI, "face": BoundarySide.MAX},
            surface_force,
        )
    }
)

# Solve problem as nonlinear
displacement = np.zeros_like(external_force)
model.compute_mf_stiffness(displacement)
start = time()
StaticElastoPlasticity().solve(model, displacement, external_force)
finish = time()
print("Problem solved in %.2f seconds" % (finish - start))

# Post processing
strain_3d = model.interpolate_strain(displacement, convert_to_3d=True)
stress_3d = material.eval_elastic_stress(strain_3d)
vmstress = material.eval_von_mises_stress(stress_3d)
np.testing.assert_almost_equal(vmstress.max(), 2.679077396499368)
np.testing.assert_almost_equal(vmstress.min(), 0.030932057760478505)

IgaPostprocessing.export_patch(
    model.part, "dual", fields={"vms": vmstress}, filename="iga_elasticity"
)


# Extra verification
def surface_force_x(args: dict):
    f = surface_force(args)
    return f[0]


def surface_force_y(args: dict):
    f = surface_force(args)
    return f[1]


constraint_nodes = model.get_free_and_constraint_nodes()[-1]
external_force_2 = model.assemble_surface_force(
    {
        DOF.UX: (
            {"direction": ParametricDirection.XI, "face": BoundarySide.MAX},
            surface_force_x,
        ),
        DOF.UY: (
            {"direction": ParametricDirection.XI, "face": BoundarySide.MAX},
            surface_force_y,
        ),
    }
)
external_force[constraint_nodes] = 0.0
external_force_2[constraint_nodes] = 0.0
err = np.linalg.norm(external_force - external_force_2) / np.linalg.norm(
    external_force_2
)
np.testing.assert_almost_equal(err, 0.0)
