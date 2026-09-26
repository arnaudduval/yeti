from yeti_iga.pymfiga.common.material import J2General
from yeti_iga.pymfiga.common.physics import IncrementalElastoPlasticity
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel
import numpy as np

# Set global variables
TRACTION = 400.0
YOUNG, POISSON = 2500, 0.25
NBSTEPS = 101


def surface_force(args: dict):
    position = args["position"]
    x = position[0]
    nnz = np.size(position, axis=1)
    tmp = np.zeros((2, nnz))
    tmp[1] = x**2 - 1 / 4
    force = np.zeros((2, nnz))
    force[1] = -TRACTION * (np.min(tmp, axis=0)) ** 2
    return force


# Create geometry
degree, nbel = 2, 8
material = J2General(
    {
        "elastic_modulus": YOUNG,
        "elastic_limit": 1,
        "poisson_ratio": POISSON,
        "iso_hardening": {"name": "linear", "Eiso": 0.0},
        "kine_hardening": {"parameters": np.array([[500, 0]])},
    }
)
geometry = GeomdlGenerator(
    filename="square", geo_args={"degree": degree, "nbel": nbel}
).export_geometry()
patch = SinglePatch(geometry, quadclass="gs")
patch.generate()

# Set Dirichlet boundaries
boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX, DOF.UY))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MIN,
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

# Solve
model = MechanicalModel(material, patch, boundary)
time_list = np.linspace(0, 1.0, NBSTEPS)
force_ref = model.assemble_surface_force(
    {
        DOF.ALL: (
            {"direction": ParametricDirection.ETA, "face": BoundarySide.MAX},
            surface_force,
        )
    }
)
external_force = np.tensordot(time_list / time_list[-1], force_ref, axes=0)
displacement = np.zeros_like(external_force)
IncrementalElastoPlasticity().solve(model, displacement, external_force)
np.testing.assert_almost_equal(np.linalg.norm(displacement[1]), 0.00022436817518163766)
np.testing.assert_almost_equal(np.linalg.norm(displacement[-1]), 3.279200678425416)
