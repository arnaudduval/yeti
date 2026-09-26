from yeti_iga.pymfiga.common.material import J2General
from yeti_iga.pymfiga.common.physics import IncrementalElastoPlasticity
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


def volume_force(args: dict):
    position = args["position"]
    force = 1e2 * np.ones_like(position)
    return force


# Global variables
DEGREE, NBEL = 2, 16
YOUNG, LENGTH = 1e6, 1
NBSTEPS = 21

# Define geometry
geometry = GeomdlGenerator(
    filename="line",
    geo_args={"degree": DEGREE, "nbel": NBEL, "parameters": {"L": LENGTH}},
).export_geometry()
patch = SinglePatch(geometry, quadclass="wq")
patch.generate()

# Define plasticity model
material = J2General(
    {
        "elastic_modulus": YOUNG,
        "elastic_limit": 1.0e1,
        "iso_hardening": {"name": "linear", "Eiso": YOUNG / 5},
    },
    is_unidimensional=True,
)

# Set boundary conditions
boundary = BoundaryCondition(patch.nbctrlpts, dofs=(DOF.UX,))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MIN,
            "dofs": (DOF.UX,),
        }
    ],
    constraint_type="dirichlet",
)

# Define model
model = MechanicalModel(material=material, patch=patch, boundary=boundary)

# Add external force
time_list = np.linspace(0, 1, NBSTEPS)
force_ref = model.assemble_volumetric_force({DOF.ALL: volume_force})
external_force = np.tensordot(time_list, force_ref, axes=0)

# Solve
displacement = np.zeros_like(external_force)
IncrementalElastoPlasticity().solve(model, displacement, external_force)
reference = [
    1.71611660e-05,
    4.91397481e-05,
    7.87745801e-05,
    1.06065662e-04,
    1.31012994e-04,
    1.53616578e-04,
    1.73876402e-04,
    1.91792514e-04,
    2.07364735e-04,
    2.20593731e-04,
    2.31477015e-04,
    2.40023875e-04,
    2.46199644e-04,
    2.50133703e-04,
    2.51343189e-04,
    2.51662869e-04,
    2.51680605e-04,
]
np.testing.assert_allclose(displacement[-1, 1:], reference)
