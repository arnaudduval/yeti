from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.physics import StaticElastoPlasticity
from yeti_iga.pymfiga.common.io import IgaPostprocessing, Postprocessing
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.single_model import MechanicalModel
from time import time
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


# Set global variables
DEGREE, NBEL = 4, 16
YOUNG, POISSON = 1e3, 0.3


def surface_force(args: dict):
    position = args["position"]
    y = position[1]
    # FOR PLATE
    mask = np.where(np.abs(y - 4.0) < 1e-8)
    force = np.zeros_like(position)
    force[1, mask] = 1.0
    return force


# Create geometry
geo_args = {
    "degree": DEGREE,
    "nbel": (NBEL, 2 * NBEL),
}
material = LinearElasticity({"elastic_modulus": YOUNG, "poisson_ratio": POISSON})
geometry = GeomdlGenerator(
    filename="nurbs_plate_hole", geo_args=geo_args
).export_geometry()
patch = SinglePatch(geometry, quadclass="gs")
patch.generate()

# Set Dirichlet boundaries
boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.UX, DOF.UY))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MAX,
            "dofs": (DOF.UX,),
        },
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MIN,
            "dofs": (DOF.UY,),
        },
    ],
    constraint_type="dirichlet",
)

# Solve elastic problem
model = MechanicalModel(material, patch, boundary)
external_force = model.assemble_surface_force(
    {
        DOF.ALL: (
            {"direction": ParametricDirection.ETA, "face": BoundarySide.MAX},
            surface_force,
        )
    }
)
displacement = np.zeros_like(external_force)

start = time()
StaticElastoPlasticity().solve(model, displacement, external_force)
finish = time()
print("Problem solved in %.2f seconds" % (finish - start))

strain_3d = model.interpolate_strain(displacement, convert_to_3d=True)
stress_3d = material.eval_elastic_stress(strain_3d)
vmstress = material.eval_von_mises_stress(stress_3d)
np.testing.assert_almost_equal(vmstress.max(), 3.1385213894458874)
np.testing.assert_almost_equal(vmstress.min(), 0.06228192216655326)

IgaPostprocessing.export_patch(
    patch,
    "dual",
    fields={"vms": vmstress},
    filename="iga_nurbs_elasticity",
    folder=RESULT_FOLDER,
)

Postprocessing.vtk2png(
    filename="iga_nurbs_elasticity",
    folder=RESULT_FOLDER,
    fieldname="vms",
    cmap="coolwarm",
    format="vts",
    cpos="xy",
    n_colors=11,
    scalar_bar_args={
        "title": "Von Mises\n",
        "n_labels": 3,
        "height": 0.05,
    },
)
