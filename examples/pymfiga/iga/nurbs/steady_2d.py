from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import SteadyHeatTransfer
from yeti_iga.pymfiga.common.io import IgaPostprocessing, Postprocessing
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, DOF, ParametricDirection, BoundarySide
from yeti_iga.pymfiga.iga.single_model import ThermalModel
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


def power_density(args: dict):
    position = args["position"]
    x = position[0]
    y = position[1]

    f = (
        x
        * y
        * (x + 4)
        * (y - 4)
        * np.sin(np.pi * x)
        * np.sin(np.pi * y)
        * (x**2 + y**2 - 1)
    )
    return f


DEGREE, NBEL = 4, 32
geo_args = {
    "degree": DEGREE,
    "nbel": (NBEL, 2 * NBEL),
}

material = ThermalMaterial()
material.add_capacity(1, is_uniform=True)
material.add_conductivity(np.array([[1, 0.5], [0.5, 2]]), is_uniform=True, ndim=2)
geometry = GeomdlGenerator(
    filename="nurbs_plate_hole", geo_args=geo_args
).export_geometry()
patch = SinglePatch(geometry, quadclass="gs").generate()

# Set Dirichlet boundaries
boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.T,))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": (ParametricDirection.XI, ParametricDirection.ETA),
            "face": BoundarySide.BOTH,
            "dofs": (DOF.T,),
        }
    ],
    constraint_type="dirichlet",
)

# Solve heat problem
model = ThermalModel(material, patch, boundary)
external_force = model.assemble_volumetric_force({DOF.ALL: power_density})
temperature = np.zeros_like(external_force)
SteadyHeatTransfer().solve(
    model,
    temperature,
    external_force,
)

IgaPostprocessing.export_patch(
    model.part,
    "primal",
    fields={"temp": temperature},
    filename="iga_nurbs_steady",
    folder=RESULT_FOLDER,
    sample_size=np.array([201, 200]),
)

Postprocessing.vtk2png(
    filename="iga_nurbs_steady",
    folder=RESULT_FOLDER,
    fieldname="temp",
    cmap="coolwarm",
    format="vts",
    cpos="xy",
    n_colors=11,
    scalar_bar_args={
        "title": "Temperature\n",
        "n_labels": 3,
        "height": 0.05,
    },
)
