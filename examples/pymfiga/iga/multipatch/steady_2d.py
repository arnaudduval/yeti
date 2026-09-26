from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import SteadyHeatTransfer
from yeti_iga.pymfiga.common.io import IgaPostprocessing
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import (
    BoundaryCondition,
    ParametricDirection,
    BoundarySide,
    TargetPriority,
    DOF,
)
from yeti_iga.pymfiga.iga.single_model import ThermalModel as SpThermalModel
from yeti_iga.pymfiga.iga.multi_model import ThermalModel as MpThermalModel
from copy import deepcopy
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


def power_density(args: dict):
    position = args["position"]
    x = position[0]
    y = position[1]

    f = np.sin(np.pi * x) * np.sin(np.pi * y)
    return f


geo_args = {
    "degree": 3,
    "nbel": 8,
    "parameters": {"Rin": 1.0, "Rex": 2.0},
}
material = ThermalMaterial()
material.add_capacity(1, is_uniform=True)
material.add_conductivity(np.array([[1, 0.5], [0.5, 2]]), is_uniform=True, ndim=2)
geometry = GeomdlGenerator(filename="square", geo_args=geo_args).export_geometry()
patch1 = SinglePatch(geometry, quadclass="gs")
patch2 = deepcopy(patch1)
patch2.reflect(plane="xz")
patch1.generate()
patch2.generate()

# Set Dirichlet boundaries
boundary1 = BoundaryCondition(nbctrlpts=patch1.nbctrlpts, dofs=(DOF.T,))
boundary1.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MAX,
            "dofs": (DOF.T,),
        }
    ],
    constraint_type="dirichlet",
)
boundary1.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MIN,
            "target_priority": TargetPriority.SLAVE,
            "target_id": 2,
        }
    ],
    constraint_type="glue",
)
print(boundary1.select_nodes_for_gluing()[-1])

boundary2 = BoundaryCondition(nbctrlpts=patch2.nbctrlpts, dofs=(DOF.T,))
boundary2.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MAX,
            "dofs": (DOF.T,),
        }
    ],
    constraint_type="dirichlet",
)
boundary2.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MIN,
            "target_priority": TargetPriority.MASTER,
            "target_id": 1,
        }
    ],
    constraint_type="glue",
)

# Solve thermal problem
model1 = SpThermalModel(material, patch1, boundary1)
model2 = SpThermalModel(material, patch2, boundary2)
multipatchmodel = MpThermalModel({1: model1, 2: model2})

# Add forces
external_force_1 = model1.assemble_volumetric_force({DOF.ALL: power_density})
external_force_2 = np.zeros_like(external_force_1)
external_force = multipatchmodel.glue({1: external_force_1, 2: external_force_2})

# Solve
temperature = np.zeros_like(external_force)
SteadyHeatTransfer().solve(
    multipatchmodel,
    temperature,
    external_force,
)

# Postprocess
temperature_cutted = multipatchmodel.cut(temperature)
IgaPostprocessing.export_patch(
    patch1,
    "primal",
    filename="iga_multipatch_steady_1",
    fields={"temp": temperature_cutted[1]},
    folder=RESULT_FOLDER,
    sample_size=201,
)
IgaPostprocessing.export_patch(
    patch2,
    "primal",
    filename="iga_multipatch_steady_2",
    fields={"temp": temperature_cutted[2]},
    folder=RESULT_FOLDER,
    sample_size=201,
)
