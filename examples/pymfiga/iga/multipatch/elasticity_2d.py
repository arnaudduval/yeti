from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.physics import StaticElastoPlasticity
from yeti_iga.pymfiga.common.io import IgaPostprocessing
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import (
    BoundaryCondition,
    ParametricDirection,
    BoundarySide,
    TargetPriority,
    DOF,
)
from yeti_iga.pymfiga.iga.single_model import MechanicalModel as SpMechanicalModel
from yeti_iga.pymfiga.iga.multi_model import MechanicalModel as MpMechanicalModel
from copy import deepcopy
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


def force(args: dict):
    position = args["position"]
    x = position[0]
    f = np.zeros((2, len(x)))
    f[1] = 1
    return f


geo_args = {
    "degree": 5,
    "nbel": 16,
    "parameters": {"Rin": 1.0, "Rex": 2.0},
}
material = LinearElasticity({"elastic_modulus": 1e9, "poisson_ratio": 0.3})
geometry = GeomdlGenerator(
    filename="nurbs_quarter_annulus", geo_args=geo_args
).export_geometry()
patch1 = SinglePatch(geometry, quadclass="wq")
patch2 = deepcopy(patch1)
patch2.reflect(plane="xz")
patch1.generate()
patch2.generate()

# Set Dirichlet boundaries
boundary1 = BoundaryCondition(nbctrlpts=patch1.nbctrlpts, dofs=(DOF.UX, DOF.UY))
boundary1.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MAX,
            "dofs": (
                DOF.UX,
                DOF.UY,
            ),
        },
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

boundary2 = BoundaryCondition(nbctrlpts=patch2.nbctrlpts, dofs=(DOF.UX, DOF.UY))
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
model1 = SpMechanicalModel(material, patch1, boundary1)
model2 = SpMechanicalModel(material, patch2, boundary2)
multipatchmodel = MpMechanicalModel({1: model1, 2: model2})

# Add forces
external_force_1 = model1.assemble_volumetric_force({DOF.ALL: force})
external_force_2 = np.zeros_like(external_force_1)
external_force = multipatchmodel.glue({1: external_force_1, 2: external_force_2})

# Solve
displacement = np.zeros_like(external_force)
StaticElastoPlasticity().solve(
    multipatchmodel,
    displacement,
    external_force,
)

# Postprocess
displacement_cutted = multipatchmodel.cut(displacement)
IgaPostprocessing.export_patch(
    patch1,
    "primal",
    filename="iga_multipatch_elasticity_1",
    fields={"disp": np.reshape(displacement_cutted[1], (2, -1))},
    folder=RESULT_FOLDER,
    sample_size=201,
)
IgaPostprocessing.export_patch(
    patch2,
    "primal",
    filename="iga_multipatch_elasticity_2",
    fields={"disp": np.reshape(displacement_cutted[2], (2, -1))},
    folder=RESULT_FOLDER,
    sample_size=201,
)
