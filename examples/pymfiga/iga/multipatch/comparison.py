from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import SteadyHeatTransfer
from yeti_iga.pymfiga.common.io import IgaPostprocessing
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from yeti_iga.pymfiga.iga.geometry.primitives import make_BSPLINE_line, create_NURBS_arc
from yeti_iga.pymfiga.iga.boundary import (
    BoundaryCondition,
    ParametricDirection,
    BoundarySide,
    DOF,
    TargetPriority,
)
from yeti_iga.pymfiga.iga.single_model import ThermalModel as SpThermalModel
from yeti_iga.pymfiga.iga.multi_model import ThermalModel as MpThermalModel
from yeti_iga.pymfiga.iga.norms import SpaceNormSinglePatch
from geomdl import NURBS
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
        3
        * np.pi**2
        * np.sin(np.pi * x)
        * np.sin(np.pi * y)
        * (x**2 + y**2 - 1)
        * (x**2 + y**2 - 4)
        - 16 * y**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
        - 6 * np.sin(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 1)
        - 6 * np.sin(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 4)
        - 8 * x * y * np.sin(np.pi * x) * np.sin(np.pi * y)
        - np.pi**2
        * np.cos(np.pi * x)
        * np.cos(np.pi * y)
        * (x**2 + y**2 - 1)
        * (x**2 + y**2 - 4)
        - 4 * x * np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 1)
        - 2 * x * np.pi * np.cos(np.pi * y) * np.sin(np.pi * x) * (x**2 + y**2 - 1)
        - 4 * x * np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 4)
        - 2 * x * np.pi * np.cos(np.pi * y) * np.sin(np.pi * x) * (x**2 + y**2 - 4)
        - 2 * y * np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 1)
        - 8 * y * np.pi * np.cos(np.pi * y) * np.sin(np.pi * x) * (x**2 + y**2 - 1)
        - 2 * y * np.pi * np.cos(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 4)
        - 8 * y * np.pi * np.cos(np.pi * y) * np.sin(np.pi * x) * (x**2 + y**2 - 4)
        - 8 * x**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    )
    return f


def exact_temperature(args: dict):
    position = args["position"]
    x = position[0]
    y = position[1]
    u = (
        np.sin(np.pi * x)
        * np.sin(np.pi * y)
        * (x**2 + y**2 - 1.0)
        * (x**2 + y**2 - 4.0)
    )
    return u


def create_nurbs_patch(
    degree_u: int, degree_v: int, nbel_u: int, nbel_v: int, **geo_args
) -> NURBS.Surface:

    assert degree_u > 1 and degree_v > 0

    RINT = geo_args["radius_int"]
    REXT = geo_args["radius_ext"]
    THET = geo_args["angle_extension"]

    # Construction of the arc
    obj_arc = create_NURBS_arc(degree_v, nbel_v, alpha_ini=THET[0], alpha_end=THET[1])
    knotvector_v = obj_arc.knotvector
    ctrlpts_arc = obj_arc.ctrlpts
    weights_arc = obj_arc.weights

    # Construction of line
    knotvector_u, ctrlpts_line = make_BSPLINE_line(degree_u, nbel_u)
    ctrlpts_line = RINT + ctrlpts_line * (REXT - RINT)

    # Construction of annulus sector
    ctrlpts = [
        [x_line * x_arc * w_arc, x_line * y_arc * w_arc, 0.0, w_arc]
        for x_line in ctrlpts_line
        for (x_arc, y_arc, _), w_arc in zip(ctrlpts_arc, weights_arc)
    ]

    # Create surface
    obj = NURBS.Surface()
    obj.degree_u = degree_u
    obj.degree_v = degree_v
    obj.set_ctrlpts(ctrlpts, len(ctrlpts_line), len(ctrlpts_arc))
    obj.knotvector_u = knotvector_u
    obj.knotvector_v = knotvector_v

    return obj


# Define material
material = ThermalMaterial()
material.add_capacity(1, is_uniform=True)
material.add_conductivity(np.array([[1, 0.5], [0.5, 2]]), is_uniform=True, ndim=2)

# ################################################################################
# # SINGLE PATCH
# ################################################################################
DEGREE = 4

# Define geometry
nbel = 16
geo_args = {"radius_int": 1.0, "radius_ext": 2.0, "angle_extension": (0, np.pi / 2)}
geometry = create_nurbs_patch(
    degree_u=DEGREE, degree_v=DEGREE, nbel_u=nbel, nbel_v=nbel, **geo_args
)
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

# Solve heat transfer problem
model = SpThermalModel(material, patch, boundary)
external_force = model.assemble_volumetric_force({DOF.ALL: power_density})
temperature = np.zeros_like(external_force)
SteadyHeatTransfer().solve(model, temperature, external_force)

relerror = SpaceNormSinglePatch(
    model, "l2", {"exact_function": exact_temperature}
).eval(temperature)[-1]
print(f"Error of single patch: {relerror:.2e}")

IgaPostprocessing.plot_patch(
    patch_list=[patch],
    filename="iga_comparison_singlepatch",
    folder=RESULT_FOLDER,
    primal_list=[temperature],
    add_colorbar=True,
)

################################################################################
# MULTI PATCH
################################################################################
nbel = 16

geo_args = {
    "radius_int": 1.0,
    "radius_ext": 1.5,
    "angle_extension": (0.0, np.pi / 4),
}
geometry = create_nurbs_patch(
    degree_u=DEGREE, degree_v=DEGREE, nbel_u=nbel, nbel_v=nbel, **geo_args
)
patch_0 = SinglePatch(geometry, quadclass="gs").generate()

geo_args = {
    "radius_int": 1.0,
    "radius_ext": 1.5,
    "angle_extension": (np.pi / 4, np.pi / 2),
}
geometry = create_nurbs_patch(
    degree_u=DEGREE, degree_v=DEGREE, nbel_u=nbel, nbel_v=nbel, **geo_args
)
patch_1 = SinglePatch(geometry, quadclass="gs").generate()

geo_args = {
    "radius_int": 1.5,
    "radius_ext": 2.0,
    "angle_extension": (0.0, np.pi / 2),
}
geometry = create_nurbs_patch(
    degree_u=DEGREE, degree_v=DEGREE, nbel_u=nbel, nbel_v=nbel, **geo_args
)
patch_2 = SinglePatch(geometry, quadclass="gs").generate()

boundary_0 = BoundaryCondition(nbctrlpts=patch_0.nbctrlpts, dofs=(DOF.T,))
boundary_0.add_constraint(
    constraint_info=[
        {
            "direction": (ParametricDirection.XI, ParametricDirection.ETA),
            "face": BoundarySide.MIN,
            "dofs": (DOF.T,),
        }
    ],
    constraint_type="dirichlet",
)

boundary_0.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MAX,
            "target_priority": TargetPriority.MASTER,
            "target_id": 1,
        },
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MAX,
            "target_priority": TargetPriority.MASTER,
            "target_id": 2,
        },
    ],
    constraint_type="glue",
)

boundary_1 = BoundaryCondition(nbctrlpts=patch_1.nbctrlpts, dofs=(DOF.T,))
boundary_1.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MIN,
            "dofs": (DOF.T,),
        },
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MAX,
            "dofs": (DOF.T,),
        },
    ],
    constraint_type="dirichlet",
)

boundary_1.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.MIN,
            "target_priority": TargetPriority.SLAVE,
            "target_id": 0,
        },
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MAX,
            "target_priority": TargetPriority.MASTER,
            "target_id": 2,
        },
    ],
    constraint_type="glue",
)

boundary_2 = BoundaryCondition(nbctrlpts=patch_2.nbctrlpts, dofs=(DOF.T,))
boundary_2.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MAX,
            "dofs": (DOF.T,),
        },
        {
            "direction": ParametricDirection.ETA,
            "face": BoundarySide.BOTH,
            "dofs": (DOF.T,),
        },
    ],
    constraint_type="dirichlet",
)

boundary_2.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MIN,
            "target_priority": TargetPriority.SLAVE,
            "target_id": 0,
        },
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MIN,
            "target_priority": TargetPriority.SLAVE,
            "target_id": 1,
        },
    ],
    constraint_type="glue",
)

# Solve thermal problem
model_0 = SpThermalModel(material, patch_0, boundary_0)
model_1 = SpThermalModel(material, patch_1, boundary_1)
model_2 = SpThermalModel(material, patch_2, boundary_2)
multipatchmodel = MpThermalModel({0: model_0, 1: model_1, 2: model_2})

# Add forces
external_force_0 = model_0.assemble_volumetric_force({DOF.ALL: power_density})
external_force_1 = model_1.assemble_volumetric_force({DOF.ALL: power_density})
external_force_2 = model_2.assemble_volumetric_force({DOF.ALL: power_density})
external_force = multipatchmodel.glue(
    {0: external_force_0, 1: external_force_1, 2: external_force_2}
)

# Solve
temperature = np.zeros_like(external_force)
SteadyHeatTransfer(tolerance_nonlinear=1e-7).solve(
    multipatchmodel,
    temperature,
    external_force,
)
temperature_cutted = multipatchmodel.cut(temperature)

relerror_0 = SpaceNormSinglePatch(
    model_0, "l2", {"exact_function": exact_temperature}
).eval(temperature_cutted[0])[-1]

relerror_1 = SpaceNormSinglePatch(
    model_1, "l2", {"exact_function": exact_temperature}
).eval(temperature_cutted[1])[-1]

relerror_2 = SpaceNormSinglePatch(
    model_2, "l2", {"exact_function": exact_temperature}
).eval(temperature_cutted[2])[-1]

print(
    f"Error of multipatch patch: {relerror_0:.2e}, {relerror_1:.2e}, {relerror_2:.2e}"
)

# Postprocess
IgaPostprocessing.plot_patch(
    patch_list=[patch_0, patch_1, patch_2],
    filename="iga_comparison_multipatch",
    folder=RESULT_FOLDER,
    primal_list=[temperature_cutted[0], temperature_cutted[1], temperature_cutted[2]],
    add_colorbar=True,
)

IgaPostprocessing.export_patch(
    patch=patch_0,
    field_type="primal",
    filename=f"iga_comparison_0",
    folder=RESULT_FOLDER,
    fields={"temp": temperature_cutted[0]},
    sample_size=201,
)
IgaPostprocessing.export_patch(
    patch=patch_1,
    field_type="primal",
    filename=f"iga_comparison_1",
    folder=RESULT_FOLDER,
    fields={"temp": temperature_cutted[1]},
    sample_size=201,
)
IgaPostprocessing.export_patch(
    patch=patch_2,
    field_type="primal",
    filename=f"iga_comparison_2",
    folder=RESULT_FOLDER,
    fields={"temp": temperature_cutted[2]},
    sample_size=201,
)
