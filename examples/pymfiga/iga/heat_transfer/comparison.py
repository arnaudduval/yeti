from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import (
    TransientHeatTransfer,
    SpaceTimeHeatTransfer,
)
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.norms import SpaceTimeNormSinglePatch
from yeti_iga.pymfiga.iga.single_model import ThermalModel, SpaceTimeThermalModel
from time import time
from numpy import sin, cos, exp, pi
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


def exact_temperature_sptm(args: dict):
    time = args["time"]
    x = args["position"][0]
    u = np.zeros((len(time), len(x)))
    for i, t in enumerate(time):
        u[i] = sin(pi * x) * t * sin(t) * exp(-t)
    return np.ravel(u)


def power_density_inc(args: dict):
    t = args["time"]
    x = args["position"]
    f = exp(-t) * (-t * sin(t) + sin(t) + t * cos(t)) * sin(pi * x) + pi**2 * sin(
        pi * x
    ) * t * sin(t) * exp(-t)
    return f


def power_density_sptm(args: dict):
    time = args["time"]
    x = args["position"][0]
    f = np.zeros((len(time), len(x)))
    for i, t in enumerate(time):
        f[i] = exp(-t) * (-t * sin(t) + sin(t) + t * cos(t)) * sin(
            pi * x
        ) + pi**2 * sin(pi * x) * t * sin(t) * exp(-t)
    return np.ravel(f)


# Create geometry
degree, nbel, length = 1, 64, 1.0
geometry = GeomdlGenerator(
    filename="line",
    geo_args={"degree": degree, "nbel": nbel, "parameters": {"L": length}},
).export_geometry()
space_patch = SinglePatch(geometry, quadclass="gs")
space_patch.generate()

# Create material
material = ThermalMaterial()
material.add_conductivity(1.0, is_uniform=True, ndim=1)
material.add_capacity(1.0, is_uniform=True)

# Create boundary condition
boundary = BoundaryCondition(nbctrlpts=space_patch.nbctrlpts, dofs=(DOF.T,))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.BOTH,
            "dofs": (DOF.T,),
        },
    ],
    constraint_type="dirichlet",
)

# SPACE TIME METHOD
start = time()
time_interval = GeomdlGenerator(
    filename="line",
    geo_args={
        "degree": 1,
        "nbel": 2 * nbel,
        "parameters": {"L": 1.0},
    },
).export_geometry()
time_patch = SinglePatch(time_interval, quadclass="gs")
time_patch.generate()
model_sptm = SpaceTimeThermalModel(material, space_patch, time_patch, boundary)
external_heat_source_sptm = model_sptm.assemble_volumetric_force(
    {DOF.ALL: power_density_sptm}
)

temperature_sptm = np.zeros_like(external_heat_source_sptm)
SpaceTimeHeatTransfer().solve(
    model_sptm,
    temperature_sptm,
    external_heat_source_sptm,
)
finish = time()
time_elapsed_sptm = finish - start

rel_error_sptm = SpaceTimeNormSinglePatch(
    model_sptm, "l2", {"exact_function": exact_temperature_sptm}
).eval(temperature_sptm)[-1]
np.testing.assert_almost_equal(rel_error_sptm, 0.0001697915499258073)


# INCREMENTAL METHOD
start = time()
time_list = np.linspace(0, 1, 2 * nbel + 1)
model_inc = ThermalModel(material, space_patch, boundary)
external_heat_source_inc = np.zeros((len(time_list), space_patch.nbctrlpts_total))
for i, t in enumerate(time_list):
    external_heat_source_inc[i] = model_inc.assemble_volumetric_force(
        {DOF.ALL: power_density_inc}, time=t
    )

temperature_inc = np.zeros_like(external_heat_source_inc)
TransientHeatTransfer().solve(
    model_inc,
    temperature_inc,
    external_heat_source_inc,
    incremental_type="alpha",
    time_list=time_list,
    alpha=1.0,
)
finish = time()

time_elapsed_inc = finish - start

rel_error_inc = SpaceTimeNormSinglePatch(
    model_sptm, "l2", {"exact_function": exact_temperature_sptm}
).eval(np.ravel(temperature_inc))[-1]
np.testing.assert_almost_equal(rel_error_inc, 0.0014662728550697626)

print(
    f"For space-time method: error {rel_error_sptm:.2e} and time {time_elapsed_sptm:.2f}"
)
print(
    f"For incremental method: error {rel_error_inc:.2e} and time {time_elapsed_inc:.2f}"
)

# Post-processing
from matplotlib import pyplot as plt
from yeti_iga.pymfiga.common.numerics.operations import BsplineOperations

fig, ax = plt.subplots(figsize=(8, 4))
knots_interp = np.linspace(0, 1, 101)
POSITION, TIME = np.meshgrid(knots_interp * length, time_list)
for quarule in space_patch.quadrule_list:
    quarule.knots_to_sample = knots_interp

diff_temp = BsplineOperations.interpolate_meshgrid(
    quadrule_list=space_patch.quadrule_list,
    u_ctrlpts=temperature_inc - np.reshape(temperature_sptm, (len(time_list), -1)),
)

im = ax.contourf(POSITION, TIME, abs(diff_temp), 5, cmap="viridis")
cbar = plt.colorbar(im)
cbar.set_label("Temperature (°C)")

ax.grid(False)
ax.set_ylabel("Time (s)")
ax.set_xlabel("Position (m)")
fig.tight_layout()
fig.savefig(f"{RESULT_FOLDER}/iga_sptm_vs_incr.png")
