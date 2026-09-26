from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import TransientHeatTransfer
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.single_model import ThermalModel
from yeti_iga.pymfiga.common.numerics.operations import BsplineOperations
from matplotlib import pyplot as plt
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


def conductivity_property(args: dict):
    temperature = args.get("temperature", np.array([]))
    conductivity = np.zeros(shape=(1, 1, *np.shape(temperature)))
    conductivity[0, 0, ...] = 3.0  # + 2.0 * np.tanh(temperature / 50)
    return conductivity


def power_density(args: dict):
    t = args["time"]
    x = args["position"]
    f = (
        np.pi * np.cos((np.pi * t) / 2) * np.sin(2 * np.pi * x)
    ) / 2 + 8 * np.pi**2 * np.sin((np.pi * t) / 2) * np.sin(2 * np.pi * x)
    return f


# Create geometry
degree, nbel, length = 4, 16, 1.0
geometry = GeomdlGenerator(
    filename="line",
    geo_args={"degree": degree, "nbel": nbel, "parameters": {"L": length}},
).export_geometry()
patch = SinglePatch(geometry, quadclass="wq")
patch.generate()

# Create material
material = ThermalMaterial()
material.add_conductivity(3.0, is_uniform=True, ndim=1)
material.add_capacity(10.0, is_uniform=True)

# Create boundary condition
boundary = BoundaryCondition(nbctrlpts=patch.nbctrlpts, dofs=(DOF.T,))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.BOTH,
            "dofs": (DOF.T,),
        }
    ],
    constraint_type="dirichlet",
)

# Set transient heat transfer model
model = ThermalModel(material, patch, boundary)

# Create external force
time_list = np.linspace(0, 1, 21)
external_heat_source = np.zeros((len(time_list), patch.nbctrlpts_total))
for i, t in enumerate(time_list):
    external_heat_source[i] = model.assemble_volumetric_force(
        {DOF.ALL: power_density}, time=t
    )

# Solve problem
temperature = np.zeros_like(external_heat_source)
TransientHeatTransfer().solve(
    model,
    temperature,
    external_heat_source,
    incremental_type="alpha",
    time_list=time_list,
    alpha=1.0,
)
np.testing.assert_almost_equal(np.linalg.norm(temperature), 5.8425304268521305)

fig, ax = plt.subplots(figsize=(8, 4))
knots_interp = np.linspace(0, 1, 101)
for quarule in patch.quadrule_list:
    quarule.knots_to_sample = knots_interp
temperature_interp = BsplineOperations.interpolate_meshgrid(
    quadrule_list=patch.quadrule_list,
    u_ctrlpts=temperature,
)
POSITION, TIME = np.meshgrid(knots_interp * length, time_list)
im = ax.contourf(POSITION, TIME, temperature_interp, 20, cmap="viridis")
cbar = plt.colorbar(im)
cbar.set_label("Temperature (°C)")

ax.grid(False)
ax.set_ylabel("Time (s)")
ax.set_xlabel("Position (m)")
fig.tight_layout()
fig.savefig(f"{RESULT_FOLDER}/iga_transient_implicit.png")

#### SOLVE USING BDF
norder = 1
# Solve problem
temperature1 = np.zeros_like(external_heat_source)
TransientHeatTransfer().solve(
    model,
    temperature1,
    external_heat_source,
    incremental_type="bdf",
    tspan=(time_list[0], time_list[-1]),
    nsteps=len(time_list) - 1,
    norder=norder,
)

fig, ax = plt.subplots(figsize=(8, 4))
temperature_interp = BsplineOperations.interpolate_meshgrid(
    quadrule_list=patch.quadrule_list,
    u_ctrlpts=temperature1,
)
POSITION, TIME = np.meshgrid(knots_interp * length, time_list)
im = ax.contourf(POSITION, TIME, temperature_interp, 20, cmap="viridis")
cbar = plt.colorbar(im)
cbar.set_label("Temperature (°C)")

ax.grid(False)
ax.set_ylabel("Time (s)")
ax.set_xlabel("Position (m)")
fig.tight_layout()
fig.savefig(f"{RESULT_FOLDER}/iga_transient_bdf_{norder}.png")

err = np.linalg.norm(temperature - temperature1)
print(f"Difference between BDF and direct method: {err}")
np.testing.assert_almost_equal(err, 0.0)
