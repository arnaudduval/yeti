from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import SteadyHeatTransfer
from yeti_iga.pymfiga.common.io import IgaPostprocessing
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.norms import SpaceNormSinglePatch
from yeti_iga.pymfiga.iga.single_model import ThermalModel
from typing import Optional
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)

from matplotlib import pyplot as plt


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


def ders_exact_temperature(args: dict):
    position = args["position"]
    x = position[0]
    y = position[1]
    uders = np.zeros((1, 2, np.size(position, axis=1)))
    uders[0, 0, :] = (
        2 * x * np.sin(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 1)
        + 2 * x * np.sin(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 4)
        + np.pi
        * np.cos(np.pi * x)
        * np.sin(np.pi * y)
        * (x**2 + y**2 - 1)
        * (x**2 + y**2 - 4)
    )
    uders[0, 1, :] = (
        2 * y * np.sin(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 1)
        + 2 * y * np.sin(np.pi * x) * np.sin(np.pi * y) * (x**2 + y**2 - 4)
        + np.pi
        * np.cos(np.pi * y)
        * np.sin(np.pi * x)
        * (x**2 + y**2 - 1)
        * (x**2 + y**2 - 4)
    )
    return uders


def simulate(degree, nbel, quad_args: Optional[dict] = None):
    geo_args = {
        "degree": degree,
        "nbel": nbel,
        "parameters": {"Rin": 1.0, "Rex": 2.0},
    }

    if quad_args is None:
        quad_args = {"quadclass": "gs", "quadtype": "legendre"}

    # Define material
    material = ThermalMaterial()
    material.add_capacity(1, is_uniform=True)
    material.add_conductivity(np.array([[1, 0.5], [0.5, 2]]), is_uniform=True, ndim=2)

    # Define geometry
    geometry = GeomdlGenerator(
        filename="quarter_annulus", geo_args=geo_args
    ).export_geometry()
    patch = SinglePatch(
        geometry, quadclass=quad_args["quadclass"], quadtype=quad_args["quadtype"]
    )
    patch.generate()

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
    model = ThermalModel(material, patch, boundary)
    external_force = model.assemble_volumetric_force({DOF.ALL: power_density})
    temperature = np.zeros_like(external_force)
    SteadyHeatTransfer().solve(model, temperature, external_force)
    return model, temperature


model, temperature = simulate(4, 8)
IgaPostprocessing.export_patch(
    model.part,
    "primal",
    fields={"temp": temperature},
    filename="iga_steady",
    folder=RESULT_FOLDER,
)

degree_list = np.arange(1, 5, dtype=int)
cuts_list = np.arange(1, 7, dtype=int)
relerror_list = np.zeros_like(cuts_list, dtype=float)
color_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]

fig, ax = plt.subplots(figsize=(5.5, 5))
gauss_plot = {"marker": "s", "linestyle": "-", "markersize": 10}
wq1_plot = {"marker": "o", "linestyle": "--", "markersize": 4}
wq2_plot = {"marker": "x", "linestyle": ":", "markersize": 4}

figname = f"{RESULT_FOLDER}/iga_convergence_steady_2d"
for quadclass, quadtype, plotops in zip(
    ["gs", "wq", "wq"], ["legendre", "1", "2"], [gauss_plot, wq1_plot, wq2_plot]
):
    quad_args = {"quadclass": quadclass, "quadtype": quadtype}
    for i, degree in enumerate(degree_list):
        color = color_list[i]
        for j, cuts in enumerate(cuts_list):
            model, temperature = simulate(degree, 2**cuts, quad_args)
            relerror_list[j] = SpaceNormSinglePatch(
                model, "l2", {"exact_function": exact_temperature}
            ).eval(temperature)[-1]

        label = f"IGA-GL deg. {degree}" if quadtype == "legendre" else None
        ax.loglog(
            2**cuts_list,
            relerror_list,
            label=label,
            color=color,
            markerfacecolor="w",
            **plotops,
        )
        fig.savefig(figname)

ax.loglog(
    [],
    [],
    color="k",
    markerfacecolor="w",
    **wq1_plot,
    label="IGA-WQ 1",
)
ax.loglog(
    [],
    [],
    color="k",
    markerfacecolor="w",
    **wq2_plot,
    label="IGA-WQ 2",
)

ax.set_ylim((1e-8, 1))
ax.set_ylabel(r"Relative $||u-u^h||_{L^2(\Omega)}$")
ax.set_xlabel("Number of elements")
ax.legend()
fig.tight_layout()
fig.savefig(figname)
