from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import SteadyHeatTransfer
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.norms import SpaceNormSinglePatch
from yeti_iga.pymfiga.iga.single_model import ThermalModel
from matplotlib import pyplot as plt
import numpy as np
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


def power_density(args: dict):
    x = args["position"]
    f = 4 * np.pi**2 * np.sin(2 * np.pi * x)
    return f


def exact_temperature(args: dict):
    x = args["position"]
    u = np.sin(2 * np.pi * x)
    return u


def simulate(degree, nbel, quadclass, quadtype):
    # Create geometry
    geometry = GeomdlGenerator(
        filename="line",
        geo_args={"degree": degree, "nbel": nbel, "parameters": {"L": 1.0}},
    ).export_geometry()
    patch = SinglePatch(geometry, quadclass=quadclass, quadtype=quadtype)
    patch.generate()

    # Create material
    material = ThermalMaterial()
    material.add_capacity(1.0, is_uniform=True)
    material.add_conductivity(1.0, is_uniform=True, ndim=1)

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
    external_force = model.assemble_volumetric_force({DOF.ALL: power_density})

    # Solve
    temperature = np.zeros_like(external_force)
    SteadyHeatTransfer().solve(model, temperature, external_force)
    return model, temperature


degree_list = np.arange(1, 5)
cuts_list = np.arange(1, 8)
relerror_list = np.zeros_like(cuts_list, dtype=float)
color_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]

fig, ax = plt.subplots(figsize=(5.5, 5))
gauss_plot = {"marker": "s", "linestyle": "-", "markersize": 10}
wq1_plot = {"marker": "o", "linestyle": "--", "markersize": 4}
wq2_plot = {"marker": "x", "linestyle": ":", "markersize": 4}

figname = f"{RESULT_FOLDER}/iga_convergence_steady_1d.pdf"
for quadclass, quadtype, plotops in zip(
    ["gs", "wq", "wq"],
    ["legendre", "1", "2"],
    [gauss_plot, wq1_plot, wq2_plot],
):
    for i, degree in enumerate(degree_list):
        color = color_list[i]
        for j, cuts in enumerate(cuts_list):
            model, temperature = simulate(degree, 2**cuts, quadclass, quadtype)
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

ax.set_ylim((1e-10, 1))
ax.set_ylabel(r"Relative $||u-u^h||_{L^2(\Omega)}$")
ax.set_xlabel("Number of elements")
ax.legend()
fig.tight_layout()
fig.savefig(figname)
