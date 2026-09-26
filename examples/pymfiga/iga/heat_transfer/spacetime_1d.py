from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import SpaceTimeHeatTransfer
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.norms import SpaceTimeNormSinglePatch
from yeti_iga.pymfiga.iga.single_model import SpaceTimeThermalModel
from numpy import sin, cos, tanh, pi, cosh
import numpy as np


def conductivity_property(args: dict):
    temperature = args.get("temperature", np.array([]))
    conductivity = np.zeros(shape=(1, 1, *np.shape(temperature)))
    conductivity[0, 0, ...] = 3.0 + 2.0 * np.tanh(temperature / 50)
    return conductivity


def exact_temperature(args: dict):
    t_list = args["time"]
    position = args["position"]
    x = position[0]
    u = np.zeros((len(t_list), len(x)))
    for i, t in enumerate(t_list):
        u[i] = sin(2 * pi * x) * sin(0.5 * pi * t) * (1.0 + 0.75 * cos(1.5 * pi * t))
    return np.ravel(u)


def power_density(args: dict):
    t_list = args["time"]
    x = args["position"][0]
    f = np.zeros((len(t_list), len(x)))
    for i, t in enumerate(t_list):
        f[i] = (
            -pi
            * (
                0.32
                * pi
                * (0.75 * cos(3 * pi * t / 2) + 1.0) ** 2
                * sin(pi * t / 2) ** 2
                * cos(2 * pi * x) ** 2
                / cosh(
                    (0.015 * cos(3 * pi * t / 2) + 0.02)
                    * sin(pi * t / 2)
                    * sin(2 * pi * x)
                )
                ** 2
                - 8
                * pi
                * (0.75 * cos(3 * pi * t / 2) + 1.0)
                * (
                    2.0
                    * tanh(
                        (0.015 * cos(3 * pi * t / 2) + 0.02)
                        * sin(pi * t / 2)
                        * sin(2 * pi * x)
                    )
                    + 3.0
                )
                * sin(pi * t / 2)
                * sin(2 * pi * x)
                - (0.75 * cos(3 * pi * t / 2) + 1.0) * sin(2 * pi * x) * cos(pi * t / 2)
                + 2.25 * sin(pi * t / 2) * sin(3 * pi * t / 2) * sin(2 * pi * x)
            )
            / 2
        )
    return np.ravel(f)


# Create geometry
DEGREE, NBEL, LENGTH = 8, 128, 1.0
geometry = GeomdlGenerator(
    filename="line",
    geo_args={"degree": DEGREE, "nbel": NBEL, "parameters": {"L": LENGTH}},
).export_geometry()
patch = SinglePatch(geometry, quadclass="wq")
patch.generate()

# Create time span
NBEL_TIME = 16
time_interval = GeomdlGenerator(
    filename="line", geo_args={"degree": 1, "nbel": NBEL_TIME}
).export_geometry()
time_patch = SinglePatch(time_interval, quadclass="gs", quadtype="legendre")
time_patch.generate()

# Create material
material = ThermalMaterial()
material.add_conductivity(conductivity_property, is_uniform=False, ndim=1)
material.add_capacity(1.0, is_uniform=True)

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
model = SpaceTimeThermalModel(material, patch, time_patch, boundary)

# Create external force
external_force = model.assemble_volumetric_force({DOF.ALL: power_density})

# Solve
temperature = np.zeros_like(external_force)
SpaceTimeHeatTransfer().solve(model, temperature, external_force)

# Post processing
error = SpaceTimeNormSinglePatch(
    model, "l2", {"exact_function": exact_temperature}
).eval(temperature)[-1]
np.testing.assert_almost_equal(error, 5.800e-04)
print(f"Relative error is {error:.3e}")
