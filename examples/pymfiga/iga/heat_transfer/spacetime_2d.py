from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.physics import SpaceTimeHeatTransfer
from yeti_iga.pymfiga.iga.geometry import GeomdlGenerator, SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.norms import SpaceTimeNormSinglePatch
from yeti_iga.pymfiga.iga.single_model import SpaceTimeThermalModel
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from numpy import sin, cos, tanh, pi, cosh
import numpy as np
import time
import os

RESULT_FOLDER = os.path.join(os.getcwd(), "results/")
if not os.path.isdir(RESULT_FOLDER):
    os.mkdir(RESULT_FOLDER)


# Global constants
CST = 50
DEGREE, NBEL = 3, 16


def exact_temperature(args: dict):
    t_list = args["time"]
    position = args["position"]
    x = position[0]
    y = position[1]
    u = np.zeros((len(t_list), len(x)))
    for i, t in enumerate(t_list):
        u[i] = (
            CST
            * tanh(pi * (1.0 - x**2 - y**2))
            * sin(pi * (x**2 + y**2 - 0.25**2))
            * sin(2 * pi * x * y)
            * sin(pi * t / 2)
        )
    return np.ravel(u)


def power_density(args: dict):
    position = args["position"]
    t_list = args["time"]
    x = position[0]
    y = position[1]
    f = np.zeros((len(t_list), len(x)))
    for i, t in enumerate(t_list):
        f[i] = (
            pi
            * CST
            * (
                -24.0
                * pi
                * x**2
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * sin(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                - 16.0
                * pi
                * x**2
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * sin(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 8.0
                * pi
                * x**2
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * cos(2 * pi * x * y)
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 16.0
                * pi
                * x**2
                * sin(pi * t / 2)
                * sin(2 * pi * x * y)
                * cos(pi * (x**2 + y**2 - 0.0625))
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 8.0
                * pi
                * x**2
                * sin(pi * t / 2)
                * cos(pi * (x**2 + y**2 - 0.0625))
                * cos(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                - 16.0
                * pi
                * x
                * y
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * sin(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                - 16.0
                * pi
                * x
                * y
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * sin(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 48.0
                * pi
                * x
                * y
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * cos(2 * pi * x * y)
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 16.0
                * pi
                * x
                * y
                * sin(pi * t / 2)
                * sin(2 * pi * x * y)
                * cos(pi * (x**2 + y**2 - 0.0625))
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 48.0
                * pi
                * x
                * y
                * sin(pi * t / 2)
                * cos(pi * (x**2 + y**2 - 0.0625))
                * cos(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                - 24.0
                * pi
                * y**2
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * sin(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                - 32.0
                * pi
                * y**2
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * sin(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 8.0
                * pi
                * y**2
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * cos(2 * pi * x * y)
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 32.0
                * pi
                * y**2
                * sin(pi * t / 2)
                * sin(2 * pi * x * y)
                * cos(pi * (x**2 + y**2 - 0.0625))
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 8.0
                * pi
                * y**2
                * sin(pi * t / 2)
                * cos(pi * (x**2 + y**2 - 0.0625))
                * cos(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                + 12.0
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * sin(2 * pi * x * y)
                / cosh(pi * (x**2 + y**2 - 1.0)) ** 2
                + 4.0
                * sin(pi * t / 2)
                * sin(pi * (x**2 + y**2 - 0.0625))
                * cos(2 * pi * x * y)
                * tanh(pi * (x**2 + y**2 - 1.0))
                + 12.0
                * sin(pi * t / 2)
                * sin(2 * pi * x * y)
                * cos(pi * (x**2 + y**2 - 0.0625))
                * tanh(pi * (x**2 + y**2 - 1.0))
                - sin(pi * (x**2 + y**2 - 0.0625))
                * sin(2 * pi * x * y)
                * cos(pi * t / 2)
                * tanh(pi * (x**2 + y**2 - 1.0))
            )
            / 2
        )
    return np.ravel(f)


# Create model
geometry = GeomdlGenerator(
    filename="quarter_annulus",
    geo_args={"name": "QA", "degree": (DEGREE, DEGREE + 1), "nbel": NBEL},
).export_geometry()
space_patch = SinglePatch(geometry, quadclass="gs", quadtype="legendre")
space_patch.generate()

# Create time span
NBEL_TIME = 16
time_interval = GeomdlGenerator(
    filename="line", geo_args={"degree": 2, "nbel": NBEL_TIME}
).export_geometry()
time_patch = SinglePatch(time_interval, quadclass="gs", quadtype="legendre")
time_patch.generate()

# Add material
material = ThermalMaterial()
material.add_capacity(1, is_uniform=True)
material.add_conductivity(np.array([[1.0, 0.5], [0.5, 2.0]]), is_uniform=True, ndim=2)

# Block boundaries
boundary = BoundaryCondition(nbctrlpts=space_patch.nbctrlpts, dofs=(DOF.T,))
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

# Define space time model
model = SpaceTimeThermalModel(material, space_patch, time_patch, boundary)
contraint_nodes = model.get_free_and_constraint_nodes()[-1]

# Add external heat force
external_force = model.assemble_volumetric_force({DOF.ALL: power_density})
external_force[contraint_nodes] = 0.0

# Solve space time problem as linear
linear_solver = LinearSolver(tolerance=1e-8, maxiters=100, cleandod=contraint_nodes)

output = model.solve_linearized_system(
    external_force,
    linear_solver_backend=linear_solver,
    spacetime_type="picard",
)

# Solve space time problem as nonlinear
start = time.time()
temperature = np.zeros_like(external_force)
SpaceTimeHeatTransfer().solve(model, temperature, external_force)
finish = time.time()

error = np.linalg.norm(output - temperature) / np.linalg.norm(temperature)
np.testing.assert_array_less(error, 1.0e-8)

# Post processing
error = SpaceTimeNormSinglePatch(
    model, "l2", {"exact_function": exact_temperature}
).eval(temperature)[-1]
np.testing.assert_array_less(error, 1e-4)
print(f"Relative error is {error:.3e} in {finish-start:.2e}s")
