from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.physics import StaticElastoPlasticity
from yeti_iga.pymfiga.iga.geometry import SinglePatch, GeomdlGenerator
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.norms import SpaceNormSinglePatch
from yeti_iga.pymfiga.iga.beam_model import TimoshenkoModel
import numpy as np

# Global variables
DEGREE, NBEL = 2, 64
YOUNG, POISSON = 1e4, 0.0
LENGTH, THICKNESS, WIDTH = 10, 1, 1
CHARGE = 1.0
S = WIDTH * THICKNESS
I = WIDTH * THICKNESS**3 / 12
G = YOUNG / (2 * (1 + POISSON))


def f2(args: dict):
    position = args["position"]
    x = position[0]
    return CHARGE * np.ones_like(x)


def exact_solution(args: dict):
    position = args["position"]
    x = position[0]
    u = np.zeros_like(x)
    w = x**2 * (6 * LENGTH**2 - 4 * LENGTH * x + x**2) / (24 * YOUNG * I) + x * (
        2 * LENGTH - x
    ) / (2 * G * S)
    t = x * (3 * LENGTH**2 - 3 * LENGTH * x + x**2) / (6 * YOUNG * I)
    return CHARGE * np.vstack((u, w, t))


# Define geometry
geometry = GeomdlGenerator(
    filename="line",
    geo_args={"degree": DEGREE, "nbel": NBEL, "parameters": {"L": LENGTH}},
).export_geometry()
patch = SinglePatch(geometry, quadclass="wq", quadtype="2")
patch.generate()

# Define plasticity model
material = LinearElasticity(
    {
        "elastic_modulus": YOUNG,
    },
    is_unidimensional=True,
)

# Set boundary conditions
boundary = BoundaryCondition(patch.nbctrlpts, dofs=(DOF.UX, DOF.UY, DOF.UZ))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MIN,
            "dofs": (DOF.UX, DOF.UY, DOF.UZ),
        }
    ],
    constraint_type="dirichlet",
)

# Define model
model = TimoshenkoModel(axis_patch=patch, material=material, boundary=boundary)
model.add_area_section(S, is_uniform=True)
model.add_inertia_section(I, is_uniform=True)

# Add external force
external_force = model.assemble_volumetric_force({DOF.UY: f2})

# Solve
displacement = np.zeros_like(external_force)
StaticElastoPlasticity().solve(model, displacement, external_force)

# Adapted from pymfiga_jcf's original `abserr, relerr = ...`: SpaceNormSinglePatch's
# eval() returns a 3rd value (the exact solution's own norm) in this vendored copy,
# needed elsewhere (benchs/pymfiga's spacetime-error tests) -- discarded here.
abserr, relerr, _ = SpaceNormSinglePatch(
    model, "l2", {"exact_function": exact_solution}
).eval(np.reshape(displacement, (3, -1)))
print(f"Relative error matrix-free: {relerr:.3e}")
np.testing.assert_almost_equal(relerr, 3.2787647336482134e-07)
