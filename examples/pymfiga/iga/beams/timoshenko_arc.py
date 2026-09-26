from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.physics import StaticElastoPlasticity
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition, ParametricDirection, BoundarySide, DOF
from yeti_iga.pymfiga.iga.norms import SpaceNormSinglePatch
from yeti_iga.pymfiga.iga.beam_model import TimoshenkoModel
from geomdl import NURBS, operations
import numpy as np

# Global variables
DEGREE, NBEL = 4, 10
YOUNG, POISSON = 1e10, 0.0
RADIUS, THICKNESS, WIDTH = 1.0, 1e-1, 1e-2


def exact_solution(args: dict):
    position = args["position"]
    x = position[0]
    y = position[1]
    r = np.sqrt(x**2 + y**2)
    sin_th = y / r
    cos_th = x / r
    theta = np.arcsin(sin_th)
    a = 12 * RADIUS**2 / (YOUNG * WIDTH)
    u_t = RADIUS / 2 * (sin_th - theta * cos_th)
    u_n = RADIUS / 2 * theta * sin_th
    th_b = sin_th
    return a * np.vstack((u_t, u_n, th_b))


def create_NURBS_quarter_circonference(degree: int, nbel: int) -> NURBS.Curve:

    assert degree > 1, "At least degree 2"

    # Create the vanilla circle using NURBS
    obj = NURBS.Curve()
    obj.degree = 2
    obj.ctrlptsw = [
        [-1, 0, 0, 1],
        [-1 / np.sqrt(2), 1 / np.sqrt(2), 0, 1 / np.sqrt(2)],
        [0, 1, 0, 1],
    ]
    obj.knotvector = [0, 0, 0, 1, 1, 1]

    # Add degree elevation
    if degree > 2:
        operations.degree_operations(obj, [degree - 2])

    # Add knot refinement
    for knot in np.linspace(0.0, 1.0, nbel + 1)[1:-1]:
        operations.insert_knot(obj, [knot], [1])

    return obj


# Define geometry
geometry = create_NURBS_quarter_circonference(DEGREE, NBEL)
patch = SinglePatch(geometry, quadclass="wq", quadtype="2")
patch.reflect(plane="yz")
patch.generate()

# Define plasticity model
material = LinearElasticity(
    {
        "elastic_modulus": YOUNG,
    },
    is_unidimensional=True,
)

# Set boundary conditions
boundary = BoundaryCondition(patch.nbctrlpts, dofs=(DOF.UX, DOF.UY, DOF.ROTZ))
boundary.add_constraint(
    constraint_info=[
        {
            "direction": ParametricDirection.XI,
            "face": BoundarySide.MIN,
            "dofs": (DOF.UX, DOF.UY, DOF.ROTZ),
        }
    ],
    constraint_type="dirichlet",
)

# Define model
model = TimoshenkoModel(axis_patch=patch, material=material, boundary=boundary)
model.add_area_section(WIDTH * THICKNESS, is_uniform=True)
model.add_inertia_section(WIDTH * THICKNESS**3 / 12, is_uniform=True)


# Add external force
def f3(args: dict):
    position = args["position"]
    x = position[0]
    y = position[1]
    r = np.sqrt(x**2 + y**2)
    sin_th = y / r
    return sin_th * THICKNESS**3


external_force = model.assemble_volumetric_force({DOF.ROTZ: f3})

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
np.testing.assert_almost_equal(relerr, 2.816679e-06)
