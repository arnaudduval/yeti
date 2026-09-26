from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.norms import FeaNorm
from yeti_iga.pymfiga.fem.model import MITC_condensed
import numpy as np

# Global constants
YOUNG = 10920.0
POISSON = 0.3
THICKNESS = 1e-1


def w_(pts):
    x, y = pts[0], pts[1]
    a = (1 / 3) * x**3 * y**3 * (x - 1) ** 3 * (y - 1) ** 3
    b1 = y**3 * (y - 1) ** 3 * x * (x - 1) * (5 * x**2 - 5 * x + 1)
    b2 = x**3 * (x - 1) ** 3 * y * (y - 1) * (5 * y**2 - 5 * y + 1)
    b = (2 * THICKNESS**2) / (5 * (1 - POISSON)) * (b1 + b2)
    return a - b


def drx_(pts):
    x, y = pts[0], pts[1]
    return x**3 * (x - 1) ** 3 * y**2 * (y - 1) ** 2 * (2 * y - 1)


def dry_(pts):
    x, y = pts[0], pts[1]
    return -(y**3 * (y - 1) ** 3 * x**2 * (x - 1) ** 2 * (2 * x - 1))


# Material definition
material = LinearElasticity({"elastic_modulus": YOUNG, "poisson_ratio": POISSON})

# Geometry definition
mesh = MeshConverter(
    MeshpyGenerator(filename="square", max_volume=0.001).mesh_geometry()
)
mesh.convert_linear_to_quadratic()

# Boundary conditions
boundary = BoundaryCondition(
    mesh,
    dofs=(DOF.UZ, DOF.ROTX, DOF.ROTY),
    dofs_funspace=[("lagrange", 1), ("lagrange", 2), ("lagrange", 2)],
)
boundary.add_constraint(
    {DOF.UZ: True, DOF.ROTX: True, DOF.ROTY: True}, constraint_type="dirichlet"
)

# Model definition
model = MITC_condensed(
    mesh, THICKNESS, material, FiniteElementQuadrature(boundary_order=2), boundary
)


# Computation of external force
def force(pts):
    def kernel(a, b):
        return (
            12
            * b
            * (b - 1)
            * (5 * a**2 - 5 * a + 1)
            * (2 * b**2 * (b - 1) ** 2 + a * (a - 1) * (5 * b * b - 5 * b + 1))
        )

    cte = YOUNG / (12.0 * (1.0 - POISSON**2)) * THICKNESS**3
    x, y = pts[0], pts[1]
    return cte * (kernel(x, y) + kernel(y, x))


external_force = model.assemble_volumetric_force({DOF.UZ: force})
displacement = model.solve_linearized_system(external_force)

# Compute error
output = model.cut_array(displacement)
nbel = len(model.part.recover_elements("lagrange", 1))
norm = FeaNorm(
    model, "l2", {"degree": "linear", "u_at_nodes": output[DOF.UZ], "fun": w_}
)
abserror, relerror = norm.evaluate()
np.testing.assert_almost_equal(relerror, 0.012585068221728451)
print("Relative error:", relerror, " with ", nbel, " elements")

# Under integrating
# 4      4.380159049586686
# 16     0.42301543740482317
# 150    0.10478326400807882
# 1539   0.009718885160921405
# 15568  0.0009534603534621327
# 155195 0.00010026898865838262

# Good integration
# 4      1.4655731355275181
# 16     0.2315489853975396
# 150    0.0713737796229132
# 1539   0.00790422128389764
# 15568  0.0008777207537882751
