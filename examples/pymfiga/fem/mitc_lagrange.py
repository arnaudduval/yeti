from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.fem.geometry import MeshpyGenerator, MeshConverter
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition, DOF
from yeti_iga.pymfiga.fem.norms import FeaNorm
from yeti_iga.pymfiga.fem.model import MITC_condensed, MITC_lagrange
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
import numpy as np

# Global constants
YOUNG = 10920
POISSON = 0.3
THICKNESS = 1e-1


def cut_matrix(output_matrix):
    LL = output_matrix[(DOF.UZ, DOF.UZ)]
    MM1 = output_matrix[(DOF.ROTX, DOF.ROTX)]
    MM2 = output_matrix[(DOF.ROTX, DOF.ROTY)]
    MM3 = output_matrix[(DOF.ROTY, DOF.ROTX)]
    MM4 = output_matrix[(DOF.ROTY, DOF.ROTY)]
    MM = np.block([[MM1, MM2], [MM3, MM4]])
    AA = np.block(
        [
            [LL, np.zeros((LL.shape[0], MM.shape[1]))],
            [np.zeros((MM.shape[0], LL.shape[1])), MM],
        ]
    )
    BB = output_matrix[(DOF.LAG1, DOF.LAG1)]
    NN1 = output_matrix[(DOF.UZ, DOF.LAG2)]
    NN2 = output_matrix[(DOF.ROTX, DOF.LAG2)]
    NN3 = output_matrix[(DOF.ROTY, DOF.LAG2)]
    PP = np.block([[NN2], [NN3]])
    CC = np.block([[NN1], [NN2], [NN3]])
    DD = output_matrix[((DOF.LAG1, DOF.LAG2))]
    return AA, BB, CC, DD


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
    MeshpyGenerator(filename="square", max_volume=0.05).mesh_geometry()
)
mesh.convert_linear_to_quadratic()
mesh.initialize_nedelec(order=1)
mesh.build_edge_connectivity()

############# CONDENSATION #############
# Boundary conditions with condensation
boundary_cond = BoundaryCondition(
    mesh,
    dofs=(DOF.UZ, DOF.ROTX, DOF.ROTY),
    dofs_funspace=[("lagrange", 1), ("lagrange", 2), ("lagrange", 2)],
)
boundary_cond.add_constraint(
    {DOF.UZ: True, DOF.ROTX: True, DOF.ROTY: True}, constraint_type="dirichlet"
)

# Model definition with condensation
model_cond = MITC_condensed(
    mesh, THICKNESS, material, FiniteElementQuadrature(boundary_order=2), boundary_cond
)
mat_cond_1 = model_cond.assemble_stiffness().toarray()

free_dofs, constraint_dofs = boundary_cond.select_nodes_for_solving()
assert isinstance(free_dofs, list)
mat = mat_cond_1[np.ix_(free_dofs, free_dofs)]
print(f"Cond of matrix (condensed): {np.linalg.cond(mat):.3e}")

############# LAGRANGE #############
# Boundary conditions with lagrange
boundary_lag = BoundaryCondition(
    mesh,
    dofs=(DOF.UZ, DOF.ROTX, DOF.ROTY, DOF.LAG1, DOF.LAG2),
    dofs_funspace=[
        ("lagrange", 1),
        ("lagrange", 2),
        ("lagrange", 2),
        ("nedelec", 1),
        ("nedelec", 1),
    ],
)
boundary_lag.add_constraint(
    {DOF.UZ: True, DOF.ROTX: True, DOF.ROTY: True},
    constraint_type="dirichlet",
)

# Model definition with lagrange
model_lag = MITC_lagrange(
    mesh, THICKNESS, material, FiniteElementQuadrature(boundary_order=2), boundary_lag
)
matrix_lag = model_lag.assemble_stiffness().toarray()
output_matrix = model_lag.cut_matrix(matrix_lag)
AA, BB, CC, DD = cut_matrix(output_matrix)

# Compute condensed matrix
CC_invDD = CC @ np.linalg.inv(DD)
mat_condensed_2 = AA + CC_invDD @ BB @ CC_invDD.T
diff = mat_condensed_2 - mat_cond_1
np.testing.assert_almost_equal(np.linalg.norm(diff), 0.0)


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


external_force = model_lag.assemble_volumetric_force({DOF.UZ: force})

# Solve with lagrange (no condensation)
free_dofs, constraint_dofs = boundary_lag.select_nodes_for_solving()
assert isinstance(free_dofs, list)
mat = matrix_lag[np.ix_(free_dofs, free_dofs)]
print(f"Cond of matrix (lagrange): {np.linalg.cond(mat):.3e}")

displacement = np.zeros_like(external_force)
displacement[free_dofs] = LinearSolver.direct(
    matrix_lag[np.ix_(free_dofs, free_dofs)], external_force[free_dofs]
)["sol"]

# Postprocessing
output = model_lag.cut_array(displacement)
nbel = len(model_lag.part.recover_elements("lagrange", 1))
norm = FeaNorm(
    model_lag, "l2", {"degree": "linear", "u_at_nodes": output[DOF.UZ], "fun": w_}
)
abserror, relerror = norm.evaluate()
np.testing.assert_almost_equal(relerror, 0.4677923062000459)
print("Relative error:", relerror, " with ", nbel, " elements")
