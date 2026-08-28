from yeti_iga.pymfiga.common.numerics.quadrature_rules import StandardGauss
from geomdl import NURBS, operations
from scipy.sparse import linalg as scsplin
from scipy import linalg as sclin
from typing import Tuple, List
import numpy as np


def make_uniform_knotvector(
    degree: int, nbel: int, multiplicity: int = 1
) -> np.ndarray:
    knotvector = np.concatenate(
        (
            np.zeros(degree + 1),
            np.repeat(np.linspace(0.0, 1.0, nbel + 1)[1:-1], multiplicity),
            np.ones(degree + 1),
        )
    )
    return knotvector


def make_periodic_knotvector(ukv, degree):
    assert len(ukv) >= 2, "Knot vector must contain at least 2 knots"
    ukv = np.asarray(ukv)
    n = len(ukv)
    m = degree // (n - 1) + 1
    kv_right = ukv[1:] + 1.0
    kv_left = ukv[-2::-1] - 1.0
    for i in range(1, m):
        kv_right = np.concatenate((kv_right, kv_right[: n - 1] + i))
        kv_left = np.concatenate((kv_left, kv_left[: n - 1] - i))
    return np.concatenate((kv_left[:degree][::-1], ukv, kv_right[:degree]))


def make_BSPLINE_line(degree: int, nbel: int) -> Tuple[np.ndarray, np.ndarray]:
    knotvector = make_uniform_knotvector(degree, nbel)
    nbctrlpts = len(knotvector) - degree - 1
    ctrlpts = np.array(
        [
            sum(knotvector[i + j + 1] for j in range(degree)) / degree
            for i in range(0, nbctrlpts)
        ]
    )
    return knotvector, ctrlpts


def make_BSPLINE_quarter_circonference(
    degree: int, nbel: int
) -> Tuple[np.ndarray, np.ndarray]:
    # Create a uniform knot vector for the given degree and number of elements
    knotvector = make_uniform_knotvector(degree, nbel)

    # Perform Gaussian quadrature on the knot vector
    quadrature = StandardGauss(degree, knotvector, quadtype="legendre")
    quadrature.export_quadrature_rules()

    # Compute the matrix using the weights and basis functions from the quadrature
    matrix = quadrature.weights[0] @ quadrature.basis[0]

    # Compute the vector using the weights and quadrature points, applying a cosine function
    rhs = quadrature.weights[0] @ np.cos(np.pi / 2 * quadrature.quadpts)

    # Extract submatrices and vectors
    Ann = matrix[1:-1, 1:-1]
    And = matrix[1:-1, [0, -1]]
    bn = rhs[1:-1]

    # Initialize solution vector
    u = np.zeros_like(rhs)
    u[0] = 1.0  # Boundary condition
    ud = u[[0, -1]]  # Dirichlet boundary values

    # Solve the linear system for the interior points
    u[1:-1] = scsplin.spsolve(Ann, bn - And @ ud)

    # Construct control points
    ctrlpts = np.vstack((u, np.flip(u))).T

    return knotvector, ctrlpts


def make_BSPLINE_periodic_circonference(
    degree: int, nbel: int, **geo_args
) -> Tuple[np.ndarray, np.ndarray]:
    assert degree >= 2, "At least quadratic spline"
    assert nbel > degree, "At least 3 circular sectors"
    radius = geo_args.get("radius", 1.0)
    center = geo_args.get("center", (0.0, 0.0))

    def solve_constraint_system(A: np.ndarray, b: np.ndarray, p: int):
        assert A.ndim == 2, "A must be a matrix"
        n, m = np.shape(A)
        assert n == m, "A must be square"
        P = np.identity(n - p)
        P = np.vstack((P, P[:p, :]))
        xp = sclin.lstsq(A @ P, b)
        assert xp is not None, "Linear system could not be solved"
        x = P @ xp[0]
        return x

    # Create a uniform knot vector for the given degree and number of elements
    ukv = np.linspace(0.0, 1.0, nbel + 1)
    knotvector = make_periodic_knotvector(ukv, degree)

    # Perform Gaussian quadrature on the knot vector
    quadrature = StandardGauss(
        degree,
        knotvector,
        quadtype="legendre",
    )
    quadrature.export_quadrature_rules()

    # Compute the matrix using the weights and basis functions from the quadrature
    matrix = quadrature.weights[0] @ quadrature.basis[0]

    # Compute the vector using the weights and quadrature points, applying a cosine function
    rhs = quadrature.weights[0] @ np.cos(2 * np.pi * quadrature.quadpts)

    # Solve the linear system for the interior points
    u = solve_constraint_system(matrix, rhs, degree)

    # Compute the vector using the weights and quadrature points, applying a cosine function
    rhs = quadrature.weights[0] @ np.sin(2 * np.pi * quadrature.quadpts)

    # Solve the linear system for the interior points
    v = solve_constraint_system(matrix, rhs, degree)

    # Construct control points
    ctrlpts = np.vstack((center[0] + radius * u, center[1] + radius * v)).T

    return knotvector, ctrlpts[:-degree]


def create_NURBS_line(
    degree: int, nbel: int, p_ini: List[float], p_end: List[float]
) -> NURBS.Curve:

    obj = NURBS.Curve()
    obj.degree = 1
    obj.ctrlptsw = [[p[0], p[1], 0.0, 1.0] for p in [p_ini, p_end]]
    obj.knotvector = [0, 0, 1, 1]

    # Add degree elevation
    if degree > 1:
        operations.degree_operations(obj, [degree - 1])

    # Add knot refinement
    for knot in np.linspace(0.0, 1.0, nbel + 1)[1:-1]:
        operations.insert_knot(obj, [knot], [1])
    return obj


def create_NURBS_arc(
    degree: int, nbel: int, alpha_ini: float, alpha_end: float
) -> NURBS.Curve:

    assert degree > 1
    # Arc parameters
    beta = 0.5 * (alpha_ini + alpha_end)
    w = np.cos((alpha_end - alpha_ini) / 2)

    P0 = [np.cos(alpha_ini), np.sin(alpha_ini), 0.0, 1.0]
    P1 = [np.cos(beta), np.sin(beta), 0.0, w]
    P2 = [np.cos(alpha_end), np.sin(alpha_end), 0, 1.0]

    # Create the vanilla circle using NURBS
    obj = NURBS.Curve()
    obj.degree = 2
    obj.ctrlptsw = [P0, P1, P2]
    obj.knotvector = [0, 0, 0, 1, 1, 1]

    # Add degree elevation
    if degree > 2:
        operations.degree_operations(obj, [degree - 2])

    # Add knot refinement
    for knot in np.linspace(0.0, 1.0, nbel + 1)[1:-1]:
        operations.insert_knot(obj, [knot], [1])

    return obj
