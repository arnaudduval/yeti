from typing import Callable, List, Tuple, Union, Optional
from scipy.optimize import minimize
import numpy as np
import sympy as sy


class CoordinateTransformation:
    """
    Class for coordinate transformations between two coordinate systems.
    Allows transformation of vectors and tensors between two bases
    (u^j-system and v^k-system), given symbolic transformation functions.
    """

    # Global variables
    _indices = "abcdefghijklmn"

    # Internal variables
    _jac_sympy = None
    _jac_numpy = None
    _func_numpy = None
    _func_inv_numpy = None

    # Properties
    _base0_name = None
    _base1_name = None
    _inverse_bounds = None
    _verbose = False

    def __init__(
        self,
        base0_vars: List[sy.Symbol],
        base1_vars: List[sy.Symbol],
        func: List[sy.Symbol],
        verbose: bool = False,
        base0_name: Optional[str] = None,
        base1_name: Optional[str] = None,
        func_inv: Optional[Callable] = None,
    ):
        """
        Initialize the coordinate transformation.

        Args:
            base0_vars: Variables of the first coordinate system (u^j).
            base1_vars: Variables of the second coordinate system (v^k).
            func: List of transformation functions f^j such that u^j = f^j(v^k).
            verbose: If True, print detailed information during computations.
        """
        assert len(base0_vars) == len(base1_vars), "Base variables must match in length"
        assert len(base0_vars) == len(func), "No. functions != no. variables"
        self._nbofvars = len(base0_vars)
        self._base0_vars = base0_vars  # Variables u^j
        self._base1_vars = base1_vars  # Variables v^k
        self._func = sy.Matrix(func)  # Functions f^j
        self._func_inv = func_inv
        self._verbose = verbose
        if self.verbose:
            print("START: Initializing coordinate transformation:")
            print("Base 0 variables:", self._base0_vars)
            print("Base 1 variables:", self._base1_vars)
            print("Transformation functions:", self._func)
        self._compute_jacobian()
        self._prepare_numeric_functions()
        self.base0_name = base0_name
        self.base1_name = base1_name
        self.inverse_bounds = [(None, None)] * self._nbofvars

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        if isinstance(value, bool):
            self._verbose = value

    @property
    def base0_name(self):
        return self._base0_name

    @base0_name.setter
    def base0_name(self, name):
        if isinstance(name, str):
            self._base0_name = name.lower()

    @property
    def base1_name(self):
        return self._base1_name

    @base1_name.setter
    def base1_name(self, name):
        if isinstance(name, str):
            self._base1_name = name.lower()

    @property
    def inverse_bounds(self):
        """
        Get the bounds for the optimization problem.
        """
        return self._inverse_bounds

    @inverse_bounds.setter
    def inverse_bounds(self, bounds):
        """
        Set the bounds for the optimization problem.
        Args:
            bounds: List of tuples (min, max) for each variable.
        """
        assert (
            len(bounds) == self._nbofvars
        ), "Bounds must match the number of variables"
        assert all(
            isinstance(b, tuple) and len(b) == 2 for b in bounds
        ), "Each bound must be a tuple (min, max)"
        self._inverse_bounds = bounds

    def _verify_base(self, base: Union[str, int]):
        assert isinstance(base, str) or isinstance(base, int)
        assert base in [self.base0_name, self.base1_name, 0, 1]

    def _get_vars_for_derivatives(self):
        """
        Get the variables with respect to which the Jacobian is computed.
        """
        return self._base1_vars

    def _apply_inverse(self, from_base: Union[str, int]) -> bool:
        return True if from_base == self.base1_name or from_base == 1 else False

    def _compute_jacobian(self):
        """
        Compute the symbolic Jacobian matrix of the transformation.
        """
        self._jac_sympy = self._func.jacobian(self._get_vars_for_derivatives())
        if self.verbose:
            print("Jacobian matrix (symbolic):\n", self._jac_sympy)

    def _prepare_numeric_functions(self):
        """
        Prepare a numerical (numpy) version of the Jacobian for fast evaluation.
        """
        self._func_numpy = sy.lambdify(
            self._get_vars_for_derivatives(), self._func, modules="numpy"
        )
        self._jac_numpy = sy.lambdify(
            self._get_vars_for_derivatives(), self._jac_sympy, modules="numpy"
        )

    def _get_transform_matrices(
        self, point: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Evaluate the Jacobian and its pseudoinverse at a given point.

        Args:
            point: Coordinates at which to evaluate.

        Returns:
            (Jacobian, Jacobian inverse) as numpy arrays.
        """
        assert self._jac_numpy is not None
        point = np.array(point, dtype=float)
        matrix = np.array(self._jac_numpy(*point), dtype=float)
        return matrix, np.linalg.inv(matrix)

    def _construct_combinations(
        self, jac: np.ndarray, jac_inv: np.ndarray, apply_inverse: bool
    ):
        """
        Build a dictionary of all possible transformation matrices
        (direct/inverse, lower/upper indices).

        Returns:
            Dictionary with keys like 'lower', 'upper', etc.
        """
        if self.verbose:
            print(f"Computing the Jacobian matrix and its inverse from base 0 to 1.")
        if apply_inverse:
            return {
                "lower": jac_inv.T,  # For lower indices (covariant, inverse)
                "upper": jac,  # For upper indices (contravariant, inverse)
            }
        else:
            return {
                "lower": jac.T,  # For lower indices (covariant)
                "upper": jac_inv,  # For upper indices (contravariant)"
            }

    def _transform_point_from_1_to_0(self, v_coord: np.ndarray) -> np.ndarray:
        """
        Compute the direct transformation numerically (from base 1 to base 0).

        Args:
            v_coord: coordinates in base 1.

        Returns:
            u_coord: exact coordinates in base 0.
        """
        assert self._func_numpy is not None

        v_coord = np.array(v_coord).astype(float)
        assert v_coord.ndim == 1, "Coordinates must be a 1D array"

        u_coord = np.ravel(self._func_numpy(*v_coord))

        if self.verbose:
            print(f"The direct transformation for point {v_coord} is {u_coord}")

        return u_coord

    def _transform_point_from_0_to_1(
        self, u_coord: np.ndarray, initial_guess: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Compute the inverse transformation numerically (from base 0 to base 1).
        Use scipy.optimize.root for solving the nonlinear problem.

        Args:
            u_coord: coordinates in base 0.
            initial_guess: optional initial guess.

        Returns:
            v_coord: approximated coordinates in base 1 such that f(v) = u.
        """
        u_coord = np.array(u_coord).astype(float)

        if callable(self._func_inv):
            # If a custom inverse function is provided, use it directly
            v_coord = self._func_inv(u_coord)

        else:

            assert u_coord.ndim == 1, "Coordinates must be a 1D array"

            def func(v):
                assert self._func_numpy is not None
                return np.linalg.norm(np.ravel(self._func_numpy(*v)) - u_coord) ** 2

            if initial_guess is None:
                initial_guess = np.ones_like(u_coord)

            sol = minimize(func, initial_guess, bounds=self.inverse_bounds)
            if not sol.success:
                raise ValueError(f"Numerical inversion failed: {sol.message}")
            v_coord = sol.x

        if self.verbose:
            print(f"The inverse transformation for point {u_coord} is {v_coord}")

        return v_coord

    def transform_point(
        self,
        point: np.ndarray,
        from_base: Union[str, int],
    ) -> np.ndarray:
        """
        Transform a point from one coordinate system to another.

        Args:
            point: Coordinates in the starting base.
            from_base: Name of the starting base.

        Returns:
            Coordinates in the target base as a numpy array.
        """
        self._verify_base(from_base)
        point = np.array(point, dtype=float)
        if self.verbose:
            print(f"Transforming point {point} from base {from_base} to target base")

        if from_base == self.base0_name or from_base == 0:
            return self._transform_point_from_0_to_1(point)
        elif from_base == self.base1_name or from_base == 1:
            return self._transform_point_from_1_to_0(point)
        else:
            raise ValueError(
                f"Transformation from base {from_base} to target base is not defined"
            )

    def transform_vector(
        self,
        vector: Union[list, tuple, np.ndarray],
        from_base: Union[str, int],
        point: Union[list, np.ndarray],
        point_base: Union[str, int],
        comp_type_in: str,
    ) -> np.ndarray:
        """
        Transform a vector between coordinate systems.

        Args:
            vector: The vector to transform.
            from_base: 0 if the vector is in base 0, 1 if in base 1.
            point: The point at which to evaluate the transformation.
            point_base: 0 if the point is in base 0, 1 if in base 1.
            comp_type_in: 'lower' for covariant or 'upper' for contravariant components.

        Returns:
            Transformed vector as a numpy array.
        """
        self._verify_base(from_base)
        point = np.array(point, dtype=float)
        vector = np.array(vector, dtype=float)
        if point_base != 1 or point_base != self.base1_name:
            point = self.transform_point(point, from_base=point_base)

        # Compute jacobian and its inverse
        jac, jac_inv = self._get_transform_matrices(point)
        combinations = self._construct_combinations(
            jac, jac_inv, self._apply_inverse(from_base)
        )

        comp_type_in = str(comp_type_in).lower()[:5]
        assert comp_type_in in ["lower", "upper"]

        # Select the matrix to compute the transformation
        mat = combinations.get(f"{comp_type_in}")
        assert mat is not None, "Unknown combination"

        if self.verbose:
            print(f"END: Vector was transformed to target base")

        # Perform the transformation using Einstein summation
        return np.einsum("ij,j...->i...", mat, vector, optimize=True)

    def transform_tensor(
        self,
        tensor: np.ndarray,
        from_base: Union[str, int],
        point: np.ndarray,
        point_base: Union[str, int],
        comp_type_in: List[str],
        tensor_rank: int,
    ) -> np.ndarray:
        """
        Transform a tensor of arbitrary rank between coordinate systems.

        Args:
            tensor: The tensor to transform.
            from_base: 0 if the tensor is in base 0, 1 if in base 1.
            point: The point at which to evaluate the transformation.
            point_base: 0 if the point is in base0, 1 if in base1.
            comp_type_in: 'lower' for covariant or 'upper' for contravariant components.
            tensor_rank: Rank of the tensor.

        Returns:
            Transformed tensor as numpy array.
        """
        assert tensor_rank <= len(comp_type_in)
        self._verify_base(from_base)
        point = np.array(point, dtype=float)
        tensor = np.array(tensor, dtype=float)
        if point_base != self.base1_name or point_base != 1:
            point = self.transform_point(point, from_base=point_base)

        # Index labels for einsum
        old_indices = list(self._indices[:tensor_rank].lower())
        new_indices = list(self._indices[:tensor_rank].upper())

        # Compute jacobian and its inverse
        jac, jac_inv = self._get_transform_matrices(point)
        combinations = self._construct_combinations(
            jac, jac_inv, self._apply_inverse(from_base)
        )

        einsum_matrices, einsum_labels = [], []
        for i, type_in in enumerate(comp_type_in[:tensor_rank]):
            type_in = str(type_in).lower()[:5]
            assert type_in in ["lower", "upper"]

            # Select the matrix to compute the transformation
            mat = combinations.get(f"{type_in}")
            assert mat is not None, "Unknown combination"
            einsum_matrices.append(mat)

            # Select the labels for einsum
            old_idx = old_indices[i]
            new_idx = new_indices[i]
            einsum_labels.append(f"{new_idx}{old_idx}")

        # Build einsum string for tensor transformation
        einsum_input = ",".join(einsum_labels) + "," + "".join(old_indices)
        einsum_output = "".join(new_indices)
        einsum_str = f"{einsum_input}...->{einsum_output}..."

        if self.verbose:
            print(f"END: tensor was transformed to target base")

        # Perform the transformation using Einstein summation
        return np.einsum(
            einsum_str,
            *einsum_matrices,
            tensor,
            optimize=True,
        )


class Cartesian2Polar(CoordinateTransformation):
    """
    Coordinate transformation from Cartesian to Polar coordinates and vice-versa.
    Transforms vectors and tensors between these two coordinate systems.
    """

    def __init__(self, verbose: bool = False):
        base0_vars = [sy.Symbol("x"), sy.Symbol("y")]
        base1_vars = [sy.Symbol("r"), sy.Symbol("phi")]
        func = [
            base1_vars[0] * sy.cos(base1_vars[1]),  # x = r * cos(phi)
            base1_vars[0] * sy.sin(base1_vars[1]),  # y = r * sin(phi)
        ]
        func_inv = lambda coords: np.array(
            [
                np.sqrt(coords[0] ** 2 + coords[1] ** 2),
                np.arctan2(coords[1], coords[0]),
            ]
        )
        super().__init__(
            base0_vars,
            base1_vars,
            func,
            verbose=verbose,
            base0_name="cartesian",
            base1_name="polar",
            func_inv=func_inv,
        )


class Cartesian2Cylindrical(CoordinateTransformation):
    """
    Coordinate transformation from Cartesian to Cylindrical coordinates and vice-versa.
    Transforms vectors and tensors between these two coordinate systems.
    """

    def __init__(self, verbose: bool = False):
        base0_vars = [sy.Symbol("x"), sy.Symbol("y"), sy.Symbol("z")]
        base1_vars = [sy.Symbol("r"), sy.Symbol("phi"), sy.Symbol("z")]
        func = [
            base1_vars[0] * sy.cos(base1_vars[1]),  # x = r * cos(phi)
            base1_vars[0] * sy.sin(base1_vars[1]),  # y = r * sin(phi)
            base1_vars[2],  # z = z
        ]
        func_inv = lambda coords: np.array(
            [
                np.sqrt(coords[0] ** 2 + coords[1] ** 2),
                np.arctan2(coords[1], coords[0]),
                coords[2],
            ]
        )
        super().__init__(
            base0_vars,
            base1_vars,
            func,
            verbose=verbose,
            base0_name="cartesian",
            base1_name="cylindrical",
            func_inv=func_inv,
        )


class Cartesian2Spherical(CoordinateTransformation):
    """
    Coordinate transformation from Cartesian to Spherical coordinates and vice-versa.
    Transforms vectors and tensors between these two coordinate systems.
    """

    def __init__(self, verbose: bool = False):
        base0_vars = [sy.Symbol("x"), sy.Symbol("y"), sy.Symbol("z")]
        base1_vars = [sy.Symbol("R"), sy.Symbol("theta"), sy.Symbol("phi")]
        func = [
            base1_vars[0]
            * sy.sin(base1_vars[1])
            * sy.cos(base1_vars[2]),  # x = R * sin(theta) * cos(phi)
            base1_vars[0]
            * sy.sin(base1_vars[1])
            * sy.sin(base1_vars[2]),  # y = R * sin(theta) * sin(phi)
            base1_vars[0] * sy.cos(base1_vars[1]),  # z = R * cos(theta)
        ]
        func_inv = lambda coords: np.array(
            [
                np.sqrt(coords[0] ** 2 + coords[1] ** 2 + coords[2] ** 2),
                np.arccos(
                    coords[2]
                    / np.sqrt(coords[0] ** 2 + coords[1] ** 2 + coords[2] ** 2)
                ),
                np.arctan2(coords[1], coords[0]),
            ]
        )
        super().__init__(
            base0_vars,
            base1_vars,
            func,
            verbose=verbose,
            base0_name="cartesian",
            base1_name="spherical",
            func_inv=func_inv,
        )


class Cylindrical2Spherical(CoordinateTransformation):
    """
    Coordinate transformation from Cylindrical to Spherical coordinates and vice-versa.
    Transforms vectors and tensors between these two coordinate systems.
    """

    def __init__(self, verbose: bool = False):
        base0_vars = [sy.Symbol("r"), sy.Symbol("phi"), sy.Symbol("z")]
        base1_vars = [sy.Symbol("R"), sy.Symbol("Theta"), sy.Symbol("Phi")]
        func = [
            [
                base1_vars[0] * sy.sin(base1_vars[1]),  # r = R * sin(Phi)
                base1_vars[2],  # phi = Phi
                base1_vars[0] * sy.cos(base1_vars[1]),  # z = R * cos(Phi)
            ],
        ]
        func_inv = lambda coords: np.array(
            [
                np.sqrt(coords[0] ** 2 + coords[2] ** 2),
                np.arctan2(coords[0], coords[2]),
                coords[1],
            ]
        )
        super().__init__(
            base0_vars,
            base1_vars,
            func,
            verbose=verbose,
            base0_name="cylindrical",
            base1_name="spherical",
            func_inv=func_inv,
        )
