from .coordinate_transformation import CoordinateTransformation
from typing import List, Tuple, Dict
from copy import deepcopy
import networkx as nx
import numpy as np


class ChainCoordinateTransformation:
    """
    Class to manage a graph of coordinate systems and transformations between them.
    Allows chaining multiple coordinate transformations using a graph structure.
    """

    _verbose = False

    def __init__(self, verbose: bool = False):
        """
        Initialize the transformation graph.

        Args:
            verbose: If True, print detailed information during computations.
        """
        self._graph = nx.DiGraph()  # Directed graph to store bases and transformations
        # Dictionary mapping (base0, base1) to coordinate_transformation objects
        self._coord_transform_list: Dict[Tuple[str, str], CoordinateTransformation] = {}
        self._verbose = verbose

    @property
    def verbose(self):
        return self._verbose

    def _verify_base(self, base: str):
        assert isinstance(base, str)
        assert any(
            base in [transform.base0_name, transform.base1_name, 0, 1]
            for transform in self._coord_transform_list.values()
        )

    def add_transformation(self, coord_transformation: CoordinateTransformation):
        """
        Add a coordinate transformation to the graph.

        Args:
            transform: coordinate_transformation object to add.
        """
        a = coord_transformation.base0_name
        b = coord_transformation.base1_name
        assert isinstance(a, str), "Base 0 name must be set"
        assert isinstance(b, str), "Base 1 name must be set"
        if self.verbose:
            print(f"Adding transformation from {a} to {b} or vice versa")
        coord_transformation.verbose = self.verbose
        # Store transformation for both directions
        a, b = a.lower(), b.lower()
        assert a != b, "Base names must be different"
        self._coord_transform_list.update(
            {(a, b): coord_transformation, (b, a): coord_transformation}
        )
        # Add directed edges for both directions
        self._graph.add_edge(a, b)
        self._graph.add_edge(b, a)

    def _find_path(self, from_base: str, to_base: str) -> List[str]:
        """
        Find the shortest path between two bases in the graph.

        Args:
            from_base: Starting base name.
            to_base: Target base name.

        Returns:
            List of base names representing the path.
        """
        try:
            return nx.shortest_path(self._graph, source=from_base, target=to_base)
        except nx.NetworkXNoPath:
            raise ValueError(
                f"There is not a path from base {from_base} to base {to_base}"
            )

    def transform_point(
        self, point: np.ndarray, from_base: str, to_base: str
    ) -> np.ndarray:
        """
        Transform a point from one base to another, following the shortest path.

        Args:
            point: Coordinates in the starting base.
            from_base: Name of the starting base.
            to_base: Name of the target base.

        Returns:
            Coordinates in the target base.
        """
        from_base = from_base.lower()
        self._verify_base(from_base)
        to_base = to_base.lower()
        self._verify_base(to_base)

        point = np.array(point, dtype=float)
        path = self._find_path(from_base, to_base)
        current_point = np.copy(point)
        if self.verbose:
            print(f"Transforming point {point} from {from_base} to {to_base}")
            print(f"Path: {path}")

        for i in range(len(path) - 1):
            old_base, new_base = path[i], path[i + 1]
            if old_base == new_base:
                continue

            transform: CoordinateTransformation = self._coord_transform_list[
                (old_base, new_base)
            ]
            current_point = transform.transform_point(
                current_point,
                from_base=old_base,
            )

        return current_point

    def transform_tensor(
        self,
        tensor,
        from_base: str,
        to_base: str,
        point: np.ndarray,
        point_base: str,
        comp_type_in: List[str],
        tensor_rank: int,
    ) -> np.ndarray:
        """
        Transform a tensor from one base to another, following the shortest path.

        Args:
            tensor: The tensor to transform.
            from_base: Name of the starting base.
            to_base: Name of the target base.
            point: Coordinates at which to evaluate the transformation.
            point_base: Name of the base in which the point is given.
            comp_type_in: List of 'lower' or 'upper' for each tensor index.
            tensor_rank: Rank of the tensor.

        Returns:
            Transformed tensor as a numpy array.
        """
        from_base = from_base.lower()
        self._verify_base(from_base)
        to_base = to_base.lower()
        self._verify_base(to_base)

        point = np.array(point, dtype=float)
        tensor = np.array(tensor, dtype=float)
        tensor_path = self._find_path(from_base, to_base)
        if self.verbose:
            print(f"Transforming tensor from {from_base} to {to_base}")
            print(f"Path: {tensor_path}")

        current_point = np.copy(point)
        current_tensor = np.copy(tensor)
        current_point_base = deepcopy(point_base)

        for i in range(len(tensor_path) - 1):
            old_base, new_base = tensor_path[i], tensor_path[i + 1]
            if old_base == new_base:
                continue
            if self.verbose:
                print(f"\nTransforming tensor from {old_base} to {new_base}\n")

            transform: CoordinateTransformation = self._coord_transform_list[
                (old_base, new_base)
            ]
            # Transform the point to the new base
            current_point = self.transform_point(
                current_point,
                from_base=str(current_point_base),
                to_base=str(transform.base1_name),
            )
            current_point_base = deepcopy(transform.base1_name)

            # Determine tensor base for transformation
            current_tensor_base = 0 if old_base == transform.base0_name else 1

            # Apply the tensor transformation
            current_tensor = transform.transform_tensor(
                tensor=current_tensor,
                from_base=current_tensor_base,
                point=current_point,
                point_base=1,
                comp_type_in=comp_type_in,
                tensor_rank=tensor_rank,
            )
        return current_tensor
