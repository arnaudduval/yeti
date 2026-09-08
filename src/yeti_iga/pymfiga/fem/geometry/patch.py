from yeti_iga.pymfiga.common.base.cls import BasePatch
from copy import deepcopy
from typing import Tuple
import logging

logger = logging.getLogger("SRC.FEA.GEOMETRY")


class MeshConverter(BasePatch):
    def __init__(self, mesh=None):
        # Geometry refinement level
        self.ndim = 2
        self.mesh_order = "linear"  # or "quadratic"

        # Quadratic mesh data (for P2 Lagrange)
        self.edge_to_midnode = {}
        self.points_quad = []
        self.elements_quad = []
        self.facets_quad = []

        # Edge connectivity (shared by all edge-based spaces)
        self.edge_connectivity = {}

        # Nedelec spaces (can have multiple orders)
        self.nedelec_spaces = {}  # {order: nedelec_data}

        # Initialize data
        self.points_line = []
        self.elements_line = []
        self.facets_line = []

        if mesh is None:
            return

        # Base linear mesh data
        self.points_line = [list(p) for p in mesh.points]
        self.elements_line = [list(e) for e in mesh.elements]
        self.facets_line = [list(f) for f in mesh.facets]
        logger.info(
            f"Initialized linear mesh: {len(self.points_line)} nodes, {len(self.elements_line)} elements"
        )

    def _initialize_quadratic(self):
        self.mesh_order = "quadratic"
        self.points_quad = deepcopy(self.points_line)
        self.elements_quad = deepcopy(self.elements_line)
        self.facets_quad = deepcopy(self.facets_line)
        self.edge_to_midnode = {}

    def _get_midnode(self, i, j):
        edge = tuple(sorted([i, j]))
        if edge in self.edge_to_midnode:
            return self.edge_to_midnode[edge]
        xi, yi = self.points_quad[i]
        xj, yj = self.points_quad[j]
        mid = [(xi + xj) / 2.0, (yi + yj) / 2.0]
        self.points_quad.append(mid)
        mid_index = len(self.points_quad) - 1
        self.edge_to_midnode[edge] = mid_index
        return mid_index

    def convert_linear_to_quadratic(self):
        """Convert linear mesh geometry to quadratic (curved edges)"""
        self._initialize_quadratic()
        quad_elements = []
        for i, j, k in self.elements_quad:
            ij = self._get_midnode(i, j)
            jk = self._get_midnode(j, k)
            ki = self._get_midnode(k, i)
            quad_elements.append([i, j, k, ij, jk, ki])
        self.elements_quad = quad_elements
        quad_facets = []
        for i, j in self.facets_quad:
            m = self._get_midnode(i, j)
            quad_facets.append([i, m, j])
        self.facets_quad = quad_facets
        logger.info(
            f"Converted linear mesh to quadratic: {len(self.points_quad)} nodes, {len(self.elements_quad)} elements"
        )

    def build_edge_connectivity(self):
        """
        Build edge connectivity from linear mesh.

        Returns edge_connectivity dict where:
        - key: edge_id (integer, 0 to num_edges-1)
        - value: dict with:
            - 'nodes': [i, j] corner nodes in canonical order
            - 'elements': list of (elem_id, local_edge_id, orientation)
            - 'is_boundary': bool
        """
        if len(self.edge_connectivity) > 0:
            return self.edge_connectivity

        elements = self.elements_line
        facets = self.facets_line

        # Build boundary edge set
        boundary_edges = set()
        for facet in facets:
            i, j = facet[0], facet[1]
            boundary_edges.add(tuple(sorted([i, j])))

        # Map from sorted edge to edge_id
        edge_to_id = {}
        edge_list = []
        edge_id = 0

        for elem_id, elem in enumerate(elements):
            n0, n1, n2 = elem[0], elem[1], elem[2]

            # Three edges in CCW order: local edge 0, 1, 2
            edges_oriented = [(n0, n1), (n1, n2), (n2, n0)]

            for local_edge_id, edge_oriented in enumerate(edges_oriented):
                i, j = edge_oriented
                edge_key = tuple(sorted([i, j]))

                if (i, j) == edge_key:
                    orientation = +1
                else:
                    orientation = -1

                # Create edge if new
                if edge_key not in edge_to_id:
                    edge_to_id[edge_key] = edge_id
                    edge_list.append(
                        {
                            "nodes": list(edge_key),
                            "elements": [],
                            "is_boundary": edge_key in boundary_edges,
                        }
                    )
                    edge_id += 1

                # Add element reference
                current_edge_id = edge_to_id[edge_key]
                edge_list[current_edge_id]["elements"].append(
                    (elem_id, local_edge_id, orientation)
                )

        self.edge_connectivity = {i: edge_list[i] for i in range(len(edge_list))}
        return self.edge_connectivity

    def initialize_nedelec(self, order=1):
        """
        Initialize Nedelec edge elements of given order.

        For Nedelec elements of order k, we have k DoF points per edge.

        Args:
            order: Order of Nedelec element (1 = lowest order)
        """
        if len(self.edge_connectivity) == 0:
            self.build_edge_connectivity()

        if order in self.nedelec_spaces:
            return  # Already initialized

        nedelec_dofs = []

        # For each edge, create DoF points
        for edge_id, edge_info in self.edge_connectivity.items():
            i, j = edge_info["nodes"]
            xi, yi = self.points_line[i]
            xj, yj = self.points_line[j]

            # Create 'order' equally spaced points along the edge
            dof_points = []
            for k in range(1, order + 1):
                t = k / (order + 1)  # Parameter from 0 to 1
                x = xi + t * (xj - xi)
                y = yi + t * (yj - yi)
                dof_points.append([x, y])

            nedelec_dofs.append(
                {
                    "edge_id": edge_id,
                    "points": dof_points,
                    "global_ids": list(
                        range(
                            len(nedelec_dofs) * order, len(nedelec_dofs) * order + order
                        )
                    ),
                }
            )

        self.nedelec_spaces[order] = nedelec_dofs
        logger.info(
            f"Initialized Nedelec space of order {order} with {len(nedelec_dofs) * order} DoFs"
        )

    def recover_points(self, space_type: str, order: int = 1):
        """
        Recover DoF points based on finite element space.

        Args:
            space_type: 'lagrange' or 'nedelec'
            order: Order of the space (1, 2, ...)

        Returns:
            List of [x, y] coordinates
        """
        if space_type.lower() == "lagrange":
            if order == 1:
                return self.points_line
            elif order == 2:
                if len(self.points_quad) == 0:
                    raise ValueError(
                        "Quadratic mesh not initialized. Call convert_linear_to_quadratic() first."
                    )
                return self.points_quad
            else:
                raise ValueError(f"Lagrange order {order} not supported")
        elif space_type.lower() == "nedelec":
            if order not in self.nedelec_spaces:
                raise ValueError(
                    f"Nedelec space of order {order} not initialized. Call initialize_nedelec({order}) first."
                )
            points = []
            for edge_dof in self.nedelec_spaces[order]:
                points.extend(edge_dof["points"])
            return points
        else:
            raise ValueError(f"Unknown space_type: {space_type}")

    def recover_elements(self, space_type: str, order: int = 1):
        """
        Recover element connectivity (global DoF indices) based on finite element space.

        Args:
            space_type: 'lagrange' or 'nedelec'
            order: Order of the space

        Returns:
            For Lagrange: List of node indices per element
            For Nedelec: List of edge DoF indices per element
                        NOTE: For Nedelec, orientation is NOT handled here!
                        Orientation must be applied during assembly by multiplying
                        by the orientation factor.
        """
        if space_type.lower() == "lagrange":
            if order == 1:
                return self.elements_line
            elif order == 2:
                if len(self.elements_quad) == 0:
                    raise ValueError(
                        "Quadratic mesh not initialized. Call convert_linear_to_quadratic() first."
                    )
                return self.elements_quad
            else:
                raise ValueError(f"Lagrange order {order} not supported")
        elif space_type.lower() == "nedelec":
            if order not in self.nedelec_spaces:
                raise ValueError(
                    f"Nedelec space of order {order} not initialized. Call initialize_nedelec({order}) first."
                )

            nedelec_dofs = self.nedelec_spaces[order]
            num_elements = len(self.elements_line)
            elements_nedelec = [[] for _ in range(num_elements)]

            # For each edge, add its DoF indices to the corresponding elements
            for edge_id, edge_info in self.edge_connectivity.items():
                dof_indices = nedelec_dofs[edge_id]["global_ids"]

                for elem_id, local_edge_id, _ in edge_info["elements"]:
                    # NOTE: orientation is not mandatory here. It is treated when assembling
                    # Store as (local_edge_id, dof_indices) for sorting
                    elements_nedelec[elem_id].append((local_edge_id, dof_indices))

            # Sort by local_edge_id and flatten
            for elem_id in range(num_elements):
                elements_nedelec[elem_id].sort(key=lambda x: x[0])
                elements_nedelec[elem_id] = [
                    dof for _, dofs in elements_nedelec[elem_id] for dof in dofs
                ]

            return elements_nedelec
        else:
            raise ValueError(f"Unknown space_type: {space_type}")

    def recover_facets(self, space_type: str, order: int = 1):
        """Recover boundary facet information"""
        if space_type.lower() == "lagrange":
            if order == 1:
                return self.facets_line
            elif order == 2:
                if len(self.facets_quad) == 0:
                    raise ValueError("Quadratic mesh not initialized.")
                return self.facets_quad
            else:
                raise ValueError(f"Lagrange order {order} not supported")
        elif space_type.lower() == "nedelec":
            if order not in self.nedelec_spaces:
                raise ValueError(f"Nedelec space of order {order} not initialized.")
            nedelec_dofs = self.nedelec_spaces[order]
            boundary_facets = []
            for edge_id, edge_info in self.edge_connectivity.items():
                if edge_info["is_boundary"]:
                    dof_indices = nedelec_dofs[edge_id]["global_ids"]
                    boundary_facets.append(dof_indices)
            return boundary_facets
        else:
            raise ValueError(f"Unknown space_type: {space_type}")

    def get_edge_connectivity(self):
        """Return edge connectivity information"""
        if len(self.edge_connectivity) == 0:
            raise ValueError(
                "Edge connectivity not built. Call build_edge_connectivity() first."
            )
        return self.edge_connectivity

    def get_num_dofs(self, space_type: str, order: int = 1):
        """Return the total number of DoFs for a given space"""
        if space_type.lower() == "lagrange":
            points = self.recover_points(space_type, order)
            return len(points)
        elif space_type.lower() == "nedelec":
            if order not in self.nedelec_spaces:
                raise ValueError(f"Nedelec space of order {order} not initialized.")
            return len(self.nedelec_spaces[order]) * order
        else:
            raise ValueError(f"Unknown space_type: {space_type}")

    def compute_global_mesh_parameter(self) -> Tuple[float, float]:
        raise NotImplementedError()

    def __repr__(self) -> str:
        message = f"""
            GEOMETRY:
            Mesh order: {self.mesh_order}
            Number of nodes (linear): {len(self.points_line)}
            Number of elements (linear): {len(self.elements_line)}
        """
        if self.mesh_order == "quadratic":
            message += f"""
            Number of nodes (quadratic): {len(self.points_quad)}
            Number of elements (quadratic): {len(self.elements_quad)}
            """
        return message
