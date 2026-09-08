from yeti_iga.pymfiga.common.base.enum import DOF
from yeti_iga.pymfiga.common.base.cls import (
    BaseBoundaryCondition,
    BoundarySide,
    ParametricDirection,
    TargetPriority,
)
from typing import Callable, Dict, List, Set, Tuple, Union, Sequence, Literal
import numpy as np
import logging

logger = logging.getLogger("SRC.FEA.BOUNDARY")


def rearrange_points(
    indices: Union[List, np.ndarray, Tuple], points: Union[List, np.ndarray, Tuple]
) -> List:
    def key(p):
        x, y = points[p]
        angle = np.atan2(y, x)
        radius = x * x + y * y
        return (angle, radius)

    return sorted(indices, key=lambda p: key(p))


class BoundaryCondition(BaseBoundaryCondition):
    """
    BoundaryCondition manages boundary constraints (Dirichlet, Neumann, Contact)
    for a 2D mesh generated with meshpy.
    """

    def __init__(
        self, mesh, dofs: Sequence[DOF], dofs_funspace: Sequence[Tuple[str, int]] = []
    ):
        logger.info("Initialize boundary conditions")
        self._verify_dofs(dofs)
        dofs_funspace = self._set_dofs_funspace(dofs_funspace)
        self._dofs_funspace = dofs_funspace[: len(dofs)]
        self._dofs_index: Dict[DOF, Tuple[int, Tuple[str, int]]] = {
            dof: (i, funspace)
            for i, (dof, funspace) in enumerate(zip(dofs, self._dofs_funspace))
        }
        self._mesh = self._verify_mesh(mesh)

        # Boundary condition storage
        self._nodes_dirichlet: Dict[DOF, Set[int]] = {}
        self._nodes_contact: Set[int] = set()
        self._segments_contact: Set[tuple] = set()
        self._nbofdofs: Dict[DOF, int] = {}
        self._limits_dofs: Dict[DOF, Tuple[int, int]] = {}
        logger.info(repr(self))

    @property
    def dofs_index(self):
        return self._dofs_index

    @property
    def nbofdofs(self):
        return self._nbofdofs

    @property
    def limits_dofs(self):
        return self._limits_dofs

    def _verify_mesh(self, mesh):
        assert hasattr(mesh, "recover_points") and hasattr(
            mesh, "recover_facets"
        ), AttributeError()
        return mesh

    def _verify_dofs(self, dofs: Sequence[DOF]):
        assert len(dofs) > 0 and all(
            dof != DOF.ALL for dof in dofs
        ), "DOFs list must be non-empty and not contain DOF.ALL."

    def _set_dofs_funspace(self, dof_funspace: Sequence[Tuple[str, int]]):
        assert isinstance(dof_funspace, (list, tuple))
        if len(dof_funspace) == 0:
            dof_funspace = [("lagrange", 1) for _ in range(len(DOF))]
        for _ in dof_funspace:
            assert _[0] in ["lagrange", "nedelec"], "Unsupported element type."
            assert _[1] in [1, 2], "Unsupported element order."
        return dof_funspace

    def _recognize_boundary(
        self, funspace: Tuple[str, int]
    ) -> Tuple[List[tuple], Set[int]]:
        """Identify all boundary segments and nodes from the mesh."""
        idx_boundary_segments = []
        idx_boundary_nodes = set()
        facets = self._mesh.recover_facets(funspace[0], funspace[1])
        points = self._mesh.recover_points(funspace[0], funspace[1])
        for facet in facets:
            rearranged = rearrange_points(facet, points)
            idx_boundary_segments.append(tuple(rearranged))
            idx_boundary_nodes.update(rearranged)

        return idx_boundary_segments, idx_boundary_nodes

    def recognize_constraint(
        self,
        cnstr: Union[bool, Callable],
        funspace: Tuple[str, int],
    ) -> Tuple[Set[int], Set[tuple]]:
        """
        Recognize nodes and segments satisfying a condition.
        Condition can be:
        - bool: select all (True) or none (False)
        - callable: takes a point -> bool
        """
        idx_boundary_segments, idx_boundary_nodes = self._recognize_boundary(funspace)
        nb_boundary_segments = len(idx_boundary_segments)
        if isinstance(cnstr, bool):
            nodes = idx_boundary_nodes if cnstr else set()
            idx_segments = [i for i in range(nb_boundary_segments)] if cnstr else []

        elif callable(cnstr):
            nodes: Set[int] = set()
            idx_segments: List[int] = []
            points = self._mesh.recover_points(funspace[0], funspace[1])
            for idx_seg in range(nb_boundary_segments):
                idx_pair = idx_boundary_segments[idx_seg]
                flags = np.array([cnstr(points[idx]) for idx in idx_pair], dtype=bool)
                if np.all(flags):
                    idx_segments.append(idx_seg)
                nodes.update({n for n, flag in zip(idx_pair, flags) if flag})

        else:
            raise TypeError("Condition must be a bool or a callable.")

        idx_nodes_segment_list = set([idx_boundary_segments[_] for _ in idx_segments])

        return nodes, idx_nodes_segment_list

    def add_constraint(
        self,
        constraint_info: Union[List, Dict],
        constraint_type: Literal["dirichlet", "contact"],
    ):
        """Add constraints to the boundary condition set."""
        t = str(constraint_type).lower()

        if t == "contact":
            assert isinstance(constraint_info, list)
            for _ in self._dofs_funspace:
                assert _ == (
                    "lagrange",
                    1,
                ), "Contact conditions only supported for linear Lagrange elements."
            for cnstr in constraint_info:
                assert callable(cnstr), "Contact conditions must be callable."
                nodes, idx_nodes_segment_list = self.recognize_constraint(
                    cnstr, ("lagrange", 1)
                )
                self._nodes_contact.update(nodes)
                # TODO: verify if working with many constraints
                self._segments_contact.update(idx_nodes_segment_list)

        elif t == "dirichlet":
            assert isinstance(
                constraint_info, dict
            ), "Dirichlet conditions must be provided as a dict 'DOF: condition'."
            for dof, cnstr in constraint_info.items():
                assert isinstance(dof, DOF), "Dirichlet keys must be DOFs."
                funspace = self._dofs_index[dof][1]
                nodes = self.recognize_constraint(cnstr, funspace)[0]
                if dof in self._nodes_dirichlet:
                    self._nodes_dirichlet[dof].update(nodes)
                else:
                    self._nodes_dirichlet[dof] = set(nodes)

        else:
            raise ValueError(f"Unknown constraint type: {constraint_type}")

    def select_nodes_for_solving(self) -> Tuple[List[int], List[int]]:
        """
        Select nodes not subject to Dirichlet conditions.

        Returns:
            Tuple[list, list]: Lists of free and constrained node indices.
        """
        if not self._nodes_dirichlet:
            logger.warning("Warning: No Dirichlet conditions set.")

        def ravel(d: dict) -> List[int]:
            return [n for nodes in d.values() for n in nodes]

        # Initialize dictionaries for free and constraint nodes
        free_nodes: Dict[DOF, List[int]] = {}
        constraint_nodes: Dict[DOF, List[int]] = {}
        self._nbofdofs = {}
        self._limits_dofs = {}
        offset_old, offset_new = 0, 0

        for dof, (_, funspace) in self._dofs_index.items():
            # Update DOF counts and limits
            nbnodes = self._mesh.get_num_dofs(funspace[0], funspace[1])
            offset_new += nbnodes
            self._nbofdofs.update({dof: nbnodes})
            self._limits_dofs.update({dof: (offset_old, offset_new)})

            # Compute free and constraint sets
            all_nodes = set(range(nbnodes))
            constraint_set = self._nodes_dirichlet.get(dof, set())
            free_set = sorted(all_nodes.difference(constraint_set))

            # Define dictionaries of free and constraint nodes
            free_nodes[dof] = [n + offset_old for n in free_set]
            constraint_nodes[dof] = [n + offset_old for n in list(constraint_set)]

            # Update offset for next DOF
            offset_old = offset_new

        return ravel(free_nodes), ravel(constraint_nodes)

    def select_nodes_for_contact(self) -> Tuple[Set[int], Set[tuple]]:
        """Return the sets of contact nodes and segments."""
        return self._nodes_contact, self._segments_contact

    def select_nodes_for_gluing(
        self,
    ) -> Tuple[
        Dict[
            Union[int, str],
            Tuple[TargetPriority, Tuple[ParametricDirection, BoundarySide, List]],
        ],
        List,
    ]:
        raise NotImplementedError()

    def __repr__(self):
        message = f"""
            Boundary condition:
            - DOFs: {list(self._dofs_index.keys())}
            - type of function space: {set(self._dofs_funspace)}
        """
        return message
