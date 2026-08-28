from yeti_iga.pymfiga.common.base.enum import DOF, ParametricDirection, BoundarySide, TargetPriority
from yeti_iga.pymfiga.common.base.cls import BaseBoundaryCondition
from typing import List, Dict, Tuple, Union, Sequence, Literal
from itertools import product
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.BOUNDARY")


def create_connectivity_table(nnz_by_direction: list) -> np.ndarray:
    """
    Creates a connectivity table for control points based on the number of points in each direction.

    Args:
        nnz_by_direction (np.ndarray): Number of control points per direction.

    Returns:
        np.ndarray: Connectivity table.
    """
    assert isinstance(nnz_by_direction, list)
    indices = np.reshape(
        np.indices(nnz_by_direction), (len(nnz_by_direction), -1), order="F"
    )
    return np.transpose(indices).astype(int)


class BoundaryCondition(BaseBoundaryCondition):
    def __init__(self, nbctrlpts: Union[np.ndarray, List[int]], dofs: Sequence[DOF]):
        """
        Initializes the BoundaryCondition object.
        Args:
                        nbctrlpts (Union[np.ndarray, List[int]]):
                                Number of control points in each direction.
                        dofs (Sequence[DOF]):
                                Sequence of degrees of freedom to consider (e.g., DOF.X, DOF.Y, etc.).
        """
        logger.info("Initialize boundary conditions")
        self._set_dof_index(dofs)
        self._set_nbctrlpts(nbctrlpts)
        logger.info(f"Patch has {len(dofs)} degree of freedom per node")

        # Extra internal / public variables
        self._nbofdofs = {dof: np.prod(self._nbctrlpts) for dof in dofs}
        self._table_dirichlet: np.ndarray = np.zeros(
            shape=(len(dofs), self._ndim, 2), dtype=bool
        )
        self._connectivity_table: np.ndarray = create_connectivity_table(
            list(self._nbctrlpts)
        )
        self._data_of_ctrlpts_dirichlet: Dict[DOF, set] = {dof: set() for dof in dofs}
        self._data_of_ctrlpts_glued: Dict[
            Union[int, str],
            Tuple[TargetPriority, Tuple[ParametricDirection, BoundarySide, List[int]]],
        ] = {}
        self._boundary_dofs: Dict[Tuple[ParametricDirection, BoundarySide], list] = (
            self._set_boundary_nodes()
        )
        logger.info(repr(self))

    def _set_dof_index(self, dofs: Sequence[DOF]):
        assert isinstance(dofs, (list, tuple)) and len(dofs) > 0
        assert all(dof != DOF.ALL for dof in dofs)
        self._dofs_index: Dict[DOF, Tuple[int, Tuple[str, str]]] = {
            dof: (i, ("spline", "high order")) for i, dof in enumerate(dofs)
        }

    def _set_nbctrlpts(self, nbctrlpts: Union[np.ndarray, List[int]]):
        assert isinstance(nbctrlpts, (np.ndarray, list))
        nbctrlpts = np.asarray(nbctrlpts)
        nbctrlpts = nbctrlpts[nbctrlpts > 1]
        self._ndim: int = len(nbctrlpts)  # dimensionality in parametric space
        self._nbctrlpts: np.ndarray = nbctrlpts

    @property
    def dofs_index(self):
        return self._dofs_index

    @property
    def nbofdofs(self):
        return self._nbofdofs

    @property
    def table_dirichlet(self):
        return self._table_dirichlet

    @property
    def boundary_dofs(self):
        return self._boundary_dofs

    def _set_boundary_nodes(
        self,
    ) -> Dict[Tuple[ParametricDirection, BoundarySide], list]:
        nbctrlpts = self._nbctrlpts
        connectivity = self._connectivity_table
        loc_dir_list = [ParametricDirection(i) for i in range(self._ndim)]
        loc_face_list = [BoundarySide.MIN, BoundarySide.MAX]
        boundary_dofs = {}
        for loc_dir, loc_face in product(loc_dir_list, loc_face_list):
            idx_dir = loc_dir.value
            if loc_face.value == 0:
                nodes = sorted(np.where(connectivity[:, idx_dir] == 0)[0])
            else:
                nodes = sorted(
                    np.where(connectivity[:, idx_dir] == nbctrlpts[idx_dir] - 1)[0]
                )
            boundary_dofs[(loc_dir, loc_face)] = nodes
        return boundary_dofs

    def recognize_constraint(self, cnstr: Dict) -> Tuple[list, np.ndarray]:
        """
        Recognizes and marks constraints based on a list of constraints.
        """
        nbctrlpts = self._nbctrlpts
        connectivity = self._connectivity_table
        table = np.zeros(shape=(self._ndim, 2), dtype=bool)
        loc_dir = cnstr.get("direction")
        loc_face = cnstr.get("face")

        if not (
            isinstance(loc_dir, ParametricDirection)
            and isinstance(loc_face, BoundarySide)
        ):
            raise Warning(
                "Format style are not supported"
                f"Direction is {type(loc_dir)} and it should be ParametricDirection"
                f"Face is {type(loc_face)} and it should be BoundarySide"
            )

        idx_dir = loc_dir.value
        idx_fc = loc_face.value
        table[idx_dir, idx_fc] = True
        if idx_fc == 0:
            nodes = sorted(np.where(connectivity[:, idx_dir] == 0)[0])
        else:
            nodes = sorted(
                np.where(connectivity[:, idx_dir] == nbctrlpts[idx_dir] - 1)[0]
            )
        return nodes, table

    def _expand_constraint(self, cnstr: Dict) -> List[Dict]:
        """
        Expands a constraint dictionary into a list of dictionaries, one for each combination
        of direction and face, handling 'ALL' and 'BOTH' enum values.

        Args:
            cnstr (Dict): Dictionary with keys 'direction', 'face', 'dofs', 'target_id', 'target_priority'.

        Returns:
            List[Dict]: List of expanded constraints dictionaries for each direction-face combination.
        """
        directions = cnstr["direction"]
        faces = cnstr["face"]
        dofs = cnstr.get("dofs", ())
        target_id = cnstr.get("target_id", None)
        target_priority = cnstr.get("target_priority", None)

        if not isinstance(directions, (list, tuple)):
            directions = [directions]

        if not isinstance(faces, (list, tuple)):
            faces = [faces]

        expanded_faces = []
        for f in faces:
            if f == BoundarySide.BOTH:
                expanded_faces.extend([BoundarySide.MIN, BoundarySide.MAX])
            else:
                expanded_faces.append(f)

        expanded_directions = []
        for d in directions:
            if d == ParametricDirection.ALL:
                expanded_directions.extend(
                    [ParametricDirection(i) for i in range(self._ndim)]
                )
            else:
                expanded_directions.append(d)

        if dofs == DOF.ALL:
            dofs = self.dofs_index.keys()

        expanded = []
        for d, f in product(expanded_directions, expanded_faces):
            expanded.append(
                {
                    "direction": d,
                    "face": f,
                    "dofs": dofs,
                    "target_id": target_id,
                    "target_priority": target_priority,
                }
            )

        return expanded

    def _update_dirichlet_table(self, idx: int, table: np.ndarray):
        self._table_dirichlet[idx] = np.logical_or(self._table_dirichlet[idx], table)

    def _update_glued_nodes(
        self,
        target_id: Union[int, str],
        target_priority: TargetPriority,
        boundary_info: Tuple[ParametricDirection, BoundarySide, List[int]],
    ):
        """
        Updates the glued control points data structure.

        Args:
            target_id (Union[int, str]): Identifier for the glued set.
            target_priority (Priority): if target is master or slave in mortar methods
            target_info (Tuple): containts the direction and side of application,
                and the list of nodes affected
        """
        if target_id in self._data_of_ctrlpts_glued.keys():
            logger.warning(f"Overwriting glued nodes for target {target_id}")
        self._data_of_ctrlpts_glued.update(
            {target_id: (target_priority, boundary_info)}
        )

    def add_constraint(
        self, constraint_info: List[Dict], constraint_type: Literal["dirichlet", "glue"]
    ):
        """
        Adds a constraint to the boundary condition object.

        Args:
            constraint_list (List[Dict]): List of constraints.
            constraint_type (str): Type of constraint ("dirichlet", "glue", etc.).
        """
        t = constraint_type.lower()
        assert t in [
            "dirichlet",
            "glue",
        ], "Unknown constraint type"

        expanded_constraints: List[Dict] = []
        for loc in constraint_info:
            expanded_constraints.extend(self._expand_constraint(loc))

        if t == "dirichlet":
            nodes, table = set(), None
            for cnstr in expanded_constraints:
                dofs_list = cnstr.get("dofs")
                if not dofs_list:
                    raise ValueError("Empty dofs")

                nodes, table = self.recognize_constraint(cnstr)
                for dof in dofs_list:
                    idx = self._dofs_index[dof][0]
                    self._data_of_ctrlpts_dirichlet[dof].update(nodes)
                    if table is not None:
                        self._update_dirichlet_table(idx, table)

        elif t == "glue":
            nodes = set()
            for cnstr in expanded_constraints:
                target_id = cnstr.get("target_id")
                target_priority = cnstr.get("target_priority")
                if not isinstance(target_id, (int, str)):
                    raise ValueError("target_id must be int or str")
                if not isinstance(target_priority, TargetPriority):
                    raise ValueError("target_priority must be Priority")
                nodes = self.recognize_constraint(cnstr)[0]
                boundary_info = (cnstr["direction"], cnstr["face"], sorted(list(nodes)))
                self._update_glued_nodes(target_id, target_priority, boundary_info)

    def select_nodes_for_gluing(
        self,
    ) -> Tuple[
        Dict[
            Union[int, str],
            Tuple[TargetPriority, Tuple[ParametricDirection, BoundarySide, List]],
        ],
        List,
    ]:
        """
        Selects nodes for gluing (identification).
        """
        data_dict = self._data_of_ctrlpts_glued
        all_glued_nodes = []
        for val in data_dict.values():
            all_glued_nodes.extend(val[1][2])
        return data_dict, all_glued_nodes

    def select_nodes_for_solving(self) -> Tuple[List[int], List[int]]:
        """
        Selects free and constrained nodes for solving the system, taking into account
        the active DOFs.

        Returns:
            Tuple ([list, list]):
            Lists of free and constrained node indices.
        """

        def ravel(offset, list_of_nodes):
            flatten_list = []
            for idx in range(len(list_of_nodes)):
                flatten_list.extend(idx * offset + np.asarray(list_of_nodes[idx]))
            return flatten_list

        nbvars = len(self.dofs_index)
        free_nodes = [[] for _ in range(nbvars)]
        constraint_nodes = [[] for _ in range(nbvars)]
        total_nodes = np.prod(self._nbctrlpts)

        for dof, (idx, _) in self.dofs_index.items():
            blocked = self._data_of_ctrlpts_dirichlet[dof]
            free_nodes[idx] = sorted(set(range(total_nodes)).difference(blocked))
            constraint_nodes[idx] = sorted(blocked)

        return ravel(total_nodes, free_nodes), ravel(total_nodes, constraint_nodes)

    def select_nodes_for_contact(self):
        raise NotImplementedError()

    def __repr__(self):
        message = f"""
            Boundary condition:
            - DOFs: {list(self._dofs_index.keys())}
        """
        return message
