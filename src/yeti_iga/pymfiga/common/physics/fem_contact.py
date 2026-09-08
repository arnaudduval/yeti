from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition
from yeti_iga.pymfiga.fem.model.mechanical import MechanicalModel
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from .core import Physics, clear_bcs
from dataclasses import dataclass
from typing import Tuple, Sequence
from scipy import sparse as sp
from copy import deepcopy
import numpy as np
import logging

logger = logging.getLogger("SRC.PHYSICS")


def convert_array_to_matrix(array, dim=2):
    return np.reshape(array, (dim, -1))


@dataclass
class ContactArgs:
    normal_penalty: float = 10 * Constants.PENALTY
    tangent_penalty: float = 10 * Constants.PENALTY
    friction: float = 0.0

    def __post_init__(self):
        if self.normal_penalty <= 0.0:
            raise ValueError("Penalty should be positive")
        if self.tangent_penalty <= 0.0:
            raise ValueError("Penalty should be positive")
        if self.friction < 0:
            raise ValueError("Friction should be > 0")


class ContactProblem(Physics):
    def __init__(
        self,
        model_master: MechanicalModel,
        model_slave: MechanicalModel,
        contact_args: dict = {},
        **solver_args,
    ):
        super().__init__(**solver_args)
        self._verify_model(model_master)
        self._verify_model(model_slave)
        self.master = model_master
        self.slave = model_slave
        self._config = ContactArgs(**contact_args)
        self.master.set_update_manager(self.update_manager)
        self.slave.set_update_manager(self.update_manager)

    def _verify_model(self, model: MechanicalModel):
        assert (
            hasattr(model.part, "recover_points")
            and hasattr(model.part, "mesh_order")
            and hasattr(model, "assemble_stiffness")
        ), AttributeError()
        assert model.part.mesh_order == "linear"

    def _compute_contact_coupling(
        self,
        contact_boundary_m: BoundaryCondition,
        contact_boundary_s: BoundaryCondition,
        disp_m_1: np.ndarray,
        disp_m_0: np.ndarray,
        disp_s_1: np.ndarray,
        disp_s_0: np.ndarray,
        size_of_blocks: Sequence[int],
    ) -> Tuple[np.ndarray, sp.csr_array]:

        master_points = self.master.part.recover_points("lagrange", 1)
        slave_points = self.slave.part.recover_points("lagrange", 1)

        size_total_nodes = np.sum(size_of_blocks)
        size_nodes_m = len(master_points)
        size_nodes_s = len(slave_points)
        vector = np.zeros(size_total_nodes)
        row_ind, col_ind, values = [], [], []

        idx_segments_m = contact_boundary_m.select_nodes_for_contact()[-1]
        idx_nodes_s = contact_boundary_s.select_nodes_for_contact()[0]

        assert idx_segments_m is not None, "Define contact boundary"
        assert idx_nodes_s is not None, "Define contact boundary"

        for idx_nodes_segment_m in idx_segments_m:
            # NOTE: by default the indices are noted counter clockwise
            # Then, following this orientation, the solid is on the left of the segment
            # However, cntelm2d assumes that the solid is on the right
            # That is why it is necessay to flip the indices
            reversed_index: list = list(idx_nodes_segment_m)
            reversed_index.reverse()
            segment_coords_m = np.array(
                [master_points[idx_nd] for idx_nd in reversed_index]
            )
            segment_disp_m_0 = np.array(
                [disp_m_0[:, idx_nd] for idx_nd in reversed_index]
            )
            segment_disp_m_1 = np.array(
                [disp_m_1[:, idx_nd] for idx_nd in reversed_index]
            )

            for idx_nd_s in idx_nodes_s:
                node_coords_s = np.array(slave_points[idx_nd_s])
                node_disp_s_0 = np.array(disp_s_0[:, idx_nd_s])
                node_disp_s_1 = np.array(disp_s_1[:, idx_nd_s])

                elxy_0 = np.transpose(
                    np.vstack(
                        (
                            node_coords_s + node_disp_s_0,
                            segment_coords_m + segment_disp_m_0,
                        )
                    )
                )
                elxy_1 = np.transpose(
                    np.vstack(
                        (
                            node_coords_s + node_disp_s_1,
                            segment_coords_m + segment_disp_m_1,
                        )
                    )
                )

                vec_local, mat_local = cntelm2d(
                    ELXY=elxy_1,
                    ELXYP=elxy_0,
                    OMEGAN=self._config.normal_penalty,
                    OMEGAT=self._config.tangent_penalty,
                    CFRI=self._config.friction,
                )

                if vec_local is None or mat_local is None:
                    continue

                idx_list = [
                    idx_nd_s + size_of_blocks[0],
                    idx_nd_s + size_of_blocks[0] + size_nodes_s,
                    idx_nodes_segment_m[0],
                    idx_nodes_segment_m[0] + size_nodes_m,
                    idx_nodes_segment_m[1],
                    idx_nodes_segment_m[1] + size_nodes_m,
                ]

                for i in range(len(idx_list)):
                    vector[idx_list[i]] += vec_local[i]
                    for j in range(len(idx_list)):
                        row_ind.append(idx_list[i])
                        col_ind.append(idx_list[j])
                        values.append(mat_local[i, j])

        matrix = sp.coo_array(
            (values, (row_ind, col_ind)), shape=(size_total_nodes, size_total_nodes)
        ).tocsr()
        matrix.eliminate_zeros()
        return vector, sp.csr_array(matrix)

    def _compute_residual(self, array_in: np.ndarray, **kwargs):
        size_block_m = self.master.get_size_of_arrays()
        size_block_s = self.slave.get_size_of_arrays()
        free_nodes_m = self.master.get_free_and_constraint_nodes()[0]
        free_nodes_s = self.slave.get_free_and_constraint_nodes()[0]
        free_nodes = np.hstack(
            (np.array(free_nodes_m), np.array(free_nodes_s) + size_block_m)
        )

        # Extract data
        force_contact: np.ndarray = kwargs["force_contact"]
        matrix_contact: sp.csr_array = kwargs["matrix_contact"]
        old_displacement: np.ndarray = kwargs["old_displacement"]
        res_args_m: dict = kwargs["master"]
        res_args_s: dict = kwargs["slave"]

        # Compute residual (without contact)
        dj_n1_m = array_in[:size_block_m]
        dj_n1_s = array_in[size_block_m:]
        res_m, extra_args_m = self.master.compute_residual(dj_n1_m, **res_args_m)
        res_s, extra_args_s = self.slave.compute_residual(dj_n1_s, **res_args_s)
        residual = np.zeros_like(array_in)
        residual[free_nodes] = np.hstack((res_m[free_nodes_m], res_s[free_nodes_s]))

        # Compute penalization due to contact
        incr_force_contact, incr_stiff_contact = self._compute_contact_coupling(
            self.master.boundary,
            self.slave.boundary,
            convert_array_to_matrix(dj_n1_m),
            convert_array_to_matrix(old_displacement[:size_block_m]),
            convert_array_to_matrix(dj_n1_s),
            convert_array_to_matrix(old_displacement[size_block_m:]),
            size_of_blocks=(size_block_m, size_block_s),
        )
        new_matrix_contact = matrix_contact + incr_stiff_contact
        new_force_contact = force_contact + incr_force_contact

        # Add penalization from contact
        residual += new_force_contact

        return residual, {
            "master": extra_args_m,
            "slave": extra_args_s,
            "force_contact": new_force_contact,
            "matrix_contact": new_matrix_contact,
        }

    def _solve_linearized_system(self, array_in: np.ndarray, **kwargs):
        size_block_m = self.master.get_size_of_arrays()
        size_block_s = self.slave.get_size_of_arrays()
        free_nodes_m = self.master.get_free_and_constraint_nodes()[0]
        free_nodes_s = self.slave.get_free_and_constraint_nodes()[0]
        free_nodes = np.hstack(
            (np.array(free_nodes_m), np.array(free_nodes_s) + size_block_m)
        )

        # Compute tangent matrix
        matrix_contact = kwargs["matrix_contact"]
        tangent_m = self.master.assemble_stiffness(**kwargs["master"])
        tangent_s = self.slave.assemble_stiffness(**kwargs["slave"])
        tangent_matrix = sp.bmat(
            [
                [tangent_m, sp.csr_array((size_block_m, size_block_s))],
                [sp.csr_array((size_block_s, size_block_m)), tangent_s],
            ]
        )

        # Add penalization from contact
        tangent_matrix += matrix_contact

        # Solve
        increment = np.zeros_like(array_in)
        increment[free_nodes] = LinearSolver.direct(
            tangent_matrix[np.ix_(free_nodes, free_nodes)], array_in[free_nodes]
        )["sol"]
        return increment

    def solve(
        self,
        disp_list_m: np.ndarray,
        disp_list_s: np.ndarray,
        force_list_m: np.ndarray,
        force_list_s: np.ndarray,
    ):
        clear_bcs(self.master, disp_list_m, force_list_m)
        clear_bcs(self.slave, disp_list_s, force_list_s)
        self.master.set_update_manager(self.update_manager)
        self.slave.set_update_manager(self.update_manager)

        size_block_m = self.master.get_size_of_arrays()
        size_block_s = self.slave.get_size_of_arrays()

        total_size = size_block_m + size_block_s
        force_contact = np.zeros(total_size)
        matrix_contact = sp.csr_array((total_size, total_size))

        # Internal variables
        plastic_vars_m = {}
        plastic_vars_s = {}

        for it in range(1, np.shape(force_list_m)[0]):

            logger.info(f"(Pseudo) Time-step: {it}")
            self.update_manager.increment_step()

            # Assemble
            new_displacement = np.hstack((disp_list_m[it], disp_list_s[it]))
            old_displacement = np.hstack((disp_list_m[it - 1], disp_list_s[it - 1]))

            res_args = {
                "master": {
                    "external_force": np.copy(force_list_m[it]),
                    "plastic_vars": plastic_vars_m,
                },
                "slave": {
                    "external_force": np.copy(force_list_s[it]),
                    "plastic_vars": plastic_vars_s,
                },
                "old_displacement": old_displacement,
                "force_contact": force_contact,
                "matrix_contact": matrix_contact,
            }

            all_extra_args = self.nonlinear_solver.solve(
                new_displacement,
                self._compute_residual,
                self._solve_linearized_system,
                residual_args=res_args,
                increment_args={},
            )
            extra_args = all_extra_args["extra_args"]

            # Save data for next step
            disp_list_m[it] = new_displacement[:size_block_m]
            disp_list_s[it] = new_displacement[size_block_m:]
            force_contact = extra_args["force_contact"].copy()
            matrix_contact = extra_args["matrix_contact"].copy()
            plastic_vars_m = deepcopy(extra_args["master"]["new_plastic_vars"])
            plastic_vars_s = deepcopy(extra_args["slave"]["new_plastic_vars"])

        return plastic_vars_m, plastic_vars_s


def cntelm2d(ELXY, ELXYP, OMEGAN, OMEGAT, CFRI, LTAN=True):
    # ********************************************************************
    # SEARCH CONTACT POINT AND RETURN STIFFNESS AND RESIDUAL FORCE
    # IF CONTACTED FOR NORMAL CONTACT
    # ********************************************************************
    #
    ZERO = 0.0
    ONE = 1.0
    EPS = 1e-6
    P05 = 0.05
    FORCE = None
    STIFF = None
    XT = ELXY[:, 2] - ELXY[:, 1]
    XLEN = np.linalg.norm(XT)
    if XLEN < EPS:
        return FORCE, STIFF
    XTP = ELXYP[:, 2] - ELXYP[:, 1]
    XLENP = np.linalg.norm(XTP)
    #
    # UNIT NORMAL AND TANGENTIAL VECTOR
    XT = XT / XLEN
    XTP = XTP / XLENP
    XN = np.array([-XT[1], XT[0]])
    # NORMAL GAP FUNCTION Gn = (X_s - X_1).N
    GAPN = np.dot((ELXY[:, 0] - ELXY[:, 1]), XN)
    if GAPN == 0.0:
        GAPN = -EPS * XLEN
    #
    # CHECK IMPENETRATION CONDITION
    if (GAPN > ZERO) or (GAPN <= -XLEN):
        return FORCE, STIFF
    #
    # NATURAL COORDINATE AT CONTACT POINT
    ALPHA = np.dot((ELXY[:, 0] - ELXY[:, 1]), XT) / XLEN
    ALPHA0 = np.dot((ELXYP[:, 0] - ELXYP[:, 1]), XTP) / XLENP
    #
    # OUT OF SEGMENT
    if (ALPHA > ONE + P05) or (ALPHA < -P05):
        return FORCE, STIFF
    #
    # CONTACT OCCURS IN THIS SEGMENT
    XLAMBN = -OMEGAN * GAPN
    XLAMBT = 0
    LFRIC = 1
    GAPT = 0
    LSLIDE = None
    if CFRI == 0:
        LFRIC = 0
    if LFRIC:
        GAPT = (ALPHA - ALPHA0) * XLENP
        XLAMBT = -OMEGAT * GAPT
        FRTOL = XLAMBN * CFRI
        LSLIDE = 0
        if abs(XLAMBT) > FRTOL:
            LSLIDE = 1
            XLAMBT = -FRTOL * np.sign(GAPT)
    #
    # DEFINE VECTORS
    NN = np.hstack((XN, -(ONE - ALPHA) * XN, -ALPHA * XN))
    TT = np.hstack((XT, -(ONE - ALPHA) * XT, -ALPHA * XT))
    PP = np.hstack((np.zeros_like(XN), -XN, XN))
    QQ = np.hstack((np.zeros_like(XT), -XT, XT))
    CN = NN - GAPN * QQ / XLEN
    CT = TT + GAPN * PP / XLEN
    #
    # CONTACT FORCE
    FORCE = XLAMBN * CN + XLAMBT * CT
    #
    # FORM STIFFNESS
    if LTAN:
        STIFF = OMEGAN * np.outer(CN, CN)
        if LFRIC:
            TMP1 = -CFRI * OMEGAN * np.sign(GAPT)
            TMP2 = -XLAMBT / XLEN
            if LSLIDE:
                STIFF = (
                    STIFF
                    + TMP1 * np.outer(CT, CN)
                    + TMP2
                    * (
                        np.outer(CT, PP)
                        + np.outer(PP, CN)
                        - np.outer(CT, QQ)
                        - np.outer(QQ, CT)
                    )
                )
            else:
                STIFF = (
                    STIFF
                    + OMEGAT * np.outer(CT, CT)
                    + TMP2
                    * (
                        np.outer(CN, PP)
                        + np.outer(PP, CN)
                        - np.outer(CT, QQ)
                        - np.outer(QQ, CT)
                    )
                )

    return FORCE, STIFF
