from yeti_iga.pymfiga.common.base.enum import DOF
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.common.material import LinearElasticity
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition
from .core import (
    SingleModel,
    dersbasisfun_linear_triangle,
    basisfun_quadratic_triangle,
    dersbasisfun_quadratic_triangle,
)
from .mitc_condensed import nedelec_1, nedelec_2, nedelec_3, basisfun_nedelec_triangle
from scipy import sparse as sp
from typing import Tuple
import numpy as np


class MITCModel(SingleModel):
    """
    Linear material
    - DZ: vertical displacement, linear
    - ROTX: rotation around x, quadratic
    - ROTY: rotation around y, quadratic
    - GAMMA: Lagrange multiplicator, linear
    - P: Lagrange multiplicator, linear
    """

    def __init__(
        self,
        mesh,
        thickness: float,
        material: LinearElasticity,
        quadrature: FiniteElementQuadrature,
        boundary: BoundaryCondition,
    ):
        super().__init__(mesh, material, quadrature, boundary)
        self._verify_mesh("quadratic")
        self._verify_boundary()

        # Internal variables
        self._thickness = thickness
        self._mass_matrix = None
        self._stiffness_matrix = None

    @property
    def material(self):
        assert isinstance(self._material, LinearElasticity)
        return self._material

    @property
    def thickness(self):
        return self._thickness

    def _verify_boundary(self):
        assert self.boundary.dofs_index[DOF.UZ][0] == 0
        assert self.boundary.dofs_index[DOF.ROTX][0] == 1
        assert self.boundary.dofs_index[DOF.ROTY][0] == 2
        assert self.boundary.dofs_index[DOF.LAG1][0] == 3
        assert self.boundary.dofs_index[DOF.LAG2][0] == 4
        assert self.boundary.dofs_index[DOF.UZ][1] == ("lagrange", 1)
        assert self.boundary.dofs_index[DOF.ROTX][1] == ("lagrange", 2)
        assert self.boundary.dofs_index[DOF.ROTY][1] == ("lagrange", 2)
        assert self.boundary.dofs_index[DOF.LAG1][1] == ("nedelec", 1)
        assert self.boundary.dofs_index[DOF.LAG2][1] == ("nedelec", 1)

    def clear_properties(self):
        if self.material.hasnonlinearmass:
            self._mass_matrix = None
        if self.material.hasnonlinearstiffness:
            self._stiffness_matrix = None

    def _assemble_bending_elem(self, elem_coords_quad) -> np.ndarray:
        # int grad theta : C : grad theta dx

        def isotropic_bending(obj: MITCModel):
            ndim = obj.part.ndim
            mat = obj.material
            t = obj.thickness
            e = mat.elastic_modulus
            nu = mat.poisson_ratio
            d = e * t**3 / (12.0 * (1.0 - nu**2))
            lame_lambda = d * nu
            lame_mu = d * (1 - nu) / 2
            lame_tensor = LinearElasticity.set_lame_tensor(lame_lambda, lame_mu, ndim)
            return lame_tensor

        ndim = self.part.ndim
        lagrange_dersfunbasis = self.quadrature.IntTriaQuad.evaluate_function(
            None, dersbasisfun_quadratic_triangle
        )[1]
        ten_rank4 = isotropic_bending(self)
        rowOfMatrices = []
        for i in range(ndim):
            columnOfMatrix = []
            for j in range(ndim):
                ten_rank2 = ten_rank4[i][j]
                det_jac, inv_jac = self.evaluate_jacobien_of_field(
                    lagrange_dersfunbasis, elem_coords_quad
                )
                prop = np.einsum(
                    "ilk,lm,jmk,k,k->ijk",
                    inv_jac,
                    ten_rank2,
                    inv_jac,
                    self.quadrature.IntTriaQuad.quadweights,
                    det_jac,
                )
                Melem: np.ndarray = np.einsum(
                    "kli,lmk,kmj->ij",
                    lagrange_dersfunbasis,
                    prop,
                    lagrange_dersfunbasis,
                )
                columnOfMatrix.append(Melem)
            rowOfMatrices.append(columnOfMatrix)
        return np.block(rowOfMatrices)

    def _assemble_shear_elem(self, elem_coords_line) -> np.ndarray:
        # int gamma . S gamma dx
        # With Nedelec functions

        def isotropic_shear(obj: MITCModel):
            mat = obj.material
            t = obj.thickness
            e = mat.elastic_modulus
            nu = mat.poisson_ratio
            k = mat.timoshenko_ratio
            return e * k * t / (2.0 * (1.0 + nu))

        nedelec_funbasis = self.quadrature.IntTriaQuad.evaluate_function(
            basisfun_nedelec_triangle, None
        )[0]
        lagrange_dersfunbasis = self.quadrature.IntTriaQuad.evaluate_function(
            None, dersbasisfun_linear_triangle
        )[1]
        det_jac, inv_jac = self.evaluate_jacobien_of_field(
            lagrange_dersfunbasis, elem_coords_line
        )
        prop = isotropic_shear(self) * np.einsum(
            "ilk,jlk,k,k->ijk",
            inv_jac,
            inv_jac,
            self.quadrature.IntTriaQuad.quadweights,
            det_jac,
        )
        Melem: np.ndarray = np.einsum(
            "kli,lmk,kmj->ij", nedelec_funbasis, prop, nedelec_funbasis
        )
        return Melem

    def _assemble_multiplier_elem(
        self,
        info_edge: dict,
        elements_quad: list,
        points_quad: list,
        elements_line: list,
        points_line: list,
    ) -> dict:
        def assemble_grad_w_p(ders_lag, nedelec, e, scaling=1.0):
            grad_lagrange_dot_tangent = np.einsum("kji,j->ki", ders_lag, TANGENT[e])
            nedelec_dot_tangent = np.einsum("kj,j->k", nedelec, TANGENT[e])
            grad_lagrange_nedelec = scaling * np.einsum(
                "ki,k,k->i",
                grad_lagrange_dot_tangent,
                nedelec_dot_tangent,
                self.quadrature.ExtTriaQuad.quadweights[e],
            )
            return grad_lagrange_nedelec

        def assemble_theta_p(jac, lag, nedelec, e, scaling=1.0):
            lagrange_nedelec = []
            for d in range(self.part.ndim):
                jac_tan = np.einsum("ijk,j->ik", jac, TANGENT[e])
                lagrange_dot_tangent = np.einsum("ki,k->ki", lag, jac_tan[d])
                nedelec_dot_tangent = np.einsum("kj,j->k", nedelec, TANGENT[e])
                lagrange_i_nedelec = -scaling * np.einsum(
                    "ki,k,k->i",
                    lagrange_dot_tangent,
                    nedelec_dot_tangent,
                    self.quadrature.ExtTriaQuad.quadweights[e],
                )
                lagrange_nedelec.extend(lagrange_i_nedelec)
            return np.array(lagrange_nedelec)

        def assemble_gamma_p(nedelec, e, scaling=1.0):
            nedelec_dot_tangent = np.einsum("kj,j->k", nedelec, TANGENT[e])
            nedelec_nedelec = -scaling * np.einsum(
                "k,k,k",
                nedelec_dot_tangent,
                nedelec_dot_tangent,
                self.quadrature.ExtTriaQuad.quadweights[e],
            )
            return nedelec_nedelec

        # Constants at the triangle of reference
        EDGE_LENGTH = np.array([1.0, np.sqrt(2), 1.0])

        TANGENT = np.array(
            [[1.0, 0.0], [-1.0 / np.sqrt(2), 1.0 / np.sqrt(2)], [0, -1.0]]
        )

        NEDELEC_BASIS = [nedelec_1, nedelec_2, nedelec_3]

        # Edge properties
        global_nodes = info_edge["nodes"]
        global_edge_tangent = np.array(points_line[global_nodes[0]]) - np.array(
            points_line[global_nodes[1]]
        )
        global_edge_length = np.linalg.norm(global_edge_tangent)

        # There are at most two elements sharing the edge
        # We decide to work with the first element in the list
        # which is the "plus" side
        elem_global_idx, edge_local_idx, edge_orientation = info_edge["elements"][0]
        global_idx_lag_1 = elements_line[elem_global_idx]
        global_idx_lag_2 = elements_quad[elem_global_idx]

        # Recover element coordinates
        elem_coords_line = np.array([points_line[_] for _ in global_idx_lag_1])
        elem_coords_quad = np.array([points_quad[_] for _ in global_idx_lag_2])

        # ---------------------------------------------
        # Add assertions to verify correct orientation
        local_edge_tangents = [
            elem_coords_line[j] - elem_coords_line[i]
            for i, j in zip([1, 2, 0], [0, 1, 2])
        ]
        local_edge_lengths = [np.linalg.norm(local_edge_tangents[i]) for i in range(3)]
        local_orientation = np.sign(
            np.dot(global_edge_tangent, local_edge_tangents[edge_local_idx])
        )
        assert np.isclose(local_edge_lengths[edge_local_idx], global_edge_length)
        assert np.isclose(local_orientation, edge_orientation)
        # ---------------------------------------------

        # Assemble matrices
        quad = self.quadrature.ExtTriaQuad
        NEDELEC = quad.evaluate_function(
            NEDELEC_BASIS[edge_local_idx], None, edge_local_idx
        )[0]
        NEDELEC *= edge_orientation  # Adjust orientation
        DERS_LAG_LINE = quad.evaluate_function(
            None, dersbasisfun_linear_triangle, edge_local_idx
        )[1]
        LAG_QUAD, DERS_LAG_QUAD = quad.evaluate_function(
            basisfun_quadratic_triangle, dersbasisfun_quadratic_triangle, edge_local_idx
        )
        JAC_QUAD = np.einsum("kil,lj->jik", DERS_LAG_QUAD, elem_coords_quad)

        factor = EDGE_LENGTH[edge_local_idx] / local_edge_lengths[edge_local_idx]
        w_p_col = assemble_grad_w_p(DERS_LAG_LINE, NEDELEC, edge_local_idx, factor)
        th_p_col = assemble_theta_p(JAC_QUAD, LAG_QUAD, NEDELEC, edge_local_idx, factor)
        C_local = np.block([w_p_col, th_p_col])
        D_edge = assemble_gamma_p(NEDELEC, edge_local_idx, factor)

        offset_list = [
            self.boundary.limits_dofs[_][0] for _ in [DOF.UZ, DOF.ROTX, DOF.ROTY]
        ]
        global_indices_z = []
        for offset, indices in zip(
            offset_list,
            [
                global_idx_lag_1,
                global_idx_lag_2,
                global_idx_lag_2,
            ],
        ):
            global_indices_z.extend([idx + offset for idx in indices])

        return {
            "C_plus": C_local,
            "indices_z_plus": global_indices_z,
            "D_edge": D_edge,
        }

    def assemble_stiffness(self) -> sp.csr_array:

        ################# INTEGRAL OVER ELEMENTS #################
        offset_list = [
            self.boundary.limits_dofs[_][0]
            for _ in [DOF.UZ, DOF.ROTX, DOF.ROTY, DOF.LAG1]
        ]

        elements_line, points_line = self._recover_mesh_property(
            order="linear", space="lagrange"
        )
        elements_quad, points_quad = self._recover_mesh_property(
            order="quadratic", space="lagrange"
        )
        elements_edge = self._recover_mesh_property(order="linear", space="nedelec")[0]
        edge_connectivity = self.part.get_edge_connectivity()

        # Ensamblar bloques AA y BB (diagonales superiores)
        row_ind, col_ind, val_ind = [], [], []
        for elem_id, (
            global_idx_lag_1,
            global_idx_lag_2,
            global_idx_ned,
        ) in enumerate(zip(elements_line, elements_quad, elements_edge)):
            global_indices = []
            for offset, indices in zip(
                offset_list,
                [
                    global_idx_lag_1,
                    global_idx_lag_2,
                    global_idx_lag_2,
                    global_idx_ned,
                ],
            ):
                global_indices.extend([idx + offset for idx in indices])

            elem_coords_line = np.array([points_line[_] for _ in global_idx_lag_1])
            elem_coords_quad = np.array([points_quad[_] for _ in global_idx_lag_2])
            # Get orientation of the nedelec functions
            orientation_list = []
            for idx in global_idx_ned:
                edge_info = edge_connectivity[idx]
                for elem_data in edge_info["elements"]:
                    if elem_data[0] == elem_id:
                        orientation_list.append(elem_data[2])
                        break
            assert len(orientation_list) == len(global_idx_ned)
            orientation_matrix = np.outer(orientation_list, orientation_list)

            th_th = self._assemble_bending_elem(elem_coords_quad)
            gm_gm = self._assemble_shear_elem(elem_coords_line)
            w_w = np.zeros((3, 3))
            w_th = np.zeros((3, 12))
            z_gm = np.zeros((15, 3))
            AA = np.block([[w_w, w_th], [w_th.T, th_th]])
            BB = gm_gm * orientation_matrix
            Aelem = np.block([[AA, z_gm], [z_gm.T, BB]])

            for ii in range(Aelem.shape[0]):
                for jj in range(Aelem.shape[1]):
                    row_ind.append(global_indices[ii])
                    col_ind.append(global_indices[jj])
                    val_ind.append(Aelem[ii, jj])

        size_mat_total = self.get_size_of_arrays()
        matrix_element = sp.coo_array(
            (val_ind, (row_ind, col_ind)), shape=(size_mat_total, size_mat_total)
        )
        matrix_element.eliminate_zeros()

        ################# INTEGRAL OVER EDGES #################
        edge_connectivity = self.part.get_edge_connectivity()
        offset_lag1 = self.boundary.limits_dofs[DOF.LAG1][0]
        offset_lag2 = self.boundary.limits_dofs[DOF.LAG2][0]
        row_ind, col_ind, val_ind = [], [], []
        for edge_global_idx, edge_info in edge_connectivity.items():
            edge_data = self._assemble_multiplier_elem(
                edge_info,
                elements_quad,
                points_quad,
                elements_line,
                points_line,
            )

            C_plus = np.ravel(edge_data["C_plus"])
            indices_plus = edge_data["indices_z_plus"]
            idx_p = [edge_global_idx + offset_lag2] * len(C_plus)

            row_ind.extend(idx_p)
            col_ind.extend(indices_plus)
            val_ind.extend(C_plus)
            row_ind.extend(indices_plus)
            col_ind.extend(idx_p)
            val_ind.extend(C_plus)

            idx_p = edge_global_idx + offset_lag2
            idx_gamma = edge_global_idx + offset_lag1
            D_edge = edge_data["D_edge"]

            row_ind.append(idx_gamma)
            col_ind.append(idx_p)
            val_ind.append(D_edge)
            row_ind.append(idx_p)
            col_ind.append(idx_gamma)
            val_ind.append(D_edge)

        matrix_edge = sp.coo_array(
            (val_ind, (row_ind, col_ind)), shape=(size_mat_total, size_mat_total)
        )
        matrix_edge.eliminate_zeros()
        mat_final = matrix_element + matrix_edge
        return sp.csr_array(mat_final.tocsr())

    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        tangent_matrix = self.assemble_stiffness(**kwargs)
        free_nodes = self.get_free_and_constraint_nodes()[0]
        increment = np.zeros_like(array_in)
        increment[free_nodes] = LinearSolver.direct(
            tangent_matrix[np.ix_(free_nodes, free_nodes)], array_in[free_nodes]
        )["sol"]

        if self.update_manager.should_update_material:
            self.clear_properties()

        return increment

    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        return np.array([]), {}
