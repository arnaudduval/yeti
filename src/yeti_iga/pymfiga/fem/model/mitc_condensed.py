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
from typing import Sequence, Tuple
from scipy import sparse as sp
import numpy as np


def nedelec_1(pts: np.ndarray):
    xi, eta = pts
    return [1 - eta, xi]


def nedelec_2(pts: np.ndarray):
    xi, eta = pts
    return [-eta, xi]


def nedelec_3(pts: np.ndarray):
    xi, eta = pts
    return [-eta, xi - 1]


def basisfun_nedelec_triangle(pts: np.ndarray) -> np.ndarray:
    return np.array([nedelec_1(pts), nedelec_2(pts), nedelec_3(pts)]).T


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
        assert self.boundary.dofs_index[DOF.UZ][1] == ("lagrange", 1)
        assert self.boundary.dofs_index[DOF.ROTX][1] == ("lagrange", 2)
        assert self.boundary.dofs_index[DOF.ROTY][1] == ("lagrange", 2)

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
            fun=None,
            dersfun=dersbasisfun_quadratic_triangle,
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
        elem_coords_line,
        elem_coords_quad,
    ) -> Sequence[np.ndarray]:

        EDGE_LENGTH = np.array([1.0, np.sqrt(2), 1.0])

        TANGENT = np.array(
            [[1.0, 0.0], [-1.0 / np.sqrt(2), 1.0 / np.sqrt(2)], [0, -1.0]]
        )

        nedelec_basis = [nedelec_1, nedelec_2, nedelec_3]

        edge_length = [
            np.linalg.norm(elem_coords_line[j] - elem_coords_line[i])
            for i, j in zip([1, 2, 0], [0, 1, 2])
        ]

        quad = self.quadrature.ExtTriaQuad

        output_1, output_2, output_3 = [], [], []
        for e, fun in enumerate(nedelec_basis):

            # Get basis functions on edge
            NEDELEC = quad.evaluate_function(fun, None, e)[0]

            DERS_LAG_LINE = quad.evaluate_function(
                None, dersbasisfun_linear_triangle, e
            )[1]

            LAG_QUAD, DERS_LAG_QUAD = quad.evaluate_function(
                basisfun_quadratic_triangle, dersbasisfun_quadratic_triangle, e
            )

            # Length of reference edge / Length of physical edge
            factor = EDGE_LENGTH[e] / edge_length[e]

            # Precompute
            JAC_QUAD = np.einsum("kil,lj->jik", DERS_LAG_QUAD, elem_coords_quad)
            nedelec_dot_tangent = np.einsum("kj,j->k", NEDELEC, TANGENT[e])
            jac_tangent = np.einsum("ijk,j->ik", JAC_QUAD, TANGENT[e])

            # assemble grad_lagrange . nedelec
            grad_lagrange_dot_tangent = np.einsum(
                "kji,j->ki", DERS_LAG_LINE, TANGENT[e]
            )
            grad_w_nedelec = factor * np.einsum(
                "ki,k,k->i",
                grad_lagrange_dot_tangent,
                nedelec_dot_tangent,
                self.quadrature.ExtTriaQuad.quadweights[e],
            )

            # assemble [lagrange, lagrange] . nedelec
            theta_nedelec = []
            for d in range(self.part.ndim):
                lagrange_dot_jac_tangent = np.einsum(
                    "ki,k->ki", LAG_QUAD, jac_tangent[d]
                )
                theta_i_nedelec = -factor * np.einsum(
                    "ki,k,k->i",
                    lagrange_dot_jac_tangent,
                    nedelec_dot_tangent,
                    self.quadrature.ExtTriaQuad.quadweights[e],
                )
                theta_nedelec.extend(theta_i_nedelec)
            theta_nedelec = np.array(theta_nedelec)

            # Assemble nedelec . nedelec
            nedelec_nedelec = -factor * np.einsum(
                "k,k,k",
                nedelec_dot_tangent,
                nedelec_dot_tangent,
                self.quadrature.ExtTriaQuad.quadweights[e],
            )

            output_1.append(grad_w_nedelec)
            output_2.append(theta_nedelec)
            output_3.append(nedelec_nedelec)

        return [np.array(output_1).T, np.array(output_2).T, np.diag(output_3)]

    def assemble_stiffness(self) -> sp.csr_array:

        offset_list = [
            self.boundary.limits_dofs[_][0] for _ in [DOF.UZ, DOF.ROTX, DOF.ROTY]
        ]
        (
            elements_line,
            points_line,
        ) = self._recover_mesh_property(order="linear", space="lagrange")
        (
            elements_quad,
            points_quad,
        ) = self._recover_mesh_property(order="quadratic", space="lagrange")
        row_ind, col_ind, val_ind = [], [], []
        for i, (global_idx_lag_1, global_idx_lag_2) in enumerate(
            zip(elements_line, elements_quad)
        ):
            global_idx_z = []
            for offset, indices in zip(
                offset_list,
                [
                    global_idx_lag_1,
                    global_idx_lag_2,
                    global_idx_lag_2,
                ],
            ):
                global_idx_z.extend([idx + offset for idx in indices])

            elem_coords_line = np.array([points_line[_] for _ in global_idx_lag_1])
            elem_coords_quad = np.array([points_quad[_] for _ in global_idx_lag_2])

            th_th: np.ndarray = self._assemble_bending_elem(elem_coords_quad)
            gm_gm: np.ndarray = self._assemble_shear_elem(elem_coords_line)
            w_p, th_p, gm_p = self._assemble_multiplier_elem(
                elem_coords_line, elem_coords_quad
            )
            w_w = np.zeros((3, 3))
            w_th = np.zeros((3, 12))

            AA = np.block([[w_w, w_th], [w_th.T, th_th]])
            BB = np.copy(gm_gm)
            CC = np.block([[w_p], [th_p]])
            DD_inv = np.diag(1.0 / np.diag(gm_p))
            Aelem = AA + CC @ DD_inv @ BB @ DD_inv @ CC.T

            for i in range(Aelem.shape[0]):
                for j in range(Aelem.shape[1]):
                    row_ind.append(global_idx_z[i])
                    col_ind.append(global_idx_z[j])
                    val_ind.append(Aelem[i, j])

        size_mat_tota = self.get_size_of_arrays()
        matrix = sp.coo_array(
            (val_ind, (row_ind, col_ind)), shape=(size_mat_tota, size_mat_tota)
        )
        matrix.eliminate_zeros()
        return sp.csr_array(matrix.tocsr())

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
