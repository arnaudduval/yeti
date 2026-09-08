from yeti_iga.pymfiga.common.base.enum import DOF
from yeti_iga.pymfiga.common.base.math import eval_inverse_and_determinant
from yeti_iga.pymfiga.common.base.cls import BaseSingleModel
from yeti_iga.pymfiga.common.material import Material
from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.common.numerics.solvers import UpdateManager
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition
from yeti_iga.pymfiga.fem.geometry import MeshConverter
from typing import Dict, Callable, Tuple, Sequence
from scipy import sparse as sp
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.FEA.MODEL")


def basisfun_linear_line(pts: np.ndarray) -> np.ndarray:
    return np.array([(1 - pts) / 2, (1 + pts) / 2])


def dersbasisfun_linear_line(pts: np.ndarray) -> np.ndarray:
    dN_dxi = np.array([-0.5, 0.5])
    return dN_dxi


def basisfun_quadratic_line(pts: np.ndarray) -> np.ndarray:
    return np.array([-pts * (1 - pts) / 2, (1 - pts) * (1 + pts), pts * (1 + pts) / 2])


def dersbasisfun_quadratic_line(pts: np.ndarray) -> np.ndarray:
    dN_dxi = np.array([(2 * pts - 1) / 2, (-2 * pts), (2 * pts + 1) / 2])
    return dN_dxi


def basisfun_linear_triangle(pts: np.ndarray) -> np.ndarray:
    return np.array([1 - pts[0] - pts[1], pts[0], pts[1]])


def dersbasisfun_linear_triangle(pts: np.ndarray) -> np.ndarray:
    dN_dxi = np.array([-1, 1, 0])
    dN_deta = np.array([-1, 0, 1])
    return np.vstack((dN_dxi, dN_deta))


def basisfun_quadratic_triangle(pts: np.ndarray) -> np.ndarray:
    xi, eta = pts
    return np.array(
        [
            (1 - xi - eta) * (1 - 2 * xi - 2 * eta),
            xi * (2 * xi - 1),
            eta * (2 * eta - 1),
            4 * xi * (1 - xi - eta),
            4 * xi * eta,
            4 * eta * (1 - xi - eta),
        ]
    )


def dersbasisfun_quadratic_triangle(pts: np.ndarray) -> np.ndarray:
    xi, eta = pts
    dN_dxi = np.array(
        [
            -(1 - 2 * xi - 2 * eta) - 2 * (1 - xi - eta),
            4 * xi - 1,
            0,
            4 - 8 * xi - 4 * eta,
            4 * eta,
            -4 * eta,
        ]
    )
    dN_deta = np.array(
        [
            -(1 - 2 * xi - 2 * eta) - 2 * (1 - xi - eta),
            0,
            4 * eta - 1,
            -4 * xi,
            4 * xi,
            4 - 8 * eta - 4 * xi,
        ]
    )
    return np.vstack((dN_dxi, dN_deta))


class SingleModel(BaseSingleModel):
    def __init__(
        self,
        mesh: MeshConverter,
        material: Material,
        quadrature: FiniteElementQuadrature,
        boundary: BoundaryCondition,
    ):
        self.TypeProblemToBeApplied = "FEA"

        # Public variables
        self._set_part(mesh)
        self._set_material(material)
        self._set_quadrature(quadrature)
        self._set_boundary(boundary)

        # Internal variables
        output = boundary.select_nodes_for_solving()
        self._free_nodes, self._constraint_nodes = output
        self._update_manager = None

    @property
    def part(self):
        return self._part

    def _set_part(self, value):
        assert isinstance(value, MeshConverter)
        self._part = value

    @property
    def quadrature(self):
        return self._quadrature

    def _set_quadrature(self, value):
        assert isinstance(value, FiniteElementQuadrature)
        self._quadrature = value

    @property
    def boundary(self):
        return self._boundary

    def _set_boundary(self, value):
        assert isinstance(value, BoundaryCondition)
        self._boundary = value

    @property
    def material(self):
        return self._material

    @property
    def constraint_nodes(self):
        return self._constraint_nodes

    @property
    def free_nodes(self):
        return self._free_nodes

    @property
    def update_manager(self):
        if not isinstance(self._update_manager, UpdateManager):
            self._update_manager = UpdateManager()
        return self._update_manager

    def _set_material(self, material: Material):
        assert isinstance(material, Material)
        self._material: Material = material

    def _verify_mesh(self, mesh_order: str):
        assert mesh_order.lower() in ["linear", "quadratic"]
        assert self.part.mesh_order == mesh_order.lower()

    def _recover_mesh_property(self, order: str, space: str) -> Tuple[list, list]:
        assert space in ["lagrange", "nedelec"]
        assert order in ["linear", "quadratic"]
        order_int = 1 if order == "linear" else 2
        elements = self.part.recover_elements(space, order_int)
        points = self.part.recover_points(space, order_int)
        return elements, points

    def prepare_lagrange_integral_on_element(
        self, order: str
    ) -> Tuple[np.ndarray, np.ndarray, list, list]:
        assert order in ["linear", "quadratic"]
        funbasis_list = {
            "linear": basisfun_linear_triangle,
            "quadratic": basisfun_quadratic_triangle,
        }
        dersfunbasis_list = {
            "linear": dersbasisfun_linear_triangle,
            "quadratic": dersbasisfun_quadratic_triangle,
        }
        funbasis, dersfunbasis = self.quadrature.IntTriaQuad.evaluate_function(
            funbasis_list[order], dersfunbasis_list[order]
        )
        elements, points = self._recover_mesh_property(order, "lagrange")
        return funbasis, dersfunbasis, elements, points

    def prepare_lagrange_integral_on_boundary(
        self, order: str
    ) -> Tuple[np.ndarray, np.ndarray, list]:
        assert order in ["linear", "quadratic"]
        funbasis_list = {
            "linear": basisfun_linear_line,
            "quadratic": basisfun_quadratic_line,
        }
        dersfunbasis_list = {
            "linear": dersbasisfun_linear_line,
            "quadratic": dersbasisfun_quadratic_line,
        }
        funbasis, dersfunbasis = self.quadrature.LineQuad.evaluate_function(
            funbasis_list[order], dersfunbasis_list[order]
        )
        points = self._recover_mesh_property(order, "lagrange")[-1]
        return funbasis, dersfunbasis, points

    def verify_fun_args(self, quad_pts: np.ndarray, idx_el: int, args: dict):
        assert isinstance(
            args, dict
        ), "allowed extra arguments should be in dictionnary"
        shape_quadpts = np.shape(quad_pts)[1:]
        default = dict(quad_pts=quad_pts, idx_el=idx_el, shape_quadpts=shape_quadpts)
        args.update(default)

    def _assemble_scalar_mass(
        self, fun: Callable, order: str, **external_args
    ) -> sp.csr_array:
        assert order in ["linear", "quadratic"]
        (
            funbasis,
            dersfunbasis,
            elements,
            points,
        ) = self.prepare_lagrange_integral_on_element(order)
        row_ind, col_ind, val_ind = [], [], []
        for idx_el, idx_nodes_element in enumerate(elements):
            elem_coords = np.array([points[idx_nd] for idx_nd in idx_nodes_element])
            det_jac = self.evaluate_jacobien_of_field(dersfunbasis, elem_coords)[0]
            quadpts = self.evaluate_field(funbasis, elem_coords)
            self.verify_fun_args(quadpts, idx_el, external_args)
            prop = (
                fun(external_args) * self.quadrature.IntTriaQuad.quadweights * det_jac
            )
            Melem: np.ndarray = np.einsum("ki,k,kj->ij", funbasis, prop, funbasis)

            for i in range(Melem.shape[0]):
                for j in range(Melem.shape[1]):
                    row_ind.append(idx_nodes_element[i])
                    col_ind.append(idx_nodes_element[j])
                    val_ind.append(Melem[i, j])
        size_mat_tota = len(points)
        matrix = sp.coo_array(
            (val_ind, (row_ind, col_ind)), shape=(size_mat_tota, size_mat_tota)
        )
        matrix.eliminate_zeros()
        return sp.csr_array(matrix.tocsr())

    def _assemble_scalar_stiffness(
        self, fun: Callable, order: str, **external_args
    ) -> sp.csr_array:
        assert order in ["linear", "quadratic"]
        (
            funbasis,
            dersfunbasis,
            elements,
            points,
        ) = self.prepare_lagrange_integral_on_element(order)
        row_ind, col_ind, val_ind = [], [], []
        for idx_el, idx_nodes_element in enumerate(elements):
            elem_coords = np.array([points[idx_nd] for idx_nd in idx_nodes_element])
            det_jac, inv_jac = self.evaluate_jacobien_of_field(
                dersfunbasis, elem_coords
            )
            quadpts = self.evaluate_field(funbasis, elem_coords)
            self.verify_fun_args(quadpts, idx_el, external_args)
            prop = np.einsum(
                "ilk,lmk,jmk,k,k->ijk",
                inv_jac,
                fun(external_args),
                inv_jac,
                self.quadrature.IntTriaQuad.quadweights,
                det_jac,
                optimize=True,
            )
            Melem: np.ndarray = np.einsum(
                "kli,lmk,kmj->ij", dersfunbasis, prop, dersfunbasis
            )

            for i in range(Melem.shape[0]):
                for j in range(Melem.shape[1]):
                    row_ind.append(idx_nodes_element[i])
                    col_ind.append(idx_nodes_element[j])
                    val_ind.append(Melem[i, j])
        size_mat_tota = len(points)
        matrix = sp.coo_array(
            (val_ind, (row_ind, col_ind)), shape=(size_mat_tota, size_mat_tota)
        )
        matrix.eliminate_zeros()
        return sp.csr_array(matrix.tocsr())

    def _assemble_scalar_volume_force(self, fun: Callable, order: str) -> np.ndarray:
        assert order in ["linear", "quadratic"]
        (
            funbasis,
            dersfunbasis,
            elements,
            points,
        ) = self.prepare_lagrange_integral_on_element(order)
        data_array = []
        max_size_elem = 0
        for idx_nodes_element in elements:
            elem_coords = np.array([points[idx_nd] for idx_nd in idx_nodes_element])
            det_jac = self.evaluate_jacobien_of_field(dersfunbasis, elem_coords)[0]
            quadpts = self.evaluate_field(funbasis, elem_coords)
            Velem: np.ndarray = np.einsum(
                "ki,...k,k,k->...i",
                funbasis,
                fun(quadpts),
                self.quadrature.IntTriaQuad.quadweights,
                det_jac,
            )
            Velem = np.atleast_2d(Velem)
            nr, nc = Velem.shape
            assert nc == 3
            max_size_elem = max(max_size_elem, nr)
            for j in range(nc):
                data_array.append([idx_nodes_element[j], Velem[:, j]])
        array = np.zeros((max_size_elem, len(points)))
        for idx, val in data_array:
            array[:, idx] += val
        return array

    def _assemble_scalar_surface_force(
        self,
        fun: Callable,
        info: Callable,
        order: str,
    ) -> np.ndarray:
        assert order in ["linear", "quadratic"]
        (
            segment_funbasis,
            segment_dersfunbasis,
            points,
        ) = self.prepare_lagrange_integral_on_boundary(order=order)
        funspace = ("lagrange", 1) if order == "linear" else ("lagrange", 2)
        idx_nodes_segment_list = self.boundary.recognize_constraint(info, funspace)[1]
        data_array = []
        max_size_elem = 0
        for idx_nodes_segment in idx_nodes_segment_list:
            segment_coords = np.array([points[idx_nd] for idx_nd in idx_nodes_segment])
            jac = np.einsum("kl,lj->jk", segment_dersfunbasis, segment_coords)
            det_jac = np.linalg.norm(jac, axis=0)
            quadpts = self.evaluate_field(segment_funbasis, segment_coords)
            Velem: np.ndarray = np.einsum(
                "ki,...k,k,k->...i",
                segment_funbasis,
                fun(quadpts),
                self.quadrature.LineQuad.quadweights,
                det_jac,
            )
            Velem = np.atleast_2d(Velem)
            nr, nc = Velem.shape
            assert nc == 2 or nc == 3  # Linear or quadratic
            max_size_elem = max(max_size_elem, nr)
            for j in range(nc):
                data_array.append([idx_nodes_segment[j], Velem[:, j]])
        array = np.zeros((max_size_elem, len(points)))
        for idx, val in data_array:
            array[:, idx] += val
        return array

    ############################# PUBLIC FUNCTIONS #############################
    def set_update_manager(self, manager):
        if isinstance(manager, UpdateManager):
            self._update_manager = manager

    def get_size_of_arrays(self) -> int:
        return int(sum(self.boundary.nbofdofs.values()))

    def get_free_and_constraint_nodes(self) -> Tuple[list, list]:
        return self.free_nodes, self.constraint_nodes

    def clear_bcs(self, array_in: np.ndarray):
        assert isinstance(self.constraint_nodes, list)
        array_in[self.constraint_nodes] = 0.0

    def cut_array(self, input: np.ndarray) -> Dict[DOF, Sequence[int]]:
        assert input.ndim == 1
        assert len(input) >= self.get_size_of_arrays()
        output = {}
        for dof in self.boundary.dofs_index.keys():
            lim = self.boundary.limits_dofs[dof]
            output.update({dof: input[lim[0] : lim[1]]})
        return output

    def cut_matrix(self, input: np.ndarray) -> Dict[Tuple[DOF, DOF], Sequence[int]]:
        assert input.ndim == 2
        assert input.shape[0] >= self.get_size_of_arrays()
        assert input.shape[1] >= self.get_size_of_arrays()
        output = {}
        for dof_i in self.boundary.dofs_index.keys():
            for dof_j in self.boundary.dofs_index.keys():
                lim_i = self.boundary.limits_dofs[dof_i]
                lim_j = self.boundary.limits_dofs[dof_j]
                indices_i = np.arange(lim_i[0], lim_i[1])
                indices_j = np.arange(lim_j[0], lim_j[1])
                mat = input[np.ix_(indices_i, indices_j)]
                output.update({(dof_i, dof_j): mat})
        return output

    def evaluate_jacobien_of_field(
        self, dersfunbasis: np.ndarray, elem_coords: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        jac = np.einsum("kil,l...->...ik", dersfunbasis, elem_coords)
        return eval_inverse_and_determinant(jac)

    def evaluate_field(
        self, funbasis: np.ndarray, elem_coords: np.ndarray
    ) -> np.ndarray:
        return np.einsum("il,l...->...i", funbasis, elem_coords)

    def assemble_volumetric_force(self, fun_info: Dict[DOF, Callable]) -> np.ndarray:
        start = time()
        assert isinstance(fun_info, dict)
        assert all(
            callable(fun) for _, fun in fun_info.items()
        ), "Insert a list of functions"
        array_out = {}
        for dof, fun in fun_info.items():
            if dof is DOF.ALL:
                order = self.part.mesh_order
            else:
                order = (
                    "linear"
                    if self.boundary.dofs_index[dof][1][1] == 1
                    else "quadratic"
                )
            output = self._assemble_scalar_volume_force(fun, order)
            array_out.update({dof: output})
        logger.info(f"Computing volume force in {time() - start:.2e} seconds")
        return super().export_force(array_out)

    def assemble_surface_force(self, fun_info: Dict[DOF, Tuple[Callable, Callable]]):
        start = time()
        assert isinstance(fun_info, dict)
        assert all(
            isinstance(dof, DOF) and callable(cond) and callable(fun)
            for dof, (cond, fun) in fun_info.items()
        )
        array_out = {}
        for dof, (info, fun) in fun_info.items():
            if dof is DOF.ALL:
                order = self.part.mesh_order
            else:
                order = (
                    "linear"
                    if self.boundary.dofs_index[dof][1][1] == 1
                    else "quadratic"
                )
            output = self._assemble_scalar_surface_force(fun, info, order)
            array_out.update({dof: output})
        logger.info(f"Computing surface force in {time() - start:.2e} seconds")
        return super().export_force(array_out)

    def _smoothing_dual_field_at_cell(
        self, field_at_cell: np.ndarray, order: str
    ) -> np.ndarray:
        assert order in ["linear", "quadratic"]
        elements, points = self._recover_mesh_property(order, "lagrange")
        array_out = np.zeros(len(points))
        multiplicity = np.zeros_like(array_out)
        for e, elem in enumerate(elements):
            for idx in elem:
                array_out[idx] += field_at_cell[e]
                multiplicity[idx] += 1
        array_out /= multiplicity
        return array_out
