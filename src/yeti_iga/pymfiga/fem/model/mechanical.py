from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.common.material.mechanical import IsotropicMat
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition
from .core import SingleModel
from typing import Literal
from scipy import sparse as sp
import numpy as np


class MechanicalModel(SingleModel):
    def __init__(
        self,
        mesh,
        material: IsotropicMat,
        quadrature: FiniteElementQuadrature,
        boundary: BoundaryCondition,
        lagrange_order: Literal["linear", "quadratic"] = "linear",
    ):
        super().__init__(mesh, material, quadrature, boundary)
        self._verify_mesh(lagrange_order)

        # Internal variables
        self._lagrange_order = lagrange_order.lower()
        self._mass_matrix = None
        self._stiffness_matrix = None

    @property
    def material(self):
        assert isinstance(self._material, IsotropicMat)
        return self._material

    def clear_properties(self):
        if self.material.hasnonlinearmass:
            self._mass_matrix = None
        if self.material.hasnonlinearstiffness:
            self._stiffness_matrix = None

    def assemble_mass(self, **kwargs) -> sp.csr_array:
        if self._mass_matrix is not None:
            return self._mass_matrix
        mass_matrix = self._assemble_scalar_mass(
            self.material.density, order=self._lagrange_order, **kwargs
        )
        self._mass_matrix = sp.csr_array(
            sp.block_diag([mass_matrix] * self.part.ndim).tocsr()
        )
        self._mass_matrix.eliminate_zeros()
        return self._mass_matrix

    def compute_mf_mass(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        if self._mass_matrix is not None:
            return self._mass_matrix @ array_in
        self._mass_matrix = self.assemble_mass(**kwargs)
        return self._mass_matrix @ array_in

    def assemble_stiffness(self, **kwargs) -> sp.csr_array:
        def fun(args):
            idx_el = args["idx_el"]
            tensor = args["tensor"]
            return tensor[..., idx_el, :]

        if self._stiffness_matrix is not None:
            return self._stiffness_matrix

        ndim = self.part.ndim
        elements = self._recover_mesh_property(self._lagrange_order, "lagrange")[0]
        tangent = kwargs.get("consistent_tangent")
        if tangent is None:
            tangent = self.material.set_linear_elastic_tensor(
                (len(elements), 1), ndim=ndim
            )
        rowOfMatrices = []
        for i in range(ndim):
            columnOfMatrix = []
            for j in range(ndim):
                tensor = tangent[i, j, :ndim, :ndim, ...]
                columnOfMatrix.append(
                    self._assemble_scalar_stiffness(
                        fun,
                        order=self._lagrange_order,
                        tensor=tensor,
                    )
                )
            rowOfMatrices.append(columnOfMatrix)
        self._stiffness_matrix = sp.csr_array(sp.bmat(rowOfMatrices).tocsr())
        self._stiffness_matrix.eliminate_zeros()
        return self._stiffness_matrix

    def compute_mf_stiffness(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        if self._stiffness_matrix is not None:
            return self._stiffness_matrix @ array_in
        self._stiffness_matrix = self.assemble_stiffness(**kwargs)
        return self._stiffness_matrix @ array_in

    def _assemble_internal_force(self, stress: np.ndarray) -> np.ndarray:
        dersfunbasis, elements, points = self.prepare_lagrange_integral_on_element(
            self._lagrange_order
        )[1:]
        ndim = self.part.ndim
        size_mat_tota = len(points)
        array_out = np.zeros((ndim, size_mat_tota))
        for idx_el, idx_nodes_element in enumerate(elements):
            elem_coords = np.array([points[idx_nd] for idx_nd in idx_nodes_element])
            det_jac, inv_jac = self.evaluate_jacobien_of_field(
                dersfunbasis, elem_coords
            )
            stress_at_el = stress[:ndim, :ndim, idx_el, :]
            prop = np.einsum(
                "ilk,ljk,k,k->ijk",
                inv_jac,
                stress_at_el,
                self.quadrature.IntTriaQuad.quadweights,
                det_jac,
            )
            Velem: np.ndarray = np.einsum("kli,ljk->ij", dersfunbasis, prop)
            for i in range(Velem.shape[0]):
                array_out[:, idx_nodes_element[i]] += Velem[i]
        return np.ravel(array_out)

    def interpolate_strain(
        self, displacement: np.ndarray, convert_to_3d: bool = False
    ) -> np.ndarray:
        _, dersfunbasis, elements, points = self.prepare_lagrange_integral_on_element(
            self._lagrange_order
        )
        ndim = self.part.ndim
        displacement_2d = np.reshape(displacement, (ndim, -1))
        size_tensor = 3 if convert_to_3d else ndim
        strain = [np.array([]) for _ in range(len(elements))]
        for idx_el, idx_nodes_element in enumerate(elements):
            elem_coords = np.array([points[idx_nd] for idx_nd in idx_nodes_element])
            disp_elem = np.array(
                [displacement_2d[:, idx_nd] for idx_nd in idx_nodes_element]
            )
            inv_jac = self.evaluate_jacobien_of_field(dersfunbasis, elem_coords)[1]
            ders_par = np.einsum("kil,lj->jik", dersfunbasis, disp_elem)
            ders_phy = np.einsum("lik,jlk->ijk", inv_jac, ders_par)
            strain_el = 0.5 * (ders_phy + np.einsum("ijk->jik", ders_phy))
            newshape = (size_tensor, size_tensor, np.shape(strain_el)[-1])
            strain_el_expand = np.zeros(newshape)
            strain_el_expand[:ndim, :ndim, :] = strain_el[:ndim, :ndim, :]
            strain[idx_el] = strain_el_expand
        strain = np.moveaxis(np.asarray(strain), 0, 2)

        if (
            self.material.TypeOfConstitutiveLaw == "PLANE_STRESS"
            and self.part.ndim == 2
            and size_tensor == 3
        ):
            # ε₃₃ = -ν/(1-ν) * (ε₁₁ + ε₂₂)
            nu = self.material.poisson_ratio
            strain[2, 2, ...] = -nu / (1 - nu) * (strain[0, 0, ...] + strain[1, 1, ...])
        return strain

    def compute_residual(self, array_in: np.ndarray, **kwargs):
        external_force: np.ndarray = kwargs["external_force"]
        plastic_vars: dict = kwargs.get("plastic_vars", {})

        # Compute strain at each quadrature point
        strain = self.interpolate_strain(array_in, convert_to_3d=True)

        # Closest point projection in perfect plasticity
        stress, mechargs = self.material.return_mapping(strain, plastic_vars)

        # Compute internal force
        internal_force = self._assemble_internal_force(stress)

        # Compute residual
        residual = external_force - internal_force
        self.clear_bcs(residual)
        return residual,internal_force, mechargs

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

    def smoothing_dual_field_at_cell(self, field_at_cell: np.ndarray) -> np.ndarray:
        return super()._smoothing_dual_field_at_cell(
            field_at_cell, self._lagrange_order
        )
