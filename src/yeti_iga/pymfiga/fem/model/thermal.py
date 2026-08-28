from yeti_iga.pymfiga.common.numerics.quadrature_rules import FiniteElementQuadrature
from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.fem.boundary import BoundaryCondition
from .core import SingleModel
from typing import Tuple, Union, Literal
from scipy import sparse as sp
import numpy as np


class ThermalModel(SingleModel):
    def __init__(
        self,
        mesh,
        material: ThermalMaterial,
        quadrature: FiniteElementQuadrature,
        boundary: BoundaryCondition,
        lagrange_order: Literal["linear", "quadratic"] = "linear",
    ):
        super().__init__(mesh, material, quadrature, boundary)
        self._verify_mesh(lagrange_order)

        # Internal variables
        self._lagrange_order: str = lagrange_order.lower()
        self._mass_matrix = None
        self._stiffness_matrix = None

    @property
    def material(self):
        assert isinstance(self._material, ThermalMaterial)
        return self._material

    def clear_properties(self):
        if self.material.hasnonlinearmass:
            self._mass_matrix = None
        if self.material.hasnonlinearstiffness:
            self._stiffness_matrix = None

    def assemble_mass(self, **kwargs):
        if self._mass_matrix is not None:
            return self._mass_matrix
        self._mass_matrix = self._assemble_scalar_mass(
            self.material.capacity, order=self._lagrange_order, **kwargs
        )
        return self._mass_matrix

    def compute_mf_mass(self, array_in, **kwargs):
        if self._mass_matrix is not None:
            return self._mass_matrix @ array_in
        self._mass_matrix = self.assemble_mass(**kwargs)
        return self._mass_matrix @ array_in

    def assemble_stiffness(self, **kwargs):
        if self._stiffness_matrix is not None:
            return self._stiffness_matrix
        self._stiffness_matrix = self._assemble_scalar_stiffness(
            self.material.conductivity, order=self._lagrange_order, **kwargs
        )
        return self._stiffness_matrix

    def compute_mf_stiffness(self, array_in, **kwargs):
        if self._stiffness_matrix is not None:
            return self._stiffness_matrix @ array_in
        self._stiffness_matrix = self.assemble_stiffness(**kwargs)
        return self._stiffness_matrix @ array_in

    def interpolate_temperature(self, temperature: np.ndarray) -> np.ndarray:
        funbasis, _, elements, _ = self.prepare_lagrange_integral_on_element(
            self._lagrange_order
        )
        temp_interp = [np.array([]) for _ in range(len(elements))]
        for idx_el, idx_nodes_element in enumerate(elements):
            temp_elem = np.array([temperature[idx_nd] for idx_nd in idx_nodes_element])
            temp_interp[idx_el] = self.evaluate_field(funbasis, temp_elem)
        return np.asarray(temp_interp)

    def _assemble_internal_force(
        self,
        temp: np.ndarray,
        flux: np.ndarray,
        scalar_coefs: Union[tuple, list],
        **kwargs
    ) -> np.ndarray:
        assert isinstance(scalar_coefs, (tuple, list))
        assert len(scalar_coefs) > 1
        array_out = np.zeros(self.get_size_of_arrays())
        if scalar_coefs[0] != 0:
            mass = self.assemble_mass(**kwargs)
            array_out += scalar_coefs[0] * mass @ flux
        if scalar_coefs[1] != 0:
            stiff = self.assemble_stiffness(**kwargs)
            array_out += scalar_coefs[1] * stiff @ temp
        return array_out

    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        external_force: np.ndarray = kwargs["external_force"]
        scalar_coefs: tuple = kwargs["scalar_coefs"]
        arr_reshaped = np.reshape(array_in, (2, -1))
        mf_args = {"temperature": self.interpolate_temperature(arr_reshaped[0])}
        internal_force = self._assemble_internal_force(
            arr_reshaped[0], arr_reshaped[1], scalar_coefs, **mf_args
        )
        residual = external_force - internal_force
        self.clear_bcs(residual)
        arr_flatten = np.hstack((residual, np.zeros_like(residual)))
        return arr_flatten, mf_args

    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        scalar_coefs: tuple = kwargs["scalar_coefs"]
        flux_factor: tuple = kwargs["flux_factor"]
        array_reshaped = np.reshape(array_in, (2, -1))

        mf_args = {"temperature": kwargs.get("temperature")}
        nrows = self.get_size_of_arrays()
        tangent_matrix = sp.csr_array((nrows, nrows))
        if scalar_coefs[0] != 0:
            mass = self.assemble_mass(**mf_args)
            tangent_matrix += scalar_coefs[0] * mass
        if scalar_coefs[1] != 0:
            stiff = self.assemble_stiffness(**mf_args)
            tangent_matrix += scalar_coefs[1] * stiff

        free_nodes = self.get_free_and_constraint_nodes()[0]
        increment = np.zeros_like(array_reshaped[0])
        increment[free_nodes] = LinearSolver.direct(
            tangent_matrix[np.ix_(free_nodes, free_nodes)],
            array_reshaped[0][free_nodes],
        )["sol"]

        if self.update_manager.should_update_material:
            self.clear_properties()

        return np.hstack((increment, increment * flux_factor))

    def smoothing_dual_field_at_cell(self, field_at_cell: np.ndarray) -> np.ndarray:
        return super()._smoothing_dual_field_at_cell(
            field_at_cell, self._lagrange_order
        )
