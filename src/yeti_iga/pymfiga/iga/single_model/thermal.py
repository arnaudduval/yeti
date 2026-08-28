from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from .cls.cspace import SingleSpatialModel
from typing import List, Tuple, Union, Optional
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MODEL")


class ThermalModel(SingleSpatialModel):
    def __init__(
        self,
        material: ThermalMaterial,
        patch: SinglePatch,
        boundary: BoundaryCondition,
    ):
        assert len(boundary.dofs_index) == 1
        super().__init__(material, patch, boundary)

        # Internal variables
        self._capacity_property: Optional[np.ndarray] = None
        self._conductivity_property: Optional[np.ndarray] = None
        self._scalar_mean_capacity: Optional[List[float]] = None
        self._scalar_mean_conductivity: Optional[List[np.ndarray]] = None

        # Compute a better approximation of scalar mean
        self.compute_capacity_property()
        self.compute_conductivity_property()

    @property
    def material(self):
        assert isinstance(self._material, ThermalMaterial)
        return self._material

    @property
    def scalar_mean_capacity(self):
        return self._scalar_mean_capacity

    @property
    def scalar_mean_conductivity(self):
        return self._scalar_mean_conductivity

    @property
    def capacity_property(self) -> np.ndarray:
        if self._capacity_property is None:
            return np.array([])
        return self._capacity_property

    @property
    def conductivity_property(self) -> np.ndarray:
        if self._conductivity_property is None:
            return np.array([])
        return self._conductivity_property

    def clear_properties(self):
        if self.material.hasnonlinearmass:
            self._capacity_property = None
            self._scalar_mean_capacity = None
        if self.material.hasnonlinearstiffness:
            self._conductivity_property = None
            self._scalar_mean_conductivity = None

    def compute_capacity_property(self, **mf_args):
        self.verify_fun_args(mf_args)
        if mf_args.get("temperature") is None:
            logger.debug("Temperature is necessary but is None")
            mf_args["temperature"] = np.zeros(self.part.nbqp_total)

        self._capacity_property = self.material.capacity(mf_args) * self.part.det_jac
        self._scalar_mean_capacity = [float(np.mean(self.capacity_property))]

    def compute_mf_mass(self, array_in: np.ndarray, **mf_args) -> np.ndarray:
        start = time()

        if self._capacity_property is None:
            self.compute_capacity_property(**mf_args)

        array_out = self.operator_engine.compute_mf_scalar_u_v(
            self.part.quadrule_list,
            self.capacity_property,
            array_in,
            allow_lumping=False,  # No lumping in heat transfer
            nurbs_weights=self.part.nurbs_weights,
        )
        logger.debug(f"Matrix free capacity in {time() - start:.2e} seconds")
        return array_out

    def compute_conductivity_property(self, **mf_args):
        self.verify_fun_args(mf_args)
        if mf_args.get("temperature") is None:
            logger.debug("Temperature is necessary but is None")
            mf_args["temperature"] = np.zeros(self.part.nbqp_total)

        self._conductivity_property = np.einsum(
            "ilk,lmk,jmk,k->ijk",
            self.part.inv_jac,
            self.material.conductivity(mf_args),
            self.part.inv_jac,
            self.part.det_jac,
            optimize=True,
        )
        self._scalar_mean_conductivity = [
            np.array(
                [np.mean(self.conductivity_property[i][i]) for i in range(self.ndim)]
            )
        ]

    def compute_mf_stiffness(self, array_in: np.ndarray, **mf_args) -> np.ndarray:
        start = time()

        if self._conductivity_property is None:
            self.compute_conductivity_property(**mf_args)

        array_out = self.operator_engine.compute_mf_scalar_gradu_gradv(
            self.part.quadrule_list,
            self.conductivity_property,
            array_in,
            nurbs_weights=self.part.nurbs_weights,
        )
        logger.debug(f"Matrix free conductivity in {time() - start:.2e} seconds")
        return array_out

    def interpolate_temperature(self, u_ctrlpts: np.ndarray) -> np.ndarray:
        # TODO: maybe use interpolate_field from ProjectionMixin ?
        u_interp = self.operator_engine.interpolate_meshgrid(
            self.part.quadrule_list,
            np.atleast_2d(u_ctrlpts),
            nurbs_weights=self.part.nurbs_weights,
        )
        return np.ravel(u_interp)

    def _assemble_internal_force(
        self,
        temp: np.ndarray,
        flux: np.ndarray,
        scalar_coefs: Union[tuple, list],
        **mf_args,
    ) -> np.ndarray:
        assert isinstance(scalar_coefs, (tuple, list))
        array_out = np.zeros(self.get_size_of_arrays())
        if scalar_coefs[0] != 0:
            array_out += scalar_coefs[0] * self.compute_mf_mass(flux, **mf_args)
        if scalar_coefs[1] != 0:
            array_out += scalar_coefs[1] * self.compute_mf_stiffness(temp, **mf_args)
        return array_out

    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        start = time()
        external_force: np.ndarray = kwargs["external_force"]
        scalar_coefs: tuple = kwargs["scalar_coefs"]
        arr_reshaped = np.reshape(array_in, (2, -1))
        mf_args = {"temperature": self.interpolate_temperature(arr_reshaped[0])}
        residual = external_force - self._assemble_internal_force(
            arr_reshaped[0], arr_reshaped[1], scalar_coefs, **mf_args
        )
        self.clear_bcs(residual)
        arr_flatten = np.hstack((residual, np.zeros_like(residual)))
        logger.debug(f"Computing residual in {time() - start:.2e} seconds")
        return arr_flatten, mf_args

    def compute_mf_tangent(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        scalar_coefs: tuple = kwargs["scalar_coefs"]
        arr_out = np.zeros(self.get_size_of_arrays())
        mf_args = {"temperature": kwargs.get("temperature")}
        if scalar_coefs[0] != 0:
            arr_out += scalar_coefs[0] * self.compute_mf_mass(array_in, **mf_args)
        if scalar_coefs[1] != 0:
            arr_out += scalar_coefs[1] * self.compute_mf_stiffness(array_in, **mf_args)
        return arr_out

    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        start = time()
        scalar_coefs: tuple = kwargs["scalar_coefs"]
        flux_factor: tuple = kwargs["flux_factor"]
        linear_solver: LinearSolver = kwargs["linear_solver_backend"]
        inner_tolerance: float = (
            kwargs.get("inner_tolerance") or linear_solver.config.tolerance
        )
        linear_solver.update(tolerance=inner_tolerance)
        linear_solver.set_cleandod(self.constraint_nodes)
        array_reshaped = np.reshape(array_in, (2, -1))

        if self.update_manager.should_update_preconditioner:
            self.preconditioner.add_scalar_space_time_correctors(
                mass_corrector=self.scalar_mean_capacity,
                stiffness_corrector=self.scalar_mean_conductivity,
            )

        self.preconditioner.update_space_eigenvalues(scalar_coefs=scalar_coefs)
        increment = linear_solver.solve(
            self.compute_mf_tangent,
            array_reshaped[0],
            Pfun=self.preconditioner.apply_spatial_preconditioner,
            **kwargs,
        )["sol"]

        if self.update_manager.should_update_material:
            self.clear_properties()

        logger.debug(f"Solving linearized system in {time() - start:.2e} seconds")
        return np.hstack((increment, increment * flux_factor))
