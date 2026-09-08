from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.common.material import ThermalMaterial
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from .cls.ctime import SingleSpaceTimeModel
from typing import List, Tuple, Optional
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MODEL")


class SpaceTimeThermalModel(SingleSpaceTimeModel):
    def __init__(
        self,
        material: ThermalMaterial,
        patch: SinglePatch,
        time_patch: SinglePatch,
        boundary: BoundaryCondition,
    ):
        assert isinstance(material, ThermalMaterial)
        assert len(boundary.dofs_index) == 1
        super().__init__(material, patch, time_patch, boundary)

        # Internal variables
        self._capacity_property: Optional[np.ndarray] = None
        self._conductivity_property: Optional[np.ndarray] = None
        self._ders_capacity_property: Optional[np.ndarray] = None
        self._ders_conductivity_property: Optional[np.ndarray] = None
        self._scalar_mean_capacity: Optional[List[float]] = None
        self._scalar_mean_conductivity: Optional[List[np.ndarray]] = None

        # Compute a better approximation of scalar mean
        self.compute_capacity_property()
        self.compute_conductivity_property()

    @property
    def material(self) -> ThermalMaterial:
        assert isinstance(self._material, ThermalMaterial)
        return self._material

    @property
    def scalar_mean_capacity(self):
        return self._scalar_mean_capacity

    @property
    def scalar_mean_conductivity(self):
        return self._scalar_mean_conductivity

    @property
    def capacity_property(self):
        if self._capacity_property is None:
            return np.array([])
        return self._capacity_property

    @property
    def ders_capacity_property(self):
        if self._ders_capacity_property is None:
            return np.array([])
        return self._ders_capacity_property

    @property
    def conductivity_property(self):
        if self._conductivity_property is None:
            return np.array([])
        return self._conductivity_property

    @property
    def ders_conductivity_property(self):
        if self._ders_conductivity_property is None:
            return np.array([])
        return self._ders_conductivity_property

    def clear_properties(self):
        if self.material.hasnonlinearmass:
            self._capacity_property = None
            self._scalar_mean_capacity = None
        if self.material.hasnonlinearstiffness:
            self._conductivity_property = None
            self._scalar_mean_conductivity = None

        # By context, the following are always nonlinear
        self._ders_capacity_property = None
        self._ders_conductivity_property = None

    def compute_capacity_property(self, **mf_args):
        self.verify_fun_args(mf_args)
        if mf_args.get("temperature") is None:
            logger.debug("Temperature is necessary but is None")
            mf_args["temperature"] = np.zeros(
                self.part.nbqp_total * self.time.nbqp_total
            )
        self._capacity_property = self.material.capacity(mf_args) * np.kron(
            np.ones_like(self.time.det_jac), self.part.det_jac
        )
        self._scalar_mean_capacity = [float(np.mean(self.capacity_property))]

    def compute_mf_sptm_mass(self, array_in: np.ndarray, **mf_args) -> np.ndarray:
        start = time()

        if self._capacity_property is None:
            self.compute_capacity_property(**mf_args)

        array_out = self.operator_engine.compute_mf_scalar_u_v(
            self.sptm_quadrature_list,
            self.capacity_property,
            array_in,
            allow_lumping=False,
            enable_spacetime=True,
            time_ders=(0, 1),
            nurbs_weights=self.sptm_nurbs_weights,
        )
        logger.debug(f"Matrix free capacity in {time() - start:.2e} seconds")
        return array_out

    def compute_ders_capacity_property(self, **mf_args):
        self.verify_fun_args(mf_args)
        if mf_args.get("temperature") is None:
            logger.debug("Temperature is necessary but is None")
            mf_args["temperature"] = np.zeros(
                self.part.nbqp_total * self.time.nbqp_total
            )
        # NOTE: Gradient is mandatory
        grad_temperature = mf_args["gradient"]
        self._ders_capacity_property = (
            self.material.ders_capacity(mf_args)
            * grad_temperature[-1]
            * np.kron(self.time.det_jac, self.part.det_jac)
        )

    def compute_mf_sptm_ders_mass(self, array_in: np.ndarray, **mf_args) -> np.ndarray:
        start = time()

        if self._ders_capacity_property is None:
            self.compute_ders_capacity_property(**mf_args)

        array_out = self.operator_engine.compute_mf_scalar_u_v(
            self.sptm_quadrature_list,
            self.ders_capacity_property,
            array_in,
            allow_lumping=False,
            enable_spacetime=True,
            time_ders=(0, 0),
            nurbs_weights=self.sptm_nurbs_weights,
        )
        logger.debug(f"Matrix free ders capacity in {time() - start:.2e} seconds")
        return array_out

    def compute_conductivity_property(self, **mf_args):
        self.verify_fun_args(mf_args)
        if mf_args.get("temperature") is None:
            logger.debug("Temperature is necessary but is None")
            mf_args["temperature"] = np.zeros(
                self.part.nbqp_total * self.time.nbqp_total
            )

        tmp1 = self.material.conductivity(mf_args) * np.kron(
            self.time.det_jac, self.part.det_jac
        )
        tmp1_reshaped = np.reshape(
            tmp1,
            (self.ndim, self.ndim, self.time.nbqp_total, -1),
        )
        tmp2 = np.einsum(
            "ilk,lmpk,jmk->ijpk",
            self.part.inv_jac,
            tmp1_reshaped,
            self.part.inv_jac,
            optimize=True,
        )
        self._conductivity_property = np.reshape(tmp2, (self.ndim, self.ndim, -1))
        self._scalar_mean_conductivity = [
            np.array(
                [np.mean(self.conductivity_property[i][i]) for i in range(self.ndim)]
            )
        ]

    def compute_mf_sptm_stiffness(self, array_in: np.ndarray, **mf_args) -> np.ndarray:
        start = time()

        if self._conductivity_property is None:
            self.compute_conductivity_property(**mf_args)

        array_out = self.operator_engine.compute_mf_scalar_gradu_gradv(
            self.sptm_quadrature_list,
            self.conductivity_property,
            array_in,
            enable_spacetime=True,
            time_ders=(0, 0),
            nurbs_weights=self.sptm_nurbs_weights,
        )
        logger.debug(f"Matrix free conductivity in {time() - start:.2e} seconds")
        return array_out

    def compute_ders_conductivity_property(self, **mf_args):
        if not callable(self.material.ders_conductivity):
            return
        self.verify_fun_args(mf_args)
        if mf_args.get("temperature") is None:
            logger.debug("Temperature is necessary but is None")
            mf_args["temperature"] = np.zeros(
                self.part.nbqp_total * self.time.nbqp_total
            )
        # NOTE: Gradient is mandatory
        grad_temperature = mf_args["gradient"]
        tmp1 = np.einsum(
            "ijk,jk,k->ik",
            self.material.ders_conductivity(mf_args),
            grad_temperature[:-1, :],
            np.kron(self.time.det_jac, self.part.det_jac),
            optimize=True,
        )
        tmp1_reshaped = np.reshape(tmp1, (self.ndim, self.time.nbqp_total, -1))
        tmp2 = np.einsum(
            "ilk,lpk->ipk", self.part.inv_jac, tmp1_reshaped, optimize=True
        )
        self._ders_conductivity_property = np.reshape(tmp2, (self.ndim, -1))

    def compute_mf_sptm_ders_stiffness(
        self, array_in: np.ndarray, **mf_args
    ) -> np.ndarray:
        start = time()

        if self._ders_conductivity_property is None:
            self.compute_ders_conductivity_property(**mf_args)

        array_out = self.operator_engine.compute_mf_scalar_gradu_v(
            self.sptm_quadrature_list,
            self.ders_conductivity_property,
            array_in,
            enable_spacetime=True,
            time_ders=(0, 0),
            nurbs_weights=self.sptm_nurbs_weights,
        )
        logger.debug(f"Matrix free ders-conductivity in {time() - start:.2e} seconds")
        return array_out

    def interpolate_sptm_temperature(
        self, u_ctrlpts: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        u_interp = np.ravel(
            self.operator_engine.interpolate_meshgrid(
                self.sptm_quadrature_list,
                np.atleast_2d(u_ctrlpts),
                nurbs_weights=self.sptm_nurbs_weights,
            ),
        )
        derstmp = self.operator_engine.eval_jacobien(
            self.sptm_quadrature_list,
            np.atleast_2d(u_ctrlpts),
            nurbs_weights=self.sptm_nurbs_weights,
        )[0]
        derstmp_reshaped = np.reshape(
            derstmp, (self.ndim + 1, self.time.nbqp_total, -1)
        )
        uders_interp = np.zeros_like(derstmp_reshaped)
        uders_interp[:-1, :, :] = np.einsum(
            "ipk,ijk->jpk",
            derstmp_reshaped[:-1, :, :],
            self.part.inv_jac,
            optimize=True,
        )
        uders_interp[-1] = np.einsum(
            "pk,p->pk",
            derstmp_reshaped[-1],
            np.ravel(self.time.inv_jac),
            optimize=True,
        )
        return u_interp, np.reshape(uders_interp, (self.ndim + 1, -1))

    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        start = time()
        external_force: np.ndarray = kwargs["external_force"]
        output = self.interpolate_sptm_temperature(array_in)
        mf_args = {"temperature": output[0], "gradient": output[1]}
        internal_force = self.compute_mf_sptm_mass(
            array_in, **mf_args
        ) + self.compute_mf_sptm_stiffness(array_in, **mf_args)
        residual = external_force - internal_force
        self.clear_bcs(residual)
        logger.debug(f"Computing residual in {time() - start:.2e} seconds")
        return residual, mf_args

    def _compute_mf_tangent(self, array_in, **kwargs):
        spacetime_type: str = kwargs["spacetime_type"]
        assert spacetime_type in ["picard", "newton"]
        mf_args = {"temperature": kwargs.get("temperature")}
        array_out = self.compute_mf_sptm_mass(
            array_in, **mf_args
        ) + self.compute_mf_sptm_stiffness(array_in, **mf_args)
        if spacetime_type == "newton":
            mf_args["gradient"] = kwargs.get("gradient")
            array_out += self.compute_mf_sptm_ders_mass(
                array_in, **mf_args
            ) + self.compute_mf_sptm_ders_stiffness(array_in, **mf_args)
        return array_out

    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        start = time()
        linear_solver: LinearSolver = kwargs["linear_solver_backend"]
        inner_tolerance: float = (
            kwargs.get("inner_tolerance") or linear_solver.config.tolerance
        )
        linear_solver.update(tolerance=inner_tolerance)
        linear_solver.set_cleandod(self.sptm_constraint_nodes)

        if self.update_manager.should_update_preconditioner:
            self.preconditioner.add_scalar_space_time_correctors(
                stiffness_corrector=self.scalar_mean_conductivity,
                advection_corrector=self.scalar_mean_capacity,
            )

        self.preconditioner.update_space_eigenvalues(scalar_coefs=(0, 1))
        sol = linear_solver.solve(
            self._compute_mf_tangent,
            array_in,
            Pfun=self.preconditioner.apply_spacetime_preconditioner,
            **kwargs,
        )["sol"]

        if self.update_manager.should_update_material:
            self.clear_properties()

        logger.debug(f"Solving linearized system in {time() - start:.2e} seconds")
        return sol
