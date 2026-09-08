from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.iga.boundary import BoundaryCondition
from yeti_iga.pymfiga.common.material.mechanical import IsotropicMat
from yeti_iga.pymfiga.common.numerics.operations import MatrixFree
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from .cls.cspace import SingleSpatialModel
from typing import List, Tuple, Optional
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MODEL")


class MechanicalModel(SingleSpatialModel):
    def __init__(
        self,
        material: IsotropicMat,
        patch: SinglePatch,
        boundary: BoundaryCondition,
    ):
        assert isinstance(material, IsotropicMat)
        assert len(boundary.dofs_index) == patch.ndim

        super().__init__(material, patch, boundary)

        # Internal variables
        self._mass_property: Optional[np.ndarray] = None
        self._stiffness_property: Optional[np.ndarray] = None
        self._scalar_mean_mass: Optional[List[float]] = None
        self._scalar_mean_stiffness: Optional[List[np.ndarray]] = None

        # Compute a better approximation of scalar mean
        self.compute_mass_property()
        self.compute_stiffness_property()

    @property
    def material(self):
        assert isinstance(self._material, IsotropicMat)
        return self._material

    @property
    def scalar_mean_mass(self):
        return self._scalar_mean_mass

    @property
    def scalar_mean_stiffness(self):
        return self._scalar_mean_stiffness

    @property
    def mass_property(self) -> np.ndarray:
        if self._mass_property is None:
            return np.array([])
        return self._mass_property

    @property
    def stiffness_property(self) -> np.ndarray:
        if self._stiffness_property is None:
            return np.array([])
        return self._stiffness_property

    def clear_properties(self):
        if self.material.hasnonlinearmass:
            self._mass_property = None
            self._scalar_mean_mass = None
        if self.material.hasnonlinearstiffness:
            self._stiffness_property = None
            self._scalar_mean_stiffness = None

    def compute_mass_property(self, **mf_args):
        self.verify_fun_args(mf_args)
        self._mass_property = self.material.density(mf_args) * self.part.det_jac
        self._scalar_mean_mass = [float(np.mean(self.mass_property))] * self.ndim

    def compute_mf_mass(self, array_in: np.ndarray, **mf_args) -> np.ndarray:
        start = time()

        if self._mass_property is None:
            self.compute_mass_property(**mf_args)

        array_in = np.reshape(array_in, (self.ndim, -1))
        array_out = np.zeros_like(array_in)
        for i in range(self.ndim):
            array_out[i] = self.operator_engine.compute_mf_scalar_u_v(
                self.part.quadrule_list,
                self.mass_property,
                array_in[i],
                allow_lumping=False,
                nurbs_weights=self.part.nurbs_weights,
            )
        logger.debug(f"Matrix free mass in {time() - start:.2e} seconds")
        return np.ravel(array_out)

    def compute_stiffness_property(self, **mf_args):

        self.verify_fun_args(mf_args)
        tangent = mf_args.get("consistent_tangent")
        if tangent is None:
            logger.debug("""
                Probably an error: either the material is elastic or forgot
                to transfer information to compute the tangent matrix
                """)
            tangent = self.material.set_linear_elastic_tensor((1,), ndim=self.ndim)

        self._stiffness_property = np.zeros(
            (
                self.ndim,
                self.ndim,
                self.ndim,
                self.ndim,
                self.part.nbqp_total,
            )
        )
        # NOTE: if elastic tangent, einsum will broadcast it automatically
        for i in range(self.ndim):
            for j in range(self.ndim):
                self._stiffness_property[i, j, ...] = np.einsum(
                    "ilk,lmk,jmk,k->ijk",
                    self.part.inv_jac,
                    tangent[i, j, : self.ndim, : self.ndim],
                    self.part.inv_jac,
                    self.part.det_jac,
                    optimize=True,
                )

        self._scalar_mean_stiffness = [
            np.array(
                [np.mean(self.stiffness_property[i][i][j][j]) for j in range(self.ndim)]
            )
            for i in range(self.ndim)
        ]

    def compute_mf_stiffness(self, array_in: np.ndarray, **mf_args) -> np.ndarray:
        "In this algorithm we only consider the linear elastic case"
        start = time()

        if self._stiffness_property is None:
            self.compute_stiffness_property(**mf_args)

        array_in = np.reshape(array_in, (self.ndim, -1))
        array_out = np.zeros_like(array_in)
        for i in range(self.ndim):
            array_out[i] = sum(
                self.operator_engine.compute_mf_scalar_gradu_gradv(
                    self.part.quadrule_list,
                    self.stiffness_property[i, j, ...],
                    array_in[j],
                    nurbs_weights=self.part.nurbs_weights,
                )
                for j in range(self.ndim)
            )
        logger.debug(f"Matrix free stiffness in {time() - start:.2e} seconds")
        return np.ravel(array_out)

    def interpolate_strain(
        self, array_in: np.ndarray, convert_to_3d: bool = False
    ) -> np.ndarray:
        "Compute strain field from displacement field"
        array_in = np.reshape(array_in, (self.ndim, -1))
        ders_par = self.operator_engine.eval_jacobien(
            self.part.quadrule_list, array_in, nurbs_weights=self.part.nurbs_weights
        )
        ders_phy = np.einsum("ijl,jkl->ikl", ders_par, self.part.inv_jac, optimize=True)
        grad_sym = 0.5 * (ders_phy + np.einsum("ijl->jil", ders_phy, optimize=True))
        size_tensor = 3 if convert_to_3d else self.ndim
        strain = np.zeros(
            (
                size_tensor,
                size_tensor,
                self.part.nbqp_total,
            )
        )
        strain[: self.ndim, : self.ndim, :] = grad_sym
        if (
            self.material.TypeOfConstitutiveLaw == "PLANE_STRESS"
            and self.part.ndim == 2
            and size_tensor == 3
        ):
            nu = self.material.poisson_ratio
            strain[2, 2, ...] = -nu / (1 - nu) * (strain[0, 0, ...] + strain[1, 1, ...])
        return strain

    def _assemble_internal_force(self, stress: np.ndarray) -> np.ndarray:
        prop = np.einsum(
            "ilk,ljk,k->ijk",
            self.part.inv_jac,
            stress[: self.ndim, : self.ndim, :],
            self.part.det_jac,
            optimize=True,
        )
        array_out = np.zeros(self.get_size_of_arrays()).reshape((self.ndim, -1))
        for i in range(self.ndim):
            for k in range(self.ndim):
                alpha_list = np.ones(self.ndim, dtype=int)
                alpha_list[k] = 3
                array_out[i] += MatrixFree.apply(
                    [
                        quadrule.weights[alpha]
                        for quadrule, alpha in zip(self.part.quadrule_list, alpha_list)
                    ],
                    prop[k, i, :],
                    is_transpose=False,
                )
        return np.ravel(array_out)

    def compute_residual(
        self, array_in: np.ndarray, **kwargs
    ) -> Tuple[np.ndarray, dict]:
        start = time()
        external_force: np.ndarray = kwargs["external_force"]
        plastic_vars: dict = kwargs.get("plastic_vars", {})
        convert_to_3d = False if self.ndim == 1 else True
        strain = self.interpolate_strain(array_in, convert_to_3d=convert_to_3d)
        stress, mf_args = self.material.return_mapping(strain, plastic_vars)
        internal_force = self._assemble_internal_force(stress)
        residual = external_force - internal_force
        self.clear_bcs(residual)
        logger.debug(f"Computing residual in {time() - start:.2e} seconds")
        return residual, internal_force, mf_args


    
    def solve_linearized_system(self, array_in: np.ndarray, **kwargs) -> np.ndarray:
        start = time()

        linear_solver: LinearSolver = kwargs["linear_solver_backend"]

        inner_tolerance: float = (
            kwargs.get("inner_tolerance") or linear_solver.config.tolerance
        )

        linear_solver.update(tolerance=inner_tolerance)
        linear_solver.set_cleandod(self.constraint_nodes)

        use_preconditioner = kwargs.get("use_preconditioner", True)
        preconditioner_type = kwargs.get("preconditioner_type", "fastdiag")

        if (not use_preconditioner) or (preconditioner_type is None):
            Pfun = None

        elif preconditioner_type == "fastdiag":

            if self.update_manager.should_update_preconditioner:
                self.preconditioner.add_scalar_space_time_correctors(
                    mass_corrector=self.scalar_mean_mass
                )

            self.preconditioner.update_space_eigenvalues(scalar_coefs=(1, 0))
            Pfun = self.preconditioner.apply_spatial_preconditioner

        elif preconditioner_type == "scaled_mass":

            Pfun = self.scaled_mass_preconditioner.apply_spatial_preconditioner

        else:
            raise ValueError(
                f"Unknown preconditioner_type: {preconditioner_type}"
            )

        sol = linear_solver.solve(
            self.compute_mf_mass,
            array_in,
            Pfun=Pfun,
            **kwargs,
        )["sol"]

        if self.update_manager.should_update_material:
            self.clear_properties()

        logger.debug(f"Solving linearized system in {time() - start:.2e} seconds")
        return sol   
    
