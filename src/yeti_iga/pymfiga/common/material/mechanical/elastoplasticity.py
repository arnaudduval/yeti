from yeti_iga.pymfiga.common.numerics.solvers import NonLinearSolver
from .core import IsotropicMat, TensorOperations
from .linear_elasticity import LinearElasticity
from .hardening import IsotropicHardening, KinematicHardening
from typing import Tuple, Dict, List, Sequence
from abc import ABC, abstractmethod
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.MATERIAL")


class IsotropicPlasticity(IsotropicMat, ABC):
    def __init__(self, mat_args: dict, is_unidimensional):
        logger.info("Isotropic elastoplastic material")
        super().__init__("3DFULL", is_unidimensional)
        self._linear_elasticity = LinearElasticity(mat_args, is_unidimensional)
        self._initialize_hardening(mat_args)

    def _initialize_hardening(self, mat_args: dict):
        elas_lim = self._linear_elasticity.elastic_limit
        iso_args = mat_args.get("iso_hardening", {})
        self._isohard = IsotropicHardening(elas_lim, iso_args)

        kine_args = mat_args.get("kine_hardening", {})
        self._kinehard = KinematicHardening(kine_args)
        self._activated_plasticity = (
            False if (len(iso_args) == 0 and len(kine_args) == 0) else True
        )

    @property
    def elastic_modulus(self):
        return self._linear_elasticity.elastic_modulus

    @property
    def poisson_ratio(self):
        return self._linear_elasticity.poisson_ratio

    @property
    def isotropic_hardening(self):
        return self._isohard

    @property
    def kinematic_hardening(self):
        return self._kinehard

    @property
    def hasnonlinearstiffness(self):
        return self._activated_plasticity

    def eval_von_mises_stress(self, stress):
        return self._linear_elasticity.eval_von_mises_stress(stress)

    def eval_elastic_stress(self, strain: np.ndarray):
        return self._linear_elasticity.eval_elastic_stress(strain)

    def set_linear_elastic_tensor(self, shape: Sequence[int], ndim: int):
        return self._linear_elasticity.set_linear_elastic_tensor(shape, ndim)

    @abstractmethod
    def return_mapping(
        self, strain_n1: np.ndarray, plastic_vars: dict
    ) -> Tuple[np.ndarray, dict]:
        pass


class J2General(IsotropicPlasticity):

    _maxiters: int = 20
    _threshold: float = 1e-6

    def __init__(self, mat_args: dict, is_unidimensional: bool = False):
        logger.info("Setting J2 PLASTIC material")

        super().__init__(mat_args=mat_args, is_unidimensional=is_unidimensional)
        self._set_TypeOfConstitutiveLaw("3DFULL")

    def _validate_inputs(self, strain: np.ndarray, plastic_vars: dict):
        # Check dimensions
        if np.ndim(strain) <= 2:
            raise ValueError(
                f"Expected strain tensor with ndim > 2, got {np.ndim(strain)}"
            )
        ndim = strain.shape[0]
        if ndim not in [1, 3]:
            raise ValueError(f"Unsupported dimension: {ndim}. Must be 1 or 3.")

        # Check for invalid values
        if not np.all(np.isfinite(strain)):
            raise ValueError("Strain contains NaN or Inf values")

        # Validate plastic variables
        if plastic_vars:
            if "plastic_strain" in plastic_vars:
                ps: np.ndarray = plastic_vars["plastic_strain"]
                if ps.shape != strain.shape:
                    raise ValueError(
                        f"Plastic strain shape {ps.shape} doesn't match "
                        f"strain shape {strain.shape}"
                    )

    def _elastic_predictor(self, strain_n1, plasticstrain_n0, back_n0, plseq_n0):
        # Compute trial stress
        strain_trial = strain_n1 - plasticstrain_n0
        stress_trial = self.eval_elastic_stress(strain_trial)

        # Compute von mises stress
        vm_trial = self.eval_von_mises_stress(stress_trial - np.sum(back_n0, axis=0))
        if self.is_unidim:
            vm_trial = np.ravel(vm_trial)

        # Check yield status
        J2_trial = vm_trial - self.isotropic_hardening.fun(plseq_n0)
        return stress_trial, J2_trial

    def _prepare_parameters(
        self,
        stress_trial: np.ndarray,
        back_n0: np.ndarray,
        plseq_n0: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[np.ndarray]]:
        # Global variables
        YOUNG = self.elastic_modulus
        LAME_MU = self._linear_elasticity.lame_mu

        def compute_residual(dg: np.ndarray):
            output = self.kinematic_hardening.sum_chaboche_terms(dg, back_n0)
            shft_stress = stress_trial - output[0]
            if not self.is_unidim:
                shft_stress = TensorOperations.compute_deviatoric(shft_stress)
            shft_stress_norm = TensorOperations.compute_norm_tensor(shft_stress)
            normal_shft_stress = (
                np.sign(shft_stress)
                if self.is_unidim
                else shft_stress / shft_stress_norm
            )
            plseq_n1 = plseq_n0 + dg
            if self.is_unidim:
                yield_fun = -shft_stress_norm + (YOUNG + output[2]) * dg
            else:
                yield_fun = (
                    -np.sqrt(1.5) * shft_stress_norm
                    + 1.5 * (2 * LAME_MU + output[2]) * dg
                )
            yield_fun += self.isotropic_hardening.fun(plseq_n1)
            if self.is_unidim:
                ders_yield_fun = TensorOperations.compute_double_contraction(
                    normal_shft_stress, output[1]
                ) - (YOUNG + output[3])
            else:
                ders_yield_fun = np.sqrt(
                    1.5
                ) * TensorOperations.compute_double_contraction(
                    normal_shft_stress, output[1]
                ) - 1.5 * (
                    2 * LAME_MU + output[3]
                )
            ders_yield_fun -= self.isotropic_hardening.ders_fun(plseq_n1)
            extra_args: Dict[str, np.ndarray] = dict(
                hat_back=output[1],
                shifted_stress_norm=shft_stress_norm,
                normal_shifted_stress=normal_shft_stress,
                ders_yield_fun=ders_yield_fun,
            )
            return yield_fun, extra_args

        def solve_linearization(yield_fun: np.ndarray, **kwargs: dict):
            ders_yield_fun = kwargs.get("ders_yield_fun")
            return yield_fun / ders_yield_fun

        nonlinearsolver = NonLinearSolver(
            maxiters=self._maxiters,
            tolerance=self._threshold,
            allow_acceleration=False,
            allow_line_search=True,
            verbose=False,
        )
        dgamma = np.zeros_like(plseq_n0)
        output: dict = nonlinearsolver.solve(
            dgamma,
            compute_residual,
            solve_linearization,
        )
        extra_args: dict = output.get("extra_args", {})
        ders_yield_fun: np.ndarray = extra_args.get("ders_yield_fun", np.array([]))
        shifted_stress_norm: np.ndarray = extra_args.get(
            "shifted_stress_norm", np.array([])
        )
        normal_shifted_stress: np.ndarray = extra_args.get(
            "normal_shifted_stress", np.array([])
        )
        hat_back: np.ndarray = extra_args.get("hat_back", np.array([]))
        theta = [np.array([]), np.array([])]
        if self.is_unidim:
            if ders_yield_fun.size > 0:
                theta[0] = -YOUNG / ders_yield_fun
        else:
            if ders_yield_fun.size > 0:
                theta[0] = -3 * LAME_MU / ders_yield_fun
            if shifted_stress_norm.size > 0:
                theta[1] = 2 * LAME_MU * dgamma * np.sqrt(1.5) / shifted_stress_norm

        return dgamma, hat_back, normal_shifted_stress, theta

    def _plastic_corrector(
        self,
        J2_trial: np.ndarray,
        stress_n0: np.ndarray,
        back_n0: np.ndarray,
        plseq_n0: np.ndarray,
        plasticstrain_n0: np.ndarray,
        update_tangent: bool = True,
    ):

        # Global variables
        linear = self._linear_elasticity
        YOUNG = linear.elastic_modulus
        LIMIT = linear.elastic_limit
        MU = linear.lame_mu
        LAMBDA = linear.lame_lambda

        # Update
        ndim = stress_n0.shape[0]
        stress_shape = stress_n0.shape[2:]

        stress_n1 = np.copy(stress_n0)
        plseq_n1 = np.copy(plseq_n0)
        back_n1 = np.copy(back_n0)
        plasticstrain_n1 = np.copy(plasticstrain_n0)
        consistent_tangent = None

        if np.any(J2_trial > self._threshold * LIMIT):

            start = time()

            # Select the quadrature points
            idx_scalar = np.nonzero(J2_trial > self._threshold * LIMIT)
            idx_ten2d = (slice(None), slice(None), *idx_scalar)
            idx_ten3d = (slice(None), slice(None), slice(None), *idx_scalar)
            idx_ten4d = (
                slice(None),
                slice(None),
                slice(None),
                slice(None),
                *idx_scalar,
            )

            # Compute plastic-strain increment
            dgamma, hatback, normal, theta = self._prepare_parameters(
                stress_n1[idx_ten2d], back_n1[idx_ten3d], plseq_n1[idx_scalar]
            )

            # Update internal hardening variable
            plseq_n1[idx_scalar] += dgamma

            # Update stress
            factor = YOUNG if self.is_unidim else 2.0 * MU * np.sqrt(1.5)
            stress_n1[idx_ten2d] -= factor * dgamma * normal

            # Update plastic strain
            factor = 1.0 if self.is_unidim else np.sqrt(1.5)
            plasticstrain_n1[idx_ten2d] += factor * dgamma * normal

            # Update backstress
            self.kinematic_hardening.update_back_stress(
                idx_scalar, back_n1, normal, dgamma, is_unidimensional=self.is_unidim
            )

            # Update stiffness tensor
            if update_tangent:
                consistent_tangent = linear.set_linear_elastic_tensor(
                    stress_shape, ndim=ndim
                )
                if self.is_unidim:
                    new_consistent_tangent = YOUNG * (1 - theta[0])
                else:
                    idnt = np.eye(ndim)
                    omega_1 = -2 * MU * (theta[0] - theta[1])
                    omega_2 = -np.sqrt(2.0 / 3.0) * theta[0] * theta[1]
                    lame_lambda = LAMBDA + 2.0 / 3.0 * MU * theta[1]
                    lame_mu = MU - MU * theta[1]
                    new_consistent_tangent = (
                        np.einsum("il,jm,...->ijlm...", idnt, idnt, lame_lambda)
                        + np.einsum("im,jl,...->ijlm...", idnt, idnt, lame_mu)
                        + np.einsum("ij,lm,...->ijlm...", idnt, idnt, lame_mu)
                        + np.einsum("il...,jm...,...->ijlm...", normal, normal, omega_1)
                        + np.einsum(
                            "il...,jm...,...->ijlm...", hatback, normal, omega_2
                        )
                        - np.einsum(
                            "il...,jm...,...->ijlm...", normal, hatback, omega_2
                        )
                    )
                consistent_tangent[idx_ten4d] = new_consistent_tangent

            logger.debug(f"Return mapping computation in {time() - start:.2e} seconds")

        return stress_n1, plseq_n1, back_n1, plasticstrain_n1, consistent_tangent

    def return_mapping(
        self, strain_n1: np.ndarray, plastic_vars: dict, update_tangent: bool = True
    ) -> Tuple[np.ndarray, dict]:
        """Return mapping algorithm for multidimensional rate-independent plasticity."""

        self._validate_inputs(strain_n1, plastic_vars)

        # Get correct shapes
        ndim = strain_n1.shape[0]
        strain_shape = strain_n1.shape[2:]
        nb_chpar = self.kinematic_hardening.nb_chpar

        # Recover last values of internal variables
        plasticstrain_n0 = plastic_vars.get("plastic_strain", np.zeros_like(strain_n1))
        plseq_n0 = plastic_vars.get("plastic_equivalent", np.zeros(strain_shape))
        back_n0 = plastic_vars.get(
            "back_stress", np.zeros((nb_chpar, ndim, ndim, *strain_shape))
        )

        # Compute trial stress and trial J2 yield function
        stress_trial, J2_trial = self._elastic_predictor(
            strain_n1, plasticstrain_n0, back_n0, plseq_n0
        )

        # Plastic corrector
        output = self._plastic_corrector(
            J2_trial,
            stress_trial,
            back_n0,
            plseq_n0,
            plasticstrain_n0,
            update_tangent=update_tangent,
        )
        stress_n1, plseq_n1, back_n1, plasticstrain_n1, consistent_tangent = output

        new_plastic_vars = {
            "plastic_strain": plasticstrain_n1,
            "plastic_equivalent": plseq_n1,
            "back_stress": back_n1,
        }
        return stress_n1, {
            "consistent_tangent": consistent_tangent,
            "new_plastic_vars": new_plastic_vars,
        }
