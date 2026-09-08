from .elastoplasticity import J2General
from typing import Tuple, Sequence
import numpy as np
import logging
from time import time

logger = logging.getLogger("SRC.MATERIAL")


class J2PlaneStress(J2General):
    """
    J2 Plasticity with plane stress constraint (σ₃₃ = σ₁₃ = σ₂₃ = 0).
    """

    def __init__(self, mat_args: dict):
        logger.info("NOTE: we apply plane stress constraints")

        # Initialize with 3D formulation
        super().__init__(mat_args=mat_args, is_unidimensional=False)
        self._set_TypeOfConstitutiveLaw("PLANE_STRESS")

    def set_linear_elastic_tensor(
        self, shape: Sequence[int], ndim: int = 3
    ) -> np.ndarray:
        """
        Create condensed elastic tensor for plane stress.

        For plane stress, the elastic relation between in-plane stresses and
        in-plane strains is:

        [σ₁₁]   E/(1-nu²) [1   nu   0      ] [ε₁₁]
        [σ₂₂] = ────────  [ni   1   0      ] [ε₂₂]
        [σ₁₂]             [0   0   (1-nu)/2] [γ₁₂]

        This is DIFFERENT from 3D because ε₃₃ is eliminated using σ₃₃ = 0.
        """
        E = self._linear_elasticity.elastic_modulus
        NU = self._linear_elasticity.poisson_ratio
        factor = E / (1 - NU**2)

        # Initialize 3x3x3x3 tensor (in-plane components only)
        tensor = np.zeros((ndim, ndim, ndim, ndim, *shape))

        # C₁₁₁₁ = C₂₂₂₂ = E/(1-ν²)
        tensor[0, 0, 0, 0, ...] = factor
        tensor[1, 1, 1, 1, ...] = factor

        # C₁₁₂₂ = C₂₂₁₁ = νE/(1-ν²)
        tensor[0, 0, 1, 1, ...] = NU * factor
        tensor[1, 1, 0, 0, ...] = NU * factor

        # C₁₂₁₂ = C₁₂₂₁ = C₂₁₁₂ = C₂₁₂₁ = G = E/(2(1+ν))
        G = E / (2 * (1 + NU))
        tensor[0, 1, 0, 1, ...] = G
        tensor[0, 1, 1, 0, ...] = G
        tensor[1, 0, 0, 1, ...] = G
        tensor[1, 0, 1, 0, ...] = G

        return tensor

    def eval_elastic_stress(self, strain: np.ndarray) -> np.ndarray:
        """
        Compute elastic stress using plane stress constitutive law.

        This uses the CONDENSED constitutive relation, not the 3D relation!
        """
        E = self._linear_elasticity.elastic_modulus
        NU = self._linear_elasticity.poisson_ratio
        factor = E / (1 - NU**2)

        stress = np.zeros_like(strain)

        # In-plane normal stresses
        stress[0, 0, ...] = factor * (strain[0, 0, ...] + NU * strain[1, 1, ...])
        stress[1, 1, ...] = factor * (strain[1, 1, ...] + NU * strain[0, 0, ...])

        # Shear stress
        G = E / (2 * (1 + NU))
        stress[0, 1, ...] = 2 * G * strain[0, 1, ...]
        stress[1, 0, ...] = stress[0, 1, ...]

        # Out-of-plane components are zero
        stress[2, 2, ...] = 0.0
        stress[0, 2, ...] = 0.0
        stress[2, 0, ...] = 0.0
        stress[1, 2, ...] = 0.0
        stress[2, 1, ...] = 0.0

        return stress

    def _elastic_predictor(
        self,
        strain_n1: np.ndarray,
        plasticstrain_n0: np.ndarray,
        back_n0: np.ndarray,
        plseq_n0: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Elastic predictor for plane stress.

        Key difference: use plane stress elastic law, not 3D elastic law.
        """
        # Compute trial strain
        strain_trial = strain_n1 - plasticstrain_n0

        # Use PLANE STRESS elastic law (not 3D!)
        stress_trial = self.eval_elastic_stress(strain_trial)

        # Compute von Mises in FULL 3D (including s₃₃ component!)
        # This is correct because s₃₃ = -(σ₁₁ + σ₂₂)/3 ≠ 0
        vm_trial = self._linear_elasticity.eval_von_mises_stress(
            stress_trial - np.sum(back_n0, axis=0)
        )

        # Check yield status
        J2_trial = vm_trial - self.isotropic_hardening.fun(plseq_n0)

        return stress_trial, J2_trial

    def _plastic_corrector(
        self,
        J2_trial: np.ndarray,
        stress_n0: np.ndarray,
        back_n0: np.ndarray,
        plseq_n0: np.ndarray,
        plasticstrain_n0: np.ndarray,
        update_tangent: bool = True,
    ):
        """
        Plastic corrector for plane stress.

        The return mapping works in FULL 3D space (because plasticity depends on
        deviatoric stress which has non-zero s₃₃ component), but we enforce
        σ₃₃ = 0 throughout.
        """
        # Initialize
        ndim = 3
        stress_shape = stress_n0.shape[2:]

        stress_n1 = np.copy(stress_n0)
        plseq_n1 = np.copy(plseq_n0)
        back_n1 = np.copy(back_n0)
        plasticstrain_n1 = np.copy(plasticstrain_n0)

        # Use plane stress elastic tangent as initial value
        consistent_tangent = None
        LIMIT = self._linear_elasticity.elastic_limit
        if np.any(J2_trial > self._threshold * LIMIT):

            start = time()

            # Select plastic points
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

            # Compute plastic multiplier and flow direction
            # NOTE: _prepare_parameters works in FULL 3D, which is correct!
            dgamma, hatback, normal, theta = self._prepare_parameters(
                stress_n1[idx_ten2d], back_n1[idx_ten3d], plseq_n1[idx_scalar]
            )

            # Update equivalent plastic strain
            plseq_n1[idx_scalar] += dgamma

            # Update stress (in full 3D)
            factor = 2.0 * self._linear_elasticity.lame_mu * np.sqrt(1.5)
            stress_n1[idx_ten2d] -= factor * dgamma * normal

            # CRITICAL: Re-enforce plane stress constraint
            stress_n1[2, 2, ...] = 0.0
            stress_n1[0, 2, ...] = 0.0
            stress_n1[2, 0, ...] = 0.0
            stress_n1[1, 2, ...] = 0.0
            stress_n1[2, 1, ...] = 0.0

            # Update plastic strain (including out-of-plane component!)
            factor = np.sqrt(1.5)
            plasticstrain_n1[idx_ten2d] += factor * dgamma * normal

            # Update back stress (in full 3D)
            self.kinematic_hardening.update_back_stress(
                idx_scalar, back_n1, normal, dgamma, is_unidimensional=False
            )

            # Update consistent tangent for plane stress
            if update_tangent:
                # Extract theta parameters (these are NOT rotations!)
                theta_1 = theta[0]  # Derivative correction: -3μ/(∂f/∂λ)
                theta_2 = theta[1]  # Geometric correction: 2μλ√(1.5)/||s||

                # Modified elastic moduli for plane stress with plasticity
                # Start with plane stress elastic moduli
                E = self._linear_elasticity.elastic_modulus
                NU = self._linear_elasticity.poisson_ratio

                # For plane stress, the effective Lamé parameters are:
                lame_lambda_ps = E * NU / (1 - NU**2)
                lame_mu_ps = E / (2 * (1 + NU))

                # Apply plastic corrections (similar to 3D but using plane stress moduli)
                # These corrections come from the consistency condition
                lame_lambda_corrected = (
                    lame_lambda_ps + 2.0 / 3.0 * lame_mu_ps * theta_2
                )
                lame_mu_corrected = lame_mu_ps - lame_mu_ps * theta_2

                # Plastic correction terms
                omega_1 = -2 * lame_mu_ps * (theta_1 - theta_2)
                omega_2 = -np.sqrt(2.0 / 3.0) * theta_1 * theta_2

                # Build consistent tangent using einsum (vectorized)
                idnt = np.eye(ndim)

                new_consistent_tangent = (
                    np.einsum(
                        "il,jm,...->ijlm...",
                        idnt,
                        idnt,
                        lame_lambda_corrected,
                        optimize=True,
                    )
                    + np.einsum(
                        "im,jl,...->ijlm...",
                        idnt,
                        idnt,
                        lame_mu_corrected,
                        optimize=True,
                    )
                    + np.einsum(
                        "ij,lm,...->ijlm...",
                        idnt,
                        idnt,
                        lame_mu_corrected,
                        optimize=True,
                    )
                    + np.einsum(
                        "il...,jm...,...->ijlm...",
                        normal,
                        normal,
                        omega_1,
                        optimize=True,
                    )
                    + np.einsum(
                        "il...,jm...,...->ijlm...",
                        hatback,
                        normal,
                        omega_2,
                        optimize=True,
                    )
                    - np.einsum(
                        "il...,jm...,...->ijlm...",
                        normal,
                        hatback,
                        omega_2,
                        optimize=True,
                    )
                )
                consistent_tangent = self.set_linear_elastic_tensor(stress_shape)
                consistent_tangent[idx_ten4d] = new_consistent_tangent

            logger.debug(f"Plane stress return mapping in {time() - start:.2e} seconds")

        return stress_n1, plseq_n1, back_n1, plasticstrain_n1, consistent_tangent

    def return_mapping(
        self, strain_n1: np.ndarray, plastic_vars: dict, update_tangent: bool = True
    ) -> Tuple[np.ndarray, dict]:
        """
        Return mapping for plane stress plasticity.

        Algorithm:
        1. Elastic predictor using PLANE STRESS elastic law
        2. Check yield in FULL 3D space (s₃₃ ≠ 0!)
        3. Plastic corrector in FULL 3D space
        4. Enforce σ₃₃ = 0 constraint
        5. Update tangent using plane stress formulation
        """
        self._validate_inputs(strain_n1, plastic_vars)

        # Get shapes
        ndim = strain_n1.shape[0]
        assert ndim == 3, "Plane stress requires 3D strain tensor"

        strain_shape = strain_n1.shape[2:]
        nb_chpar = self.kinematic_hardening.nb_chpar

        # Recover internal variables
        plasticstrain_n0 = plastic_vars.get("plastic_strain", np.zeros_like(strain_n1))
        plseq_n0 = plastic_vars.get("plastic_equivalent", np.zeros(strain_shape))
        back_n0 = plastic_vars.get(
            "back_stress", np.zeros((nb_chpar, ndim, ndim, *strain_shape))
        )

        # Elastic predictor (using plane stress elastic law)
        stress_trial, J2_trial = self._elastic_predictor(
            strain_n1, plasticstrain_n0, back_n0, plseq_n0
        )

        # Plastic corrector (works in full 3D but enforces σ₃₃ = 0)
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
