from yeti_iga.pymfiga.common.base.enum import Constants
from .core import IsotropicMat, TensorOperations
from typing import Tuple, Sequence
import numpy as np
import logging

logger = logging.getLogger("SRC.MATERIAL")


class LinearElasticity(IsotropicMat):
    _timoshenko_ratio = _elastic_modulus = _poisson_ratio = None
    _elastic_limit = _lame_lambda = _lame_mu = None

    def __init__(self, mat_args: dict, is_unidimensional: bool = False):
        super().__init__("3DFULL", is_unidimensional)
        # Public variables
        self._initialize_properties(mat_args)
        logger.info(repr(self))

    @property
    def timoshenko_ratio(self) -> float:
        if self._timoshenko_ratio is None:
            raise ValueError("Timoshenko ratio not initialized")
        return self._timoshenko_ratio

    @property
    def elastic_modulus(self) -> float:
        if self._elastic_modulus is None:
            raise ValueError("Elastic modulus not initialized")
        return self._elastic_modulus

    @property
    def poisson_ratio(self) -> float:
        if self._poisson_ratio is None:
            raise ValueError("Poisson ratio not initialized")
        return self._poisson_ratio

    @property
    def elastic_limit(self) -> float:
        if self._elastic_limit is None:
            raise ValueError("Elastic limit not initialized")
        return self._elastic_limit

    @property
    def lame_lambda(self) -> float:
        if self._lame_lambda is None:
            lam = (
                self.poisson_ratio
                * self.elastic_modulus
                / ((1 + self.poisson_ratio) * (1 - 2 * self.poisson_ratio))
            )
            self._lame_lambda = lam
        return self._lame_lambda

    @property
    def lame_mu(self) -> float:
        if self._lame_mu is None:
            mu = self.elastic_modulus / (2 * (1 + self.poisson_ratio))
            self._lame_mu = mu
        return self._lame_mu

    @property
    def hasnonlinearstiffness(self):
        # By definition, linear elasticity is not nonlinear
        return False

    def _initialize_properties(self, mat_args: dict):
        self._timoshenko_ratio = mat_args.get("timoshenko_ratio", 5.0 / 6.0)
        self._elastic_modulus = mat_args.get("elastic_modulus", 0.0)
        self._poisson_ratio = mat_args.get("poisson_ratio", 0.0)
        self._elastic_limit = mat_args.get("elastic_limit", Constants.INFTY)

    @staticmethod
    def set_lame_tensor(coef_1: float, coef_2: float, ndim: int):
        idnt = np.eye(ndim)
        tensor = coef_1 * np.einsum("ik,jl->ijkl", idnt, idnt) + coef_2 * (
            np.einsum("il,jk->ijkl", idnt, idnt) + np.einsum("ij,kl->ijkl", idnt, idnt)
        )
        return tensor

    def set_linear_elastic_tensor(self, shape: Sequence[int], ndim: int):
        tensor = self.elastic_modulus * np.ones((ndim, ndim, ndim, ndim, *shape))
        if not self.is_unidim:
            lame_tensor = LinearElasticity.set_lame_tensor(
                self.lame_lambda, self.lame_mu, ndim
            )
            indices = np.indices((ndim, ndim, ndim, ndim), dtype=int)
            for i, j, l, m in zip(*[arr.flatten() for arr in indices]):
                tensor[i, j, l, m, ...] = lame_tensor[i, j, l, m]
        return tensor

    def eval_elastic_stress(self, strain: np.ndarray) -> np.ndarray:
        if self.is_unidim:
            stress = self.elastic_modulus * strain
        else:
            trace_strain = TensorOperations.compute_trace(strain)
            stress = 2 * self.lame_mu * strain
            for i in range(strain.shape[0]):
                stress[i, i, ...] += self.lame_lambda * trace_strain
        return stress

    def eval_von_mises_stress(self, stress: np.ndarray) -> np.ndarray:
        if self.is_unidim:
            return TensorOperations.compute_norm_tensor(stress)
        stress_dev = TensorOperations.compute_deviatoric(stress)
        return np.sqrt(1.5) * TensorOperations.compute_norm_tensor(stress_dev)

    def return_mapping(
        self, strain_n1: np.ndarray, plastic_vars: dict
    ) -> Tuple[np.ndarray, dict]:
        stress = self.eval_elastic_stress(strain_n1)
        return stress, {}

    def __repr__(self) -> str:
        message = f""""
            Linear elastic material with properties:
            - Young modulus: {self.elastic_modulus:.2e}
            - Poisson ratio: {self.poisson_ratio}
        """
        return message
