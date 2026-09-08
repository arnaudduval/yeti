from ..core import Material
from typing import Tuple, Sequence, Literal
from abc import ABC, abstractmethod
import numpy as np


class IsotropicMat(Material, ABC):
    _is_unidim = False
    _TypeOfConstitutiveLaw: Literal["3DFULL", "PLANE_STRESS"] = "3DFULL"

    def __init__(self, TypeOfConstitutiveLaw: str, is_unidim: bool):
        super().__init__()
        self._set_TypeOfConstitutiveLaw(TypeOfConstitutiveLaw)
        self._set_is_unidim(is_unidim)

    @property
    def is_unidim(self) -> bool:
        return self._is_unidim

    def _set_is_unidim(self, value):
        assert isinstance(value, bool)
        self._is_unidim = value

    @property
    def TypeOfConstitutiveLaw(self) -> Literal["3DFULL", "PLANE_STRESS"]:
        return self._TypeOfConstitutiveLaw

    def _set_TypeOfConstitutiveLaw(self, value):
        assert value in ["3DFULL", "PLANE_STRESS"]
        self._TypeOfConstitutiveLaw = value

    @property
    @abstractmethod
    def elastic_modulus(self) -> float:
        pass

    @property
    @abstractmethod
    def poisson_ratio(self) -> float:
        pass

    @abstractmethod
    def set_linear_elastic_tensor(self, shape: Sequence[int], ndim: int) -> np.ndarray:
        pass

    @abstractmethod
    def return_mapping(
        self, strain_n1: np.ndarray, plastic_vars: dict
    ) -> Tuple[np.ndarray, dict]:
        pass


class TensorOperations:
    @staticmethod
    def compute_double_contraction(
        tensor_1: np.ndarray, tensor_2: np.ndarray
    ) -> np.ndarray:
        return np.einsum("ij...,ij...->...", tensor_1, tensor_2, optimize=True)

    @staticmethod
    def compute_norm_tensor(tensor: np.ndarray) -> np.ndarray:
        return np.sqrt(TensorOperations.compute_double_contraction(tensor, tensor))

    @staticmethod
    def compute_trace(tensor: np.ndarray) -> np.ndarray:
        return np.einsum("ii...->...", tensor, optimize=True)

    @staticmethod
    def compute_deviatoric(tensor: np.ndarray) -> np.ndarray:
        trace = TensorOperations.compute_trace(tensor) / 3.0
        deviatoric = tensor.copy()
        for i in range(tensor.shape[0]):
            deviatoric[i, i, ...] -= trace
        return deviatoric
