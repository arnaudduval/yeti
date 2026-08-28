from .core import Material
from typing import Union, Callable
import numpy as np
import logging

logger = logging.getLogger("SRC.MATERIAL")


class ThermalMaterial(Material):
    """
    In our work we consider nonlinar materials, i.e. its thermal properties
    could change depending on the position, temperature, etc.
    """

    _capacity = _conductivity = _ders_capacity = _ders_conductivity = None

    def __init__(self):
        super().__init__()
        logger.info("Setting THERMAL material")

        # Private variables
        self._has_uniform_capacity: bool = True
        self._has_uniform_conductivity: bool = True

    @property
    def capacity(self) -> Callable:
        if not callable(self._capacity):
            raise ValueError("Capacity has not been initialized")
        return self._capacity

    @property
    def conductivity(self) -> Callable:
        if not callable(self._conductivity):
            raise ValueError("Conductivity has not been initialized")
        return self._conductivity

    @property
    def ders_capacity(self) -> Callable:
        if not callable(self._ders_capacity):
            raise ValueError("Derivative of capacity has not been initialized")
        return self._ders_capacity

    @property
    def ders_conductivity(self) -> Callable:
        if not callable(self._ders_conductivity):
            raise ValueError("Derivative of conductivity has not been initialized")
        return self._ders_conductivity

    @property
    def hasnonlinearmass(self):
        return not self._has_uniform_capacity

    @property
    def hasnonlinearstiffness(self):
        return not self._has_uniform_conductivity

    def add_capacity(self, inpt: Union[Callable, float], is_uniform: bool):
        self._has_uniform_capacity = True if is_uniform else False
        self._capacity = self.set_scalar_property(inpt, is_uniform=is_uniform)

    def add_conductivity(
        self,
        inpt: Union[Callable, float, np.ndarray],
        is_uniform: bool,
        ndim: int,
    ):
        self._has_uniform_conductivity = True if is_uniform else False
        self._conductivity = self.set_tensor_property(
            inpt, ndim=ndim, is_uniform=is_uniform
        )

    def add_ders_capacity(self, inpt: Union[Callable, float], is_uniform: bool):
        self._ders_capacity = self.set_scalar_property(inpt, is_uniform=is_uniform)

    def add_ders_conductivity(
        self,
        inpt: Union[Callable, float, np.ndarray],
        is_uniform: bool,
        ndim: int,
    ):
        self._ders_conductivity = self.set_tensor_property(
            inpt, ndim=ndim, is_uniform=is_uniform
        )
