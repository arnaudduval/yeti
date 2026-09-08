from yeti_iga.pymfiga.common.base.cls import BaseMaterial
from typing import Union, Callable
import numbers
from abc import ABC, abstractmethod
import numpy as np


class Material(BaseMaterial, ABC):

    _density = None

    def __init__(self):
        # Private variable:
        self._has_uniform_density = True

        # NOTE: by default we say that the material has a uniform density
        self.add_density(1.0, is_uniform=True)

    @property
    def density(self) -> Callable:
        if not callable(self._density):
            raise ValueError("Density has not been initialized")
        return self._density

    @property
    def hasnonlinearmass(self):
        return not self._has_uniform_density

    @property
    @abstractmethod
    def hasnonlinearstiffness(self) -> bool:
        pass

    def set_scalar_property(
        self, inpt: Union[Callable, float], is_uniform: bool = False
    ) -> Callable:
        if is_uniform:
            if isinstance(inpt, numbers.Real):
                func: Callable = lambda args: inpt * np.ones(
                    shape=args["shape_quadpts"]
                )
            elif callable(inpt):
                func: Callable = lambda args: inpt(args)
            else:
                raise NotImplementedError("Not implemented")
        else:
            assert callable(inpt), NotImplementedError("Not implemented")
            func: Callable = lambda args: inpt(args)
        return func

    def set_tensor_property(
        self,
        inpt: Union[Callable, float, np.ndarray],
        ndim: int = 2,
        is_uniform: bool = False,
    ) -> Callable:
        def broadcast(inpt: np.ndarray, shape_tensor: int, shape_quadpts: np.ndarray):
            tensor = np.zeros(shape=(shape_tensor, shape_tensor, *shape_quadpts))
            for i in range(shape_tensor):
                for j in range(shape_tensor):
                    tensor[i, j, ...] = inpt[i, j]
            return tensor

        if is_uniform:
            if isinstance(inpt, numbers.Real):
                newinpt: np.ndarray = np.asarray(inpt * np.eye(ndim))
                func: Callable = lambda args: broadcast(
                    newinpt, ndim, args["shape_quadpts"]
                )
            elif isinstance(inpt, np.ndarray):
                func: Callable = lambda args: broadcast(
                    inpt, ndim, args["shape_quadpts"]
                )
            elif callable(inpt):
                func: Callable = lambda args: inpt(args)
            else:
                raise NotImplementedError("Not implemented")
        else:
            assert callable(inpt), NotImplementedError("Not implemented")
            func: Callable = lambda args: inpt(args)
        return func

    def add_density(self, inpt: Union[Callable, float], is_uniform: bool):
        self._has_uniform_density = True if is_uniform else False
        self._density = self.set_scalar_property(inpt, is_uniform=is_uniform)
