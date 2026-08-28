from typing import Tuple, Callable, Sequence
import numpy as np
import logging

logger = logging.getLogger("SRC.MATERIAL")


class IsotropicHardening:
    _name = _fun = _ders_fun = None

    def __init__(self, elastic_limit: float, iso_args: dict):
        # Private variables
        self._elaslim: float = elastic_limit

        # Public variables
        self._name, (self._fun, self._ders_fun) = self._select_model(iso_args)
        logger.info(repr(self))

    @property
    def fun(self) -> Callable:
        if not callable(self._fun):
            raise ValueError("Function not defined")
        return self._fun

    @property
    def ders_fun(self) -> Callable:
        if not callable(self._ders_fun):
            raise ValueError("Derivative of function not defined")
        return self._ders_fun

    @property
    def name(self):
        if not isinstance(self._name, str):
            raise ValueError("Hardening not defined")
        return self._name

    def _select_model(self, iso_args: dict) -> Tuple[str, Tuple[Callable, Callable]]:
        models = {
            "none": self._model_none,
            "linear": self._model_linear,
            "swift": self._model_swift,
            "voce": self._model_voce,
        }
        model_name = str(iso_args.get("name", "none")).lower()
        if model_name not in models:
            raise ValueError("Unknown hardening model")
        return model_name, models[model_name](**iso_args)

    def _model_none(self, **args) -> Tuple[Callable, Callable]:
        factor = 1e8 * self._elaslim  # Numerical infinity
        return lambda a: factor * np.ones_like(a), lambda a: np.zeros_like(a)

    def _model_linear(self, **args) -> Tuple[Callable, Callable]:
        Eiso = args.get("Eiso")
        SY = self._elaslim
        return lambda a: SY + Eiso * a, lambda a: Eiso * np.ones_like(a)

    def _model_swift(self, **args) -> Tuple[Callable, Callable]:
        e0 = args.get("e0")
        n = args.get("n", 1.5)
        SY = self._elaslim
        return (
            lambda a: SY * (1 + a / e0) ** n,
            lambda a: SY * n / e0 * (1 + a / e0) ** (n - 1),
        )

    def _model_voce(self, **args) -> Tuple[Callable, Callable]:
        ssat = args.get("ssat")
        beta = args.get("beta", 1.0)
        SY = self._elaslim
        return (
            lambda a: SY + ssat * (1 - np.exp(-beta * a)),
            lambda a: ssat * beta * np.exp(-beta * a),
        )

    def __repr__(self) -> str:
        message = f"""
            Isotropic hardening model: {self.name}
            Using elastic limit: {self._elaslim}
        """
        return message


class KinematicHardening:
    _chaboche_table = None

    def __init__(self, kine_args: dict):
        # Private variable
        chaboche_table: np.ndarray = kine_args.get("parameters", np.array([[0, 0]]))
        self._set_chaboche_table(chaboche_table)
        logger.info(repr(self))

    @property
    def nb_chpar(self):
        return self.chaboche_table.shape[0]

    @property
    def chaboche_table(self):
        if not isinstance(self._chaboche_table, np.ndarray):
            raise ValueError("Chaboche table not defined")
        return self._chaboche_table

    def _set_chaboche_table(self, chaboche_table: np.ndarray):
        assert isinstance(chaboche_table, np.ndarray) and chaboche_table.ndim == 2
        self._chaboche_table = chaboche_table

    def sum_chaboche_terms(
        self, dgamma: np.ndarray, back: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        sum_back, hat_back = np.zeros_like(back[0, ...]), np.zeros_like(back[0, ...])
        const_1, const_2 = np.zeros_like(dgamma), np.zeros_like(dgamma)

        for i in range(self.nb_chpar):
            [c, d] = self.chaboche_table[i]
            term = 1 / (1 + d * dgamma)
            sum_back += back[i, ...] / term
            hat_back += d * back[i, ...] / term**2
            const_1 += c * term
            const_2 += c / term**2
        return sum_back, hat_back, const_1, const_2

    def update_back_stress(
        self,
        idx_scalar: Sequence[np.ndarray],
        back_n1: np.ndarray,
        normal: np.ndarray,
        dgamma: np.ndarray,
        is_unidimensional: bool = True,
    ):
        for i in range(self.nb_chpar):
            [c, d] = self.chaboche_table[i]
            idx = (slice(i), slice(None), slice(None), *idx_scalar)
            factor = c if is_unidimensional else c * np.sqrt(1.5)
            oldvalue = back_n1[idx] + factor * normal * dgamma
            back_n1[idx] = oldvalue / (1.0 + d * dgamma)

    def __repr__(self) -> str:
        message = f""""
            Initialize kinematic hardening model
            Number of parameters: {self.nb_chpar}
            Chaboche parameters: {self.chaboche_table}
        """
        return message
