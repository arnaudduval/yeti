from yeti_iga.pymfiga.fem.model.core import SingleModel
from typing import Callable, Optional, Literal
import numpy as np
import logging

logger = logging.getLogger("SRC.FEA.NORM")


class FeaNorm:
    def __init__(
        self,
        model: SingleModel,
        norm_type: Literal["l2", "h1", "semih1"],
        norm_args: dict,
    ):
        # Internal variables
        assert isinstance(model, SingleModel)
        assert norm_type in ["l2", "h1", "semih1"]
        self._model = model
        self._norm_type = norm_type
        self._norm_args = norm_args

    @property
    def model(self):
        return self._model

    @property
    def norm_type(self):
        return self._norm_type

    def _calculate_norm(
        self,
        u: np.ndarray,
        uders: np.ndarray,
        weights: np.ndarray,
    ) -> float:

        u2_l2, u2_sh1 = 0.0, 0.0
        if self.norm_type in ["l2", "h1"]:
            u2_l2 += np.sum(np.einsum("...k, k->...", u**2, weights))
        if self.norm_type in ["h1", "semih1"]:
            u2_sh1 += np.sum(np.einsum("...k, k->...", uders**2, weights))

        return np.sqrt(u2_l2 + u2_sh1)

    def evaluate(self):
        degree = self._norm_args.get("degree", "linear")
        uu: Optional[np.ndarray] = self._norm_args.get("u_at_nodes")
        assert isinstance(uu, np.ndarray), "The argument 'uu' must be a numpy"
        fun: Optional[Callable] = self._norm_args.get("fun")
        assert callable(fun), "The argument 'fun' must be a callable function."
        dersfun: Optional[Callable] = self._norm_args.get("dersfun")

        (
            funbasis,
            dersfunbasis,
            elements,
            points,
        ) = self.model.prepare_lagrange_integral_on_element(degree)
        eval_fun, eval_ders_fun = np.array([]), np.array([])
        eval_dif, eval_ders_dif = np.array([]), np.array([])
        abserror, norm_exact = 0.0, 0.0
        for _, idx_nodes_element in enumerate(elements):
            elem_coords = np.array([points[idx_nd] for idx_nd in idx_nodes_element])
            quadpts = self.model.evaluate_field(funbasis, elem_coords)
            det_jac, inv_jac = self.model.evaluate_jacobien_of_field(
                dersfunbasis, elem_coords
            )
            weights = self.model.quadrature.IntTriaQuad.quadweights * det_jac
            elem_uu = np.array([uu[idx_nd] for idx_nd in idx_nodes_element])

            if callable(fun):
                eval_fun = np.asarray(fun(quadpts))
                eval_uu = self.model.evaluate_field(funbasis, elem_uu)
                assert eval_fun.shape == eval_uu.shape
                eval_dif = eval_uu - eval_fun

            if callable(dersfun):
                eval_ders_fun = np.asarray(dersfun(quadpts))
                eval_ders_uu = np.einsum(
                    "lik,jlk->ijk",
                    inv_jac,
                    np.einsum("kil,lj->jik", dersfunbasis, elem_uu),
                )
                assert eval_ders_fun.shape == eval_ders_uu.shape
                eval_ders_dif = eval_ders_uu - eval_ders_fun

            abserror += self._calculate_norm(eval_dif, eval_ders_dif, weights)
            norm_exact += self._calculate_norm(eval_fun, eval_ders_fun, weights)

        relerror = abserror / norm_exact if norm_exact != 0 else abserror
        if norm_exact == 0:
            logger.warning("Warning: Dividing by zero")

        return abserror, relerror
