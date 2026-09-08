from typing import Tuple, Callable, Optional
from .data import legendre_table
from abc import ABC, abstractmethod
from typing import Tuple
import numpy as np


class FiniteElementQuadrature:
    class QuadratureRule(ABC):
        def __init__(self):
            self._quadpoints, self._quadweights = np.array([]), np.array([])

        @property
        def quadpoints(self):
            return self._quadpoints

        @quadpoints.setter
        def quadpoints(self, value: np.ndarray):
            assert isinstance(value, np.ndarray)
            self._quadpoints = value

        @property
        def quadweights(self):
            return self._quadweights

        @quadweights.setter
        def quadweights(self, value: np.ndarray):
            assert isinstance(value, np.ndarray)
            self._quadweights = value

        def evaluate_function(
            self, fun: Optional[Callable], dersfun: Optional[Callable], **kwargs
        ) -> Tuple[np.ndarray, np.ndarray]:
            evalbasisfun, evaldersbasisfun = [], []
            if callable(fun):
                evalbasisfun = [fun(pts) for pts in self.quadpoints]
            if callable(dersfun):
                evaldersbasisfun = [dersfun(pts) for pts in self.quadpoints]
            return np.array(evalbasisfun), np.array(evaldersbasisfun)

        @abstractmethod
        def _set_quadpos_weights(self, order: int) -> None:
            pass

    class ParametricLine(QuadratureRule):
        def __init__(self, order: int):
            super().__init__()
            # NOTE: quadrature pos and weights: (k,)
            self._set_quadpos_weights(order)

        def _set_quadpos_weights(self, order: int):
            self.quadpoints, self.quadweights = legendre_table(order)

    class ParametricIntTriangle(QuadratureRule):
        def __init__(self, order: int):
            super().__init__()
            # NOTE: quadrature pos and weights: (k,)
            self._set_quadpos_weights(order)

        def _set_quadpos_weights(self, order: int):
            if order == 1:
                self.quadpoints = np.array([[1 / 3, 1 / 3]])
                self.quadweights = np.array([1 / 2])
            else:
                points, weights = [], []
                for pts_i, wtg_i in zip(*legendre_table(order)):
                    for pts_j, wtg_j in zip(*legendre_table(order)):
                        pts_l = (1 + pts_i) / 2
                        pts_r = (1 - pts_i) * (1 + pts_j) / 4
                        wgt = (1 - pts_i) * wtg_i * wtg_j / 8
                        points.append([pts_l, pts_r])
                        weights.append(wgt)
                self.quadpoints = np.array(points)
                self.quadweights = np.array(weights)

    class ParametricExtTriangle(QuadratureRule):

        # Points:
        # 0: (0, 0), 1: (1, 0), 2: (0, 1)
        _COORDS = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        # Lines:
        # 0-1: 0<xi<1, nu=0, 1-2: 0<xi<1, nu=1-xi, 2-0: 0<nu<1, xi=0
        _CONNECTIVITY = [[0, 1], [1, 2], [2, 0]]

        def __init__(self, order: int):
            super().__init__()
            # NOTE: quadrature pos and weights: (e, k,)
            self._set_quadpos_weights(order)

        def _set_quadpos_weights(self, order: int):
            points, weights = [], []
            for idxleft, idxright in self._CONNECTIVITY:
                pleft = self._COORDS[idxleft]
                pright = self._COORDS[idxright]
                points_edge, weights_edge = [], []
                for pts, wgt in zip(*legendre_table(order)):
                    points_edge.append(0.5 * ((pright - pleft) * pts + pleft + pright))
                    weights_edge.append(0.5 * np.linalg.norm(pright - pleft) * wgt)
                points.append(points_edge)
                weights.append(weights_edge)
            self.quadpoints = np.array(points)
            self.quadweights = np.array(weights)

        def evaluate_function(
            self,
            fun: Optional[Callable],
            dersfun: Optional[Callable],
            edge_local_idx: int = 0,
            **kwargs
        ) -> Tuple[np.ndarray, np.ndarray]:
            evalbasisfun, evaldersbasisfun = None, None
            if callable(fun):
                evalbasisfun = [fun(pts) for pts in self.quadpoints[edge_local_idx]]
            if callable(dersfun):
                evaldersbasisfun = [
                    dersfun(pts) for pts in self.quadpoints[edge_local_idx]
                ]
            return np.array(evalbasisfun), np.array(evaldersbasisfun)

    def __init__(self, interior_order: int = 2, boundary_order: int = 1):
        self._LineQuad = self.ParametricLine(boundary_order)
        self._IntTriaQuad = self.ParametricIntTriangle(interior_order)
        self._ExtTriaQuad = self.ParametricExtTriangle(boundary_order)

    @property
    def IntTriaQuad(self):
        return self._IntTriaQuad

    @property
    def LineQuad(self):
        return self._LineQuad

    @property
    def ExtTriaQuad(self):
        return self._ExtTriaQuad
