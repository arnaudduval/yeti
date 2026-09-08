from yeti_iga.pymfiga.iga.boundary import BoundaryCondition
from yeti_iga.pymfiga.iga.geometry import SinglePatch
from yeti_iga.pymfiga.common.material import J2General
from yeti_iga.pymfiga.common.numerics.operations import MatrixFree
from yeti_iga.pymfiga.common.numerics.quadrature_rules import (
    IGAQuadratureRule,
    StandardGauss,
)
from yeti_iga.pymfiga.iga.single_model.cls.cspace import SingleSpatialModel
from typing import Callable, List, Tuple, Union
from scipy import sparse as sp
import numpy as np


class PipeModel(SingleSpatialModel):

    # Only for 1D straight beams
    _nbDoFs = 3

    def __init__(
        self,
        axis_patch: SinglePatch,
        circum_patch: SinglePatch,
        radial_quadrature: StandardGauss,
        boundary: BoundaryCondition,
        material: J2General,
    ):
        assert (
            axis_patch.ndim == 1 and circum_patch.ndim == 1
        ), "Patches can only be one dimensional"
        assert isinstance(radial_quadrature, IGAQuadratureRule), "Set a quadrature rule"
        super().__init__(material, axis_patch, boundary)
        assert (
            self.nbvars == self._nbDoFs
        ), f"Only for geometries with {self._nbDoFs} DoFs per node"

        self.circum_part: SinglePatch = circum_patch
        self.material: J2General = material
        self._pipe_quadrature_list: List[IGAQuadratureRule] = (
            self.part.quadrule_list[0]
            + self.circum_part.quadrule_list[0]
            + radial_quadrature
        )
        self.__sin_theta = None
        self.__cos_theta = None

    def _compute_strain_shell(
        self,
        a: float,
        U_0: np.ndarray,
        V_0: np.ndarray,
        u0: np.ndarray,
        v0: np.ndarray,
        w0: np.ndarray,
    ):
        """
        U0 and V0 are axis displacements (only depends on x)
        and u0, v0, w0 are mid-surface displacements (depends on x and theta).
        The first element of quadrule_list constaints the quadrature rule following x and
        the second element the quadrature rule following theta
        """
        axis_nbqp = self.part.nbqp_total
        circ_nbqp = self.circum_part.nbqp_total
        one_like_circ = np.ones(circ_nbqp)

        dU_0dx = self._pipe_quadrature_list[0].basis[1] @ U_0
        d2V_0dx2 = self._pipe_quadrature_list[0].basis[2] @ V_0

        pv0 = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[0],
            ],
            v0,
            is_transpose=False,
        )
        pw0 = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[0],
            ],
            w0,
            is_transpose=False,
        )
        du0dx = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[1],
                self._pipe_quadrature_list[1].basis[0],
            ],
            u0,
            is_transpose=False,
        )
        dv0dx = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[1],
                self._pipe_quadrature_list[1].basis[0],
            ],
            v0,
            is_transpose=False,
        )
        du0dth = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[1],
            ],
            u0,
            is_transpose=False,
        )
        dv0dth = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[1],
            ],
            v0,
            is_transpose=False,
        )

        eps = np.zeros((2, 2, axis_nbqp * circ_nbqp))
        eps[0, 0, :] = (
            np.kron(one_like_circ, dU_0dx)
            + du0dx
            - np.kron(self.__sin_theta, d2V_0dx2) * (a + pw0)
            - np.kron(self.__cos_theta, d2V_0dx2) * pv0
        )
        eps[1, 1, :] = 1 / a * (pw0 + dv0dth)
        eps[0, 1, :] = 1 / 2 * (1 / a * du0dth + dv0dx)
        eps[1, 0, :] = eps[0, 1, :]
        return eps

    def _compute_strain_z_linear(
        self,
        a: float,
        U_0: np.ndarray,
        V_0: np.ndarray,
        u0: np.ndarray,
        v0: np.ndarray,
        w0: np.ndarray,
    ):

        axis_nbqp = self.part.nbqp_total
        circ_nbqp = self.circum_part.nbqp_total
        one_like_circ = np.ones(circ_nbqp)

        d2V_0dx2 = self._pipe_quadrature_list[0].basis[2] @ V_0

        pv0 = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[0],
            ],
            v0,
            is_transpose=False,
        )
        pw0 = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[0],
            ],
            w0,
            is_transpose=False,
        )
        dv0dx = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[1],
                self._pipe_quadrature_list[1].basis[0],
            ],
            v0,
            is_transpose=False,
        )
        du0dth = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[1],
            ],
            u0,
            is_transpose=False,
        )
        dw0dth = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[1],
            ],
            w0,
            is_transpose=False,
        )
        d2w0dx2 = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[2],
                self._pipe_quadrature_list[1].basis[0],
            ],
            w0,
            is_transpose=False,
        )
        d2w0dth2 = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[0],
                self._pipe_quadrature_list[1].basis[2],
            ],
            w0,
            is_transpose=False,
        )
        d2w0dxdth = MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].basis[1],
                self._pipe_quadrature_list[1].basis[1],
            ],
            w0,
            is_transpose=False,
        )

        eps = np.zeros((2, 2, axis_nbqp * circ_nbqp))
        eps[0, 0, :] = (
            -np.kron(self.__sin_theta, d2V_0dx2)
            - np.kron(one_like_circ, d2w0dx2)
            + 1 / a * np.kron(self.__cos_theta, d2V_0dx2) * (dw0dth - pv0)
        )
        eps[1, 1, :] = -1 / a**2 * (d2w0dth2 + pw0)
        eps[0, 1, :] = -1 / (2 * a) * (2 * d2w0dxdth + 1 / a * du0dth - dv0dx)
        eps[1, 0, :] = eps[0, 1, :]
        return eps

    def _compute_mf_inertia(self, array_in, a: float):
        "Here u represents the acceleration at control points"
        axis_nbqp = self.part.nbqp_total
        circ_nbqp = self.circum_part.nbqp_total
        one_like_axis = np.ones(axis_nbqp)
        one_like_circ = np.ones(circ_nbqp)
        t_0, t_1, t_2, t_3, t_4, t_5 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

        t_0 += (
            self._pipe_quadrature_list[0].weights[0]
            @ ...
            @ (self._pipe_quadrature_list[0].basis[0] @ array_in[0])
        )
        #
        t_1 += (
            self._pipe_quadrature_list[0].weights[0]
            @ ...
            @ (self._pipe_quadrature_list[0].basis[0] @ array_in[1])
        )
        t_1 += -(
            self._pipe_quadrature_list[0].weights[0]
            @ ...
            @ (self._pipe_quadrature_list[0].basis[2] @ array_in[1])
        )
        t_1 += (
            self._pipe_quadrature_list[0].weights[0]
            @ ...
            @ (
                np.kron(self.__sin_theta, one_like_axis)
                * MatrixFree.apply(
                    [
                        self._pipe_quadrature_list[0].basis[0],
                        self._pipe_quadrature_list[1].basis[0],
                    ],
                    array_in[4],
                    is_transpose=False,
                )
                + np.kron(self.__cos_theta, one_like_axis)
                * MatrixFree.apply(
                    [
                        self._pipe_quadrature_list[0].basis[0],
                        self._pipe_quadrature_list[1].basis[0],
                    ],
                    array_in[3],
                    is_transpose=False,
                )
            )
        )
        #
        t_2 += MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[0],
            ],
            ...
            @ MatrixFree.apply(
                [
                    self._pipe_quadrature_list[0].basis[0],
                    self._pipe_quadrature_list[1].basis[0],
                ],
                array_in[2],
                is_transpose=False,
            ),
        )
        #
        t_3 += MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[0],
            ],
            ...
            @ MatrixFree.apply(
                [
                    self._pipe_quadrature_list[0].basis[0],
                    self._pipe_quadrature_list[1].basis[0],
                ],
                array_in[3],
                is_transpose=False,
            ),
        )
        t_3 += MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[0],
            ],
            ... @ (self._pipe_quadrature_list[0].basis[0] @ array_in[1]),
        )
        t_3 += -MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[0],
            ],
            ...
            @ (
                MatrixFree.apply(
                    [
                        self._pipe_quadrature_list[0].basis[0],
                        self._pipe_quadrature_list[1].basis[1],
                    ],
                    array_in[4],
                    is_transpose=False,
                )
                - MatrixFree.apply(
                    [
                        self._pipe_quadrature_list[0].basis[0],
                        self._pipe_quadrature_list[1].basis[0],
                    ],
                    array_in[3],
                    is_transpose=False,
                )
            ),
        )
        #
        t_4 += MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[0],
            ],
            ...
            @ MatrixFree.apply(
                [
                    self._pipe_quadrature_list[0].basis[0],
                    self._pipe_quadrature_list[1].basis[0],
                ],
                array_in[4],
                is_transpose=False,
            ),
        )
        t_4 += MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[0],
            ],
            ... @ (self._pipe_quadrature_list[0].basis[0] @ array_in[1]),
        )
        t_4 += -MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[0],
            ],
            ...
            @ MatrixFree.apply(
                [
                    self._pipe_quadrature_list[0].basis[2],
                    self._pipe_quadrature_list[1].basis[0],
                ],
                array_in[4],
                is_transpose=False,
            ),
        )
        t_4 += -MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[0],
            ],
            ...
            @ (
                MatrixFree.apply(
                    [
                        self._pipe_quadrature_list[0].basis[0],
                        self._pipe_quadrature_list[1].basis[2],
                    ],
                    array_in[4],
                    is_transpose=False,
                )
                - MatrixFree.apply(
                    [
                        self._pipe_quadrature_list[0].basis[0],
                        self._pipe_quadrature_list[1].basis[1],
                    ],
                    array_in[3],
                    is_transpose=False,
                )
            ),
        )
        return [t_0, t_1, t_2, t_3, t_4, t_5]

    def _assemble_internal_force(self, stress_m, stress_k, array_in):
        axis_nbqp = self.part.nbqp_total
        circ_nbqp = self.circum_part.nbqp_total
        one_like_axis = np.ones(axis_nbqp)
        one_like_circ = np.ones(circ_nbqp)
        t_0, t_1, t_2, t_3, t_4, t_5 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

        t_0 += MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[1],
                self._pipe_quadrature_list[1].weights[0],
            ],
            (... @ stress_m[0, 0, :]),
        )
        #
        t_1 += -MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[2],
                self._pipe_quadrature_list[1].weights[0],
            ],
            (... @ (stress_m[0, 0, :] + stress_k[0, 0, :])),
        )
        #
        t_2 += MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[1],
                self._pipe_quadrature_list[1].weights[0],
            ],
            (... @ stress_m[0, 0, :]),
        )
        t_2 += MatrixFree.apply(
            [
                self._pipe_quadrature_list[0].weights[0],
                self._pipe_quadrature_list[1].weights[1],
            ],
            (... @ stress_m[0, 1, :] - ... @ stress_k[0, 1, :]),
        )
        return
