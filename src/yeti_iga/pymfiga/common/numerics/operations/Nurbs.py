from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from .core import Operations
from .Bsplines import BsplineOperations
from .matrix_free import MatrixFree
from typing import List, Tuple, Sequence, Literal
from scipy import sparse as sp
import numpy as np


class NurbsOperations(Operations):
    @staticmethod
    def _project_nurbs_weights(
        quadrule_list: List[IGAQuadratureRule],
        nurbs_weights: np.ndarray,
        nders: int = 0,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        nurbs_weights = np.atleast_2d(nurbs_weights)
        weights_proj = np.zeros([])
        jac_weights_proj = np.zeros([])
        hess_weights_proj = np.zeros([])

        if nders >= 0:
            weights_proj = BsplineOperations.interpolate_meshgrid(
                quadrule_list, nurbs_weights
            )[0]
        if nders >= 1:
            jac_weights_proj = BsplineOperations.eval_jacobien(
                quadrule_list, nurbs_weights
            )[0]
        if nders >= 2:
            hess_weights_proj = BsplineOperations.eval_hessian(
                quadrule_list, nurbs_weights
            )[0]

        if nders == 0:
            return 1 / weights_proj, np.zeros([]), np.zeros([])
        elif nders == 1:
            return 1 / weights_proj, jac_weights_proj / weights_proj**2, np.zeros([])
        elif nders == 2:
            return (
                1 / weights_proj,
                jac_weights_proj / weights_proj**2,
                hess_weights_proj / weights_proj**3,
            )
        else:
            raise NotImplementedError()

    @staticmethod
    def _apply_scaling(array_in, weights, depth=1) -> np.ndarray:
        """
        A broadcasting version of *-product between a n-rank tensor and 1-rank array.
        In this class it is used with depth=0, 1, 2.
        """
        idx = "ijklmnopq"
        einsum_text = f"{idx[:depth]}...,...->{idx[:depth]}..."
        return np.einsum(einsum_text, array_in, weights, optimize=True)

    @staticmethod
    def eval_hessian(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        ndim = len(quadrule_list)
        basis_list = [q.basis_to_sample for q in quadrule_list]
        w_proj, z_proj, h_proj = NurbsOperations._project_nurbs_weights(
            quadrule_list, nurbs_weights, nders=2
        )
        u_ctrlpts_scaled = NurbsOperations._apply_scaling(
            u_ctrlpts, nurbs_weights, depth=1
        )
        hess_first = BsplineOperations.eval_hessian(quadrule_list, u_ctrlpts_scaled)
        hess_first = NurbsOperations._apply_scaling(hess_first, w_proj, depth=3)
        hess = np.zeros_like(hess_first, dtype=float)

        for i in range(hess.shape[0]):
            for j in range(hess.shape[1]):
                alpha_list = np.zeros(ndim, dtype=int)
                alpha_list[j] = 1
                for k in range(hess.shape[2]):
                    beta_list = np.zeros(ndim, dtype=int)
                    beta_list[k] = 1
                    hess_second_1_ijk = (
                        MatrixFree.apply(
                            [basis[a] for basis, a in zip(basis_list, alpha_list)],
                            u_ctrlpts_scaled[i],
                            is_transpose=False,
                        )
                        * z_proj[k]
                    )
                    hess_second_2_ijk = (
                        MatrixFree.apply(
                            [basis[b] for basis, b in zip(basis_list, beta_list)],
                            u_ctrlpts_scaled[i],
                            is_transpose=False,
                        )
                        * z_proj[j]
                    )
                    hess_third_1_ij = (
                        MatrixFree.apply(
                            [basis[0] for basis in basis_list],
                            u_ctrlpts_scaled[i],
                            is_transpose=False,
                        )
                        * z_proj[j]
                        * z_proj[k]
                        / w_proj
                    )
                    hess_third_2_ij = (
                        MatrixFree.apply(
                            [basis[0] for basis in basis_list],
                            u_ctrlpts_scaled[i],
                            is_transpose=False,
                        )
                        * h_proj[j][k]
                        / w_proj
                    )
                    hess[i][j][k] = (
                        hess_first[i][j][k]
                        - hess_second_1_ijk
                        - hess_second_2_ijk
                        + 2 * hess_third_1_ij
                        - hess_third_2_ij
                    )
        return hess

    @staticmethod
    def eval_jacobien(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        basis_list = [q.basis_to_sample for q in quadrule_list]
        w_proj, z_proj, _ = NurbsOperations._project_nurbs_weights(
            quadrule_list, nurbs_weights, nders=1
        )
        u_ctrlpts_scaled = NurbsOperations._apply_scaling(
            u_ctrlpts, nurbs_weights, depth=1
        )
        jac_first = BsplineOperations.eval_jacobien(quadrule_list, u_ctrlpts_scaled)
        jac_first = NurbsOperations._apply_scaling(jac_first, w_proj, depth=2)
        jac = np.zeros_like(jac_first, dtype=float)
        for i in range(jac.shape[0]):
            for j in range(jac.shape[1]):
                jac_second_ij = (
                    MatrixFree.apply(
                        [basis[0] for basis in basis_list],
                        u_ctrlpts_scaled[i],
                        is_transpose=False,
                    )
                    * z_proj[j]
                )
                jac[i][j] = jac_first[i][j] - jac_second_ij
        return jac

    @staticmethod
    def interpolate_meshgrid(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        w_proj = NurbsOperations._project_nurbs_weights(
            quadrule_list, nurbs_weights, nders=0
        )[0]
        u_ctrlpts_scaled = NurbsOperations._apply_scaling(
            u_ctrlpts, nurbs_weights, depth=1
        )
        u_interp = BsplineOperations.interpolate_meshgrid(
            quadrule_list, u_ctrlpts_scaled
        )
        return NurbsOperations._apply_scaling(u_interp, w_proj, depth=1)

    @staticmethod
    def spkron_product_quadrature(
        quadrule_list: List[IGAQuadratureRule],
        idx_list: list,
        product_type: Literal["basis", "weights"],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> sp.csr_array:
        return BsplineOperations.spkron_product_quadrature(
            quadrule_list, idx_list, product_type, nurbs_weights
        )

    @staticmethod
    def assemble_scalar_u_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        allow_lumping: bool,
        nurbs_weights: np.ndarray = np.array([]),
    ):
        raise NotImplementedError("Not implemented")

    @staticmethod
    def assemble_scalar_gradu_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ):
        raise NotImplementedError("Not implemented")

    @staticmethod
    def assemble_scalar_u_force(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        w_proj = NurbsOperations._project_nurbs_weights(
            quadrule_list, nurbs_weights, nders=0
        )[0]
        coefficients_copy = NurbsOperations._apply_scaling(
            coefficients, w_proj, depth=1
        )
        array_out = BsplineOperations.assemble_scalar_u_force(
            quadrule_list, coefficients_copy
        )
        return NurbsOperations._apply_scaling(array_out, nurbs_weights, depth=1)

    @staticmethod
    def compute_mf_scalar_u_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        allow_lumping: bool,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        w_proj = NurbsOperations._project_nurbs_weights(
            quadrule_list, nurbs_weights, nders=0
        )[0]
        if allow_lumping:
            array_in_copy = np.copy(array_in)
            coefficients_copy = NurbsOperations._apply_scaling(
                coefficients, w_proj, depth=0
            )
        else:
            array_in_copy = NurbsOperations._apply_scaling(
                array_in, nurbs_weights, depth=0
            )
            coefficients_copy = NurbsOperations._apply_scaling(
                coefficients, w_proj**2, depth=0
            )

        array_out = BsplineOperations.compute_mf_scalar_u_v(
            quadrule_list,
            coefficients_copy,
            array_in_copy,
            allow_lumping,
            enable_spacetime=enable_spacetime,
            time_ders=time_ders,
        )
        if allow_lumping:
            return array_out
        else:
            return NurbsOperations._apply_scaling(array_out, nurbs_weights, depth=0)

    @staticmethod
    def compute_mf_scalar_gradu_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        w_proj, z_proj, _ = NurbsOperations._project_nurbs_weights(
            quadrule_list, nurbs_weights, nders=1
        )
        array_in_copy = NurbsOperations._apply_scaling(array_in, nurbs_weights, depth=0)
        # First term
        coefficients_copy = NurbsOperations._apply_scaling(
            coefficients, w_proj**2, depth=2
        )
        array_out = BsplineOperations.compute_mf_scalar_gradu_gradv(
            quadrule_list,
            coefficients_copy,
            array_in_copy,
            enable_spacetime=enable_spacetime,
            time_ders=time_ders,
        )
        # Second term
        coefficients_copy = np.einsum(
            "ij...,i...,j...->...", coefficients, z_proj, z_proj, optimize=True
        )
        array_out += BsplineOperations.compute_mf_scalar_u_v(
            quadrule_list,
            coefficients_copy,
            array_in_copy,
            allow_lumping=False,
            enable_spacetime=enable_spacetime,
            time_ders=time_ders,
        )
        # Third term
        coefficients_copy = np.einsum(
            "ij...,i...,...->j...", coefficients, z_proj, w_proj, optimize=True
        )
        array_out -= BsplineOperations.compute_mf_scalar_u_gradv(
            quadrule_list,
            coefficients_copy,
            array_in_copy,
            enable_spacetime=enable_spacetime,
            time_ders=time_ders,
        )
        # Fourth term
        coefficients_copy = np.einsum(
            "ij...,j...,...->i...", coefficients, z_proj, w_proj, optimize=True
        )
        array_out -= BsplineOperations.compute_mf_scalar_gradu_v(
            quadrule_list,
            coefficients_copy,
            array_in_copy,
            enable_spacetime=enable_spacetime,
            time_ders=time_ders,
        )
        return NurbsOperations._apply_scaling(array_out, nurbs_weights, depth=0)

    @staticmethod
    def _compute_mf_scalar_gradu_v_or_u_grad_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
        mf_type: str = "",
    ) -> np.ndarray:
        assert mf_type in ["gradu_v", "u_gradv"]
        w_proj, z_proj, _ = NurbsOperations._project_nurbs_weights(
            quadrule_list, nurbs_weights, nders=1
        )
        array_in_copy = NurbsOperations._apply_scaling(array_in, nurbs_weights, depth=0)
        oper = (
            BsplineOperations.compute_mf_scalar_gradu_v
            if mf_type == "gradu_v"
            else BsplineOperations.compute_mf_scalar_u_gradv
        )
        # First term
        coefficients_copy = NurbsOperations._apply_scaling(
            coefficients, w_proj**2, depth=1
        )
        array_out = oper(
            quadrule_list,
            coefficients_copy,
            array_in_copy,
            enable_spacetime=enable_spacetime,
            time_ders=time_ders,
        )
        # Second term
        coefficients_copy = np.einsum(
            "i...,i...,...->...", coefficients, z_proj, w_proj, optimize=True
        )
        array_out -= BsplineOperations.compute_mf_scalar_u_v(
            quadrule_list,
            coefficients_copy,
            array_in_copy,
            allow_lumping=False,
            enable_spacetime=enable_spacetime,
            time_ders=time_ders,
        )
        return NurbsOperations._apply_scaling(array_out, nurbs_weights, depth=0)

    @staticmethod
    def compute_mf_scalar_gradu_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        return NurbsOperations._compute_mf_scalar_gradu_v_or_u_grad_v(
            quadrule_list,
            coefficients,
            array_in,
            enable_spacetime,
            time_ders,
            nurbs_weights,
            mf_type="gradu_v",
        )

    @staticmethod
    def compute_mf_scalar_u_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        return NurbsOperations._compute_mf_scalar_gradu_v_or_u_grad_v(
            quadrule_list,
            coefficients,
            array_in,
            enable_spacetime,
            time_ders,
            nurbs_weights,
            mf_type="u_gradv",
        )
