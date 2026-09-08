from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from .matrix_free import MatrixFree
from .core import Operations
from typing import List, Sequence, Literal, Union
from scipy import sparse as sp
import numpy as np


class BsplineOperations(Operations):
    @staticmethod
    def _get_nbofrows(
        quadrule_list: List[IGAQuadratureRule],
    ) -> np.ndarray:
        return np.array([q.nbctrlpts for q in quadrule_list], dtype=int)

    @staticmethod
    def _get_nbofcols(
        quadrule_list: List[IGAQuadratureRule], is_sample: bool
    ) -> np.ndarray:
        if is_sample:
            return np.array([len(q.knots_to_sample) for q in quadrule_list], dtype=int)
        else:
            return np.array([q.nbquadpts for q in quadrule_list], dtype=int)

    @staticmethod
    def eval_jacobien(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        basis_list = [q.basis_to_sample for q in quadrule_list]

        # Number of scalar fields
        nm = np.shape(u_ctrlpts)[0]
        # Number of parameteric variables xi = (xi_1, xi_2, ...)
        ndim = len(quadrule_list)

        nc_list = BsplineOperations._get_nbofcols(quadrule_list, is_sample=True)
        jac = np.zeros((nm, ndim, np.prod(nc_list)))  # J_ij = d u_i / d xi_j
        for i in range(nm):
            for j in range(ndim):
                beta_list = np.zeros(ndim, dtype=int)
                beta_list[j] = 1
                jac[i][j] = MatrixFree.apply(
                    [basis[beta] for basis, beta in zip(basis_list, beta_list)],
                    u_ctrlpts[i],
                    is_transpose=False,
                )

        return jac

    @staticmethod
    def eval_hessian(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        basis_list = [q.eval_basis(q.knots_to_sample, nders=2) for q in quadrule_list]

        # Number of scalar fields
        nm = np.shape(u_ctrlpts)[0]
        # Number of parameteric variables xi = (xi_1, xi_2, ...)
        ndim = len(quadrule_list)

        nc_list = BsplineOperations._get_nbofcols(quadrule_list, is_sample=True)
        hess = np.zeros(
            (nm, ndim, ndim, np.prod(nc_list))
        )  # H_ijk = d^2 u_i / d xi_j / d xi_k
        for i in range(nm):
            for j in range(ndim):
                alpha_list = np.zeros(ndim, dtype=int)
                alpha_list[j] = 1
                for k in range(ndim):
                    beta_list = np.zeros(ndim, dtype=int)
                    beta_list[k] = 1
                    zeta_list = alpha_list + beta_list
                    hess[i][j][k] = MatrixFree.apply(
                        [basis[zeta] for basis, zeta in zip(basis_list, zeta_list)],
                        u_ctrlpts[i],
                        is_transpose=False,
                    )
        return hess

    @staticmethod
    def interpolate_meshgrid(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        basis_list = [q.basis_to_sample for q in quadrule_list]

        # Number of scalar fields
        nm = np.shape(u_ctrlpts)[0]
        nc_list = BsplineOperations._get_nbofcols(quadrule_list, is_sample=True)
        u_interp = np.zeros((nm, np.prod(nc_list)))
        for i in range(nm):
            u_interp[i] = MatrixFree.apply(
                [basis[0] for basis in basis_list], u_ctrlpts[i], is_transpose=False
            )

        return u_interp

    @staticmethod
    def spkron_product_quadrature(
        quadrule_list: List[IGAQuadratureRule],
        idx_list: list,
        product_type: Literal["basis", "weights"],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> sp.csr_array:
        assert len(quadrule_list) == len(idx_list) and len(idx_list) > 0
        if product_type == "basis":
            matrix_list = [q.basis[idx] for q, idx in zip(quadrule_list, idx_list)]
        elif product_type == "weights":
            matrix_list = [q.weights[idx] for q, idx in zip(quadrule_list, idx_list)]
        else:
            raise NotImplementedError()
        matrix = matrix_list[0]
        for curr in matrix_list[1:]:
            matrix = sp.kron(curr, matrix)
        return sp.csr_array(matrix.tocsr())

    @staticmethod
    def assemble_scalar_u_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        allow_lumping: bool,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> Union[np.ndarray, sp.csr_array]:
        ndim = len(quadrule_list)
        if allow_lumping:
            matrix = MatrixFree.apply(
                [q.weights[0] for q in quadrule_list],
                coefficients,
                is_transpose=False,
            )
        else:
            zero_list = np.zeros(ndim, dtype=int).tolist()
            tmp1 = BsplineOperations.spkron_product_quadrature(
                quadrule_list, zero_list, product_type="basis"
            )
            tmp2 = sp.diags(coefficients) @ tmp1
            matrix = sp.csr_array(
                BsplineOperations.spkron_product_quadrature(
                    quadrule_list, zero_list, product_type="weights"
                )
                @ tmp2
            )
        return matrix

    @staticmethod
    def assemble_scalar_gradu_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> sp.csr_array:
        ndim = len(quadrule_list)
        nr_list = BsplineOperations._get_nbofrows(quadrule_list)
        matrix = sp.csr_array((np.prod(nr_list), np.prod(nr_list)))
        for j in range(ndim):
            beta_list = np.zeros(ndim, dtype=int)
            beta_list[j] = 1
            tmp1 = BsplineOperations.spkron_product_quadrature(
                quadrule_list, beta_list.tolist(), product_type="basis"
            )
            for i in range(ndim):
                alpha_list = np.zeros(ndim, dtype=int)
                alpha_list[i] = 1
                zeta_list = beta_list + 2 * alpha_list
                tmp2 = sp.diags(coefficients[i][j]) @ tmp1
                matrix += (
                    BsplineOperations.spkron_product_quadrature(
                        quadrule_list, zeta_list.tolist(), product_type="weights"
                    )
                    @ tmp2
                )
        return matrix

    @staticmethod
    def assemble_scalar_u_force(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        nm = np.shape(coefficients)[0]
        nr_list = BsplineOperations._get_nbofrows(quadrule_list)
        array_out = np.zeros((nm, np.prod(nr_list)))
        for i in range(nm):
            array_out[i] = MatrixFree.apply(
                [q.weights[0] for q in quadrule_list],
                coefficients[i],
                is_transpose=False,
            )
        return array_out

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

        if allow_lumping:
            matrix_lumped = MatrixFree.apply(
                [q.weights[0] for q in quadrule_list],
                coefficients,
                is_transpose=False,
            )
            array_out = matrix_lumped * array_in
        else:
            if enable_spacetime:
                assert isinstance(
                    time_ders, (tuple, list)
                ), "Derivatives of time functions should be tuple"
            ndim = len(quadrule_list)
            zero_list = np.zeros(ndim, dtype=int)
            if enable_spacetime:
                zero_list[-1] = time_ders[1]
            array_tmp = MatrixFree.apply(
                [q.basis[beta] for q, beta in zip(quadrule_list, zero_list)],
                array_in,
                is_transpose=False,
            )

            if enable_spacetime:
                zero_list[-1] = time_ders[0]
            array_out = MatrixFree.apply(
                [q.weights[alpha] for q, alpha in zip(quadrule_list, zero_list)],
                array_tmp * coefficients,
                is_transpose=False,
            )
        return array_out

    @staticmethod
    def compute_mf_scalar_gradu_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        if enable_spacetime:
            assert isinstance(
                time_ders, (tuple, list)
            ), "Derivatives of time functions should be tuple"
        ndim = len(quadrule_list)
        sp_ndim = ndim - 1 if enable_spacetime else ndim
        array_out = np.zeros_like(array_in, dtype=float)
        for j in range(sp_ndim):
            beta_list = np.zeros(ndim, dtype=int)
            beta_list[j] = 1
            if enable_spacetime:
                beta_list[-1] = time_ders[1]
            array_tmp = MatrixFree.apply(
                [q.basis[beta] for q, beta in zip(quadrule_list, beta_list)],
                array_in,
                is_transpose=False,
            )
            for i in range(sp_ndim):
                alpha_list = np.zeros(ndim, dtype=int)
                alpha_list[i] = 1
                if enable_spacetime:
                    alpha_list[-1] = time_ders[0]
                zeta_list = beta_list + 2 * alpha_list
                array_out += MatrixFree.apply(
                    [q.weights[zeta] for q, zeta in zip(quadrule_list, zeta_list)],
                    array_tmp * coefficients[i][j],
                    is_transpose=False,
                )
        return array_out

    @staticmethod
    def compute_mf_scalar_gradu_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        if enable_spacetime:
            assert isinstance(
                time_ders, (tuple, list)
            ), "Derivatives of time functions should be tuple"
        ndim = len(quadrule_list)
        sp_ndim = ndim - 1 if enable_spacetime else ndim
        array_out = np.zeros_like(array_in, dtype=float)

        beta_list = np.zeros(ndim, dtype=int)
        if enable_spacetime:
            beta_list[-1] = time_ders[1]
        array_tmp = MatrixFree.apply(
            [q.basis[beta] for q, beta in zip(quadrule_list, beta_list)],
            array_in,
            is_transpose=False,
        )
        for i in range(sp_ndim):
            alpha_list = np.zeros(ndim, dtype=int)
            alpha_list[i] = 1
            if enable_spacetime:
                alpha_list[-1] = time_ders[0]
            zeta_list = beta_list + 2 * alpha_list
            array_out += MatrixFree.apply(
                [q.weights[zeta] for q, zeta in zip(quadrule_list, zeta_list)],
                array_tmp * coefficients[i],
                is_transpose=False,
            )
        return array_out

    @staticmethod
    def compute_mf_scalar_u_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        if enable_spacetime:
            assert isinstance(
                time_ders, (tuple, list)
            ), "Derivatives of time functions should be tuple"
        ndim = len(quadrule_list)
        nbcols = BsplineOperations._get_nbofcols(quadrule_list, is_sample=False)
        sp_ndim = ndim - 1 if enable_spacetime else ndim
        array_tmp = np.zeros(np.prod(nbcols))
        for i in range(sp_ndim):
            beta_list = np.zeros(ndim, dtype=int)
            beta_list[i] = 1
            if enable_spacetime:
                beta_list[-1] = time_ders[1]
            array_tmp += (
                MatrixFree.apply(
                    [q.basis[beta] for q, beta in zip(quadrule_list, beta_list)],
                    array_in,
                    is_transpose=False,
                )
                * coefficients[i]
            )

        alpha_list = np.zeros(ndim, dtype=int)
        if enable_spacetime:
            alpha_list[-1] = time_ders[0]
        array_out = MatrixFree.apply(
            [q.weights[alpha] for q, alpha in zip(quadrule_list, alpha_list)],
            array_tmp,
            is_transpose=False,
        )
        return array_out
    
    @staticmethod
    def _element_quad_indices_1d(q: IGAQuadratureRule, elem: int) -> np.ndarray:
        nq_per_elem = q.nbquadpts // q._nbelem
        return np.arange(elem * nq_per_elem, (elem + 1) * nq_per_elem)


    @staticmethod
    def _element_active_ctrlpts_1d(q: IGAQuadratureRule, elem: int) -> np.ndarray:
        rows = BsplineOperations._element_quad_indices_1d(q, elem)
        B0 = q.basis[0][rows, :]
        return np.unique(B0.nonzero()[1]).astype(int)


    @staticmethod
    def extract_1d_element_basis(
        q: IGAQuadratureRule,
        elem: int,
        der: int,
    ) -> sp.csr_array:
        rows = BsplineOperations._element_quad_indices_1d(q, elem)
        active = BsplineOperations._element_active_ctrlpts_1d(q, elem)

        B_e = q.basis[der][rows, :][:, active]

        return sp.csr_array(B_e.tocsr())


    @staticmethod
    def extract_1d_element_weights(
        q: IGAQuadratureRule,
        elem: int,
        der: int,
    ) -> sp.csr_array:
        rows = BsplineOperations._element_quad_indices_1d(q, elem)
        active = BsplineOperations._element_active_ctrlpts_1d(q, elem)

        W_e = q.weights[der][active, :][:, rows]

        return sp.csr_array(W_e.tocsr())


    @staticmethod
    def element_connectivity(
        quadrule_list: List[IGAQuadratureRule],
        elem_idx: tuple,
    ) -> np.ndarray:

        conn_1d = [
            BsplineOperations._element_active_ctrlpts_1d(q, elem)
            for q, elem in zip(quadrule_list, elem_idx)
        ]

        grids = np.meshgrid(*conn_1d, indexing="ij")

        nr_list = BsplineOperations._get_nbofrows(quadrule_list)

        return np.ravel_multi_index(
            [g.ravel(order="F") for g in grids],
            dims=tuple(nr_list),
            order="F",
        )


    @staticmethod
    def spkron_product_quadrature_element(
        quadrule_list: List[IGAQuadratureRule],
        elem_idx: tuple,
        der_list: list,
        product_type: Literal["basis", "weights"],
    ) -> sp.csr_array:

        matrix_list = []

        for q, elem, der in zip(quadrule_list, elem_idx, der_list):

            if product_type == "basis":
                matrix_list.append(
                    BsplineOperations.extract_1d_element_basis(q, elem, der)
                )

            elif product_type == "weights":
                matrix_list.append(
                    BsplineOperations.extract_1d_element_weights(q, elem, der)
                )

            else:
                raise NotImplementedError()

        matrix = matrix_list[0]

        for curr in matrix_list[1:]:
            matrix = sp.kron(curr, matrix, format="csr")

        return sp.csr_array(matrix.tocsr())


    @staticmethod
    def extract_element_coefficients(
        quadrule_list: List[IGAQuadratureRule],
        elem_idx: tuple,
        coefficients: np.ndarray,
    ) -> np.ndarray:

        coefficients = np.asarray(coefficients)

        nbquadpts_list = [q.nbquadpts for q in quadrule_list]
        ndim = len(quadrule_list)

        rows_1d = [
            BsplineOperations._element_quad_indices_1d(q, elem)
            for q, elem in zip(quadrule_list, elem_idx)
        ]

        grids = np.meshgrid(*rows_1d, indexing="ij")

        idx_quad = np.ravel_multi_index(
            [g.ravel(order="F") for g in grids],
            dims=tuple(nbquadpts_list),
            order="F",
        )
        

        if coefficients.shape == tuple(nbquadpts_list):
            return coefficients.ravel(order="F")[idx_quad]

        if coefficients.size == np.prod(nbquadpts_list):
            return coefficients.ravel(order="F")[idx_quad]

        if coefficients.ndim > ndim and coefficients.shape[-ndim:] == tuple(nbquadpts_list):
            return coefficients.reshape(-1, np.prod(nbquadpts_list), order="F")[:, idx_quad]

        raise ValueError(
            f"Shape coefficients incompatible : {coefficients.shape}. "
            f"Attendu {tuple(nbquadpts_list)} ou taille {np.prod(nbquadpts_list)}."
        )

    @staticmethod
    def assemble_scalar_u_v_elementwise(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> sp.csr_array:

        ndim = len(quadrule_list)
        nr_list = BsplineOperations._get_nbofrows(quadrule_list)

        matrix = sp.lil_array((np.prod(nr_list), np.prod(nr_list)))

        nb_elements_list = [q._nbelem for q in quadrule_list]
        zero_list = np.zeros(ndim, dtype=int).tolist()

        for elem_idx in np.ndindex(*nb_elements_list):

            idx = BsplineOperations.element_connectivity(quadrule_list, elem_idx)

            B_e = BsplineOperations.spkron_product_quadrature_element(
                quadrule_list,
                elem_idx,
                zero_list,
                product_type="basis",
            )

            W_e = BsplineOperations.spkron_product_quadrature_element(
                quadrule_list,
                elem_idx,
                zero_list,
                product_type="weights",
            )

            coeff_e = BsplineOperations.extract_element_coefficients(
                quadrule_list,
                elem_idx,
                coefficients,
            )

            coeff_e = np.asarray(coeff_e).ravel(order="C")

            M_e = W_e @ (sp.diags(coeff_e) @ B_e)

            for a, I in enumerate(idx):
                for b, J in enumerate(idx):
                    matrix[I, J] += M_e[a, b]

        return matrix.tocsr()
    

    @staticmethod
    def assemble_scalar_gradu_gradv_elementwise(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> sp.csr_array:

        ndim = len(quadrule_list)
        nr_list = BsplineOperations._get_nbofrows(quadrule_list)

        matrix = sp.lil_array((np.prod(nr_list), np.prod(nr_list)))

        nb_elements_list = [q._nbelem for q in quadrule_list]

        for elem_idx in np.ndindex(*nb_elements_list):

            idx = BsplineOperations.element_connectivity(quadrule_list, elem_idx)

            for j in range(ndim):

                beta_list = np.zeros(ndim, dtype=int)
                beta_list[j] = 1

                B_e = BsplineOperations.spkron_product_quadrature_element(
                    quadrule_list,
                    elem_idx,
                    beta_list.tolist(),
                    product_type="basis",
                )

                for i in range(ndim):

                    alpha_list = np.zeros(ndim, dtype=int)
                    alpha_list[i] = 1

                    zeta_list = beta_list + 2 * alpha_list

                    W_e = BsplineOperations.spkron_product_quadrature_element(
                        quadrule_list,
                        elem_idx,
                        zeta_list.tolist(),
                        product_type="weights",
                    )

                    coeff_global = coefficients[i][j]

                    coeff_e = BsplineOperations.extract_element_coefficients(
                        quadrule_list,
                        elem_idx,
                        coeff_global,
                    )

                    coeff_e = np.asarray(coeff_e).ravel(order="C")

                    K_e = W_e @ (sp.diags(coeff_e) @ B_e)

                    for a, I in enumerate(idx):
                        for b, J in enumerate(idx):
                            matrix[I, J] += K_e[a, b]

        return matrix.tocsr()

    @staticmethod
    def assemble_scalar_u_force_elementwise(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:

        coefficients = np.asarray(coefficients)

        if coefficients.ndim == len(quadrule_list):
            coefficients = coefficients.reshape(1, *coefficients.shape)

        nm = coefficients.shape[0]
        ndim = len(quadrule_list)
        nr_list = BsplineOperations._get_nbofrows(quadrule_list)

        array_out = np.zeros((nm, np.prod(nr_list)))

        nb_elements_list = [q._nbelem for q in quadrule_list]
        zero_list = np.zeros(ndim, dtype=int).tolist()

        for elem_idx in np.ndindex(*nb_elements_list):

            idx = BsplineOperations.element_connectivity(quadrule_list, elem_idx)

            W_e = BsplineOperations.spkron_product_quadrature_element(
                quadrule_list,
                elem_idx,
                zero_list,
                product_type="weights",
            )

            for i in range(nm):

                coeff_e = BsplineOperations.extract_element_coefficients(
                    quadrule_list,
                    elem_idx,
                    coefficients[i],
                )

                coeff_e = np.asarray(coeff_e).ravel(order="C")

                F_e = W_e @ coeff_e

                for a, I in enumerate(idx):
                    array_out[i, I] += F_e[a]

        return array_out

