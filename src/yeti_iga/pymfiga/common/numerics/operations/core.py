from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from typing import List, Sequence, Literal, Union
from abc import ABC, abstractmethod
from scipy import sparse as sp
import numpy as np


class Operations(ABC):

    @staticmethod
    @abstractmethod
    def eval_hessian(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        """
        Let say the 1-rank tensor field u is represented as u = sum_A N_A u_A,
        where N_A are NURBS (or B-spline) functions in the parametric space
        and u_A are the 'control points' of the field on the physical space.
        The Hessian is defined as H_ABC = d^2 u_A / d x_B / d x_C. Here x_B and x_C
        are the components of the parametric space.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            u_ctrlpts (np.ndarray): the 'control points' of the vector field u
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis
        Returns:
            The Hessian of the 1-rank field u evaluated the points given in
            the list of quadrature rules (by default quadrature points)
        """
        pass

    @staticmethod
    @abstractmethod
    def eval_jacobien(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        """
        Let say the 1-rank tensor field u is represented as u = sum_A N_A u_A,
        where N_A are NURBS (or B-spline) functions in the parametric space
        and u_A are the 'control points' of the field on the physical space.
        The jacobien is defined as J_AB = d u_A / d x_B. Here x_B are the
        components of the parametric space.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            u_ctrlpts (np.ndarray): the 'control points' of the vector field u
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis
        Returns:
            The jacobien of the 1-rank field u evaluated the points given in
            the list of quadrature rules (by default quadrature points)
        """
        pass

    @staticmethod
    @abstractmethod
    def interpolate_meshgrid(
        quadrule_list: List[IGAQuadratureRule],
        u_ctrlpts: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        """
        Let say the 1-rank tensor field u is represented as u = sum_A N_A u_A,
        where N_A are NURBS (or B-spline) functions in the parametric space
        and u_A are the 'control points' of the field on the physical space.
        The interpolation is just the projection of u on a set of given points.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            u_ctrlpts (np.ndarray): the 'control points' of the vector field u
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis
        Returns:
            The interpolation of the 1-rank field u at the points given in
            the list of quadrature rules (by default quadrature points)
        """
        pass

    @staticmethod
    @abstractmethod
    def spkron_product_quadrature(
        quadrule_list: List[IGAQuadratureRule],
        idx_list: list,
        product_type: Literal["basis", "weights"],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> sp.csr_array:
        """
        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            idx_list (List[int]): it helps to select the derivative or the function in
                quadrature list
            product_type (str): {"basis", "weights"}
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis

        Notes:
        ------
        This function is not optimize since it computes the kron product directly and
        don't use matrix-free algorithms. Avoid using it if possible.
        """
        pass

    @staticmethod
    @abstractmethod
    def assemble_scalar_u_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        allow_lumping: bool,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> Union[np.ndarray, sp.csr_array]:
        """
        Assembles the matrix M (also called mass matrix) that results from

        'int N_A(x) c(x) N_B(x) dx' <- in the hypercube [0, 1]^d

        where N_A and N_B are basis functions in the same parametric space.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            coefficients (np.ndarray): contains geometry and material properties.
                Note: it is a 0-rank tensor (or scalar) field
            allow_lumping (bool): if computes the whole matrix or it makes a lumping
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis

        Notes:
        ------
        This function is not optimize since it computes the kron product directly and
        don't use matrix-free algorithms. Avoid using it if possible. The only exception
        is when allow_lumping is True, then it call matrix-free algorithms.
        """
        pass

    @staticmethod
    @abstractmethod
    def assemble_scalar_gradu_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> sp.csr_array:
        """
        Assembles the matrix M (also called stiffness matrix) that results from

        'int grad(N_A(x)) [c(x) . grad(N_B(x))] dx' <- in the hypercube [0, 1]^d

        where N_A and N_B are basis functions in the same parametric space.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            coefficients (np.ndarray): contains geometry and material properties.
                Note: it is a d-rank tensor field
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis

        Notes:
        ------
        This function is not optimize since it computes the kron product directly and
        don't use matrix-free algorithms. Avoid using it if possible.
        """
        pass

    @staticmethod
    @abstractmethod
    def assemble_scalar_u_force(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        """
        Computes the force-like array F using matrix-free algorithms.
        Here the terms of F results from

        'int N_A(x) c(x) dx' <- in the hypercube [0, 1]^d

        where N_A are basis functions in the same parametric space.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            coefficients (np.ndarray): contains the evaluation of the force function
                at the quadrature points. It is a 0-rank tensor (or scalar) field
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis
        """
        pass

    @staticmethod
    @abstractmethod
    def compute_mf_scalar_u_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        allow_lumping: bool,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        """
        Computes the matrix-vector product M @ v using matrix-free algorithms.
        Here the terms of M (also called mass matrix) results from

        'int N_A(x) c(x) N_B(x) dx' <- in the hypercube [0, 1]^d

        where N_A and N_B are basis functions in the same parametric space.
        To generalize the method, basis funcitons may also include time deirvatives
        N_A(x, t) = N^p_t(t) x N_1(x_1) x N_2(x_2) x ... x N_d(x_d)
        if p = 0, there is no derivative, p = 1, it has been derived once, and so on.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            coefficients (np.ndarray): contains geometry and material properties.
                Note: it is a 0-rank tensor (or scalar) field
            array_in (np.ndarray): the vector to be multiplied
            enable_spacetime: True if the basis functions are space-time splines
            time_ders (Sequence[int]): Tuple of size 2, the first element sets the
                derivative of N_A, the second elements sets the derivative of N_B
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis
        """
        pass

    @staticmethod
    @abstractmethod
    def compute_mf_scalar_gradu_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        """
        Computes the matrix-vector product M @ v using matrix-free algorithms.
        Here the terms of M (also called stiffness matrix) results from

        'int grad(N_A(x)) [c(x) . grad(N_B(x))] dx' <- in the hypercube [0, 1]^d

        where N_A and N_B are basis functions in the same parametric space.
        To generalize the method, basis funcitons may also include time deirvatives
        N_A(x, t) = N^p_t(t) x N_1(x_1) x N_2(x_2) x ... x N_d(x_d)
        if p = 0, there is no derivative, p = 1, it has been derived once, and so on.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            coefficients (np.ndarray): contains geometry and material properties
                Note: it is a 2-rank tensor field
            array_in (np.ndarray): the vector to be multiplied
            enable_spacetime: True if the basis functions are space-time splines
            time_ders (Sequence[int]): Tuple of size 2, the first element sets the
                derivative of N_A, the second elements sets the derivative of N_B
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis
        """
        pass

    @staticmethod
    @abstractmethod
    def compute_mf_scalar_gradu_v(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        """
        Computes the matrix-vector product M @ v using matrix-free algorithms.
        Here the terms of M (also called advection matrix) from

        'int [grad(N_A(x)) . c(x)] . N_B(x) dx' <- in the hypercube [0, 1]^d

        where N_A and N_B are basis functions in the same parametric space.
        To generalize the method, basis funcitons may also include time deirvatives
        N_A(x, t) = N^p_t(t) x N_1(x_1) x N_2(x_2) x ... x N_d(x_d)
        if p = 0, there is no derivative, p = 1, it has been derived once, and so on.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            coefficients (np.ndarray): contains geometry and material properties
                Note: it is a 1-rank tensor field
            array_in (np.ndarray): the vector to be multiplied
            enable_spacetime: True if the basis functions are space-time splines
            time_ders (Sequence[int]): Tuple of size 2, the first element sets the
                derivative of N_A, the second elements sets the derivative of N_B
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis
        """
        pass

    @staticmethod
    @abstractmethod
    def compute_mf_scalar_u_gradv(
        quadrule_list: List[IGAQuadratureRule],
        coefficients: np.ndarray,
        array_in: np.ndarray,
        enable_spacetime: bool = False,
        time_ders: Sequence[int] = [],
        nurbs_weights: np.ndarray = np.array([]),
    ) -> np.ndarray:
        """
        Computes the matrix-vector product M @ v using matrix-free algorithms.
        Here the terms of M (also called advection matrix transposed) results from

        'int N_A(x) [c(x) . grad(N_B(x))] dx' <- in the hypercube [0, 1]^d

        where N_A and N_B are basis functions in the same parametric space.
        To generalize the method, basis funcitons may also include time deirvatives
        N_A(x, t) = N^p_t(t) x N_1(x_1) x N_2(x_2) x ... x N_d(x_d)
        if p = 0, there is no derivative, p = 1, it has been derived once, and so on.

        Args:
            quadrule_list (List[IGAQuadratureRule]): list of quadrature rules
                to build the basis functions and to compute the integral
            coefficients (np.ndarray): contains geometry and material properties
                Note: it is a 1-rank tensor field
            array_in (np.ndarray): the vector to be multiplied
            enable_spacetime: True if the basis functions are space-time splines
            time_ders (Sequence[int]): Tuple of size 2, the first element sets the
                derivative of N_A, the second elements sets the derivative of N_B
            nurbs_weights (np.ndarray): the weights of the B-spline basis to compute
                the NURBS basis
        """
        pass
