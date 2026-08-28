# -*- coding: utf-8 -*-

import numpy as np
from scipy import sparse as sp
from scipy.sparse import linalg as scsplin
from yeti_iga.pymfiga.common.numerics.operations import MatrixFree


class ScaledMassPreconditioner:
    """
    Préconditionneur de masse :
        M_tilde = D^{1/2} Dhat^{-1/2} Mhat Dhat^{-1/2} D^{1/2}

    Dans le solveur, on applique M_tilde^{-1}.
    """

    def __init__(self):
        
        """
        Initializes the scaled mass preconditioner.

        The constructor allocates the attributes required to build the
        preconditioner. No computation is performed during initialization.
        All quantities are computed afterwards by calling ``compute(model)``.

        Attributes
        ----------
        mass_matrices_1d : list
            One-dimensional parametric mass matrices M_k.

        diag_parametric_1d : list
            Diagonals of the one-dimensional parametric mass matrices,
            Dhat_k = diag(M_k).

        diag_true_mass : ndarray
            Diagonal of the physical mass matrix,
            D = diag(M).

        A_factors : list
            Normalized one-dimensional matrices
            A_k = Dhat_k^{-1/2} M_k Dhat_k^{-1/2}.

        A_solvers : list
            Sparse factorizations associated with the matrices A_k,
            used to efficiently apply A_k^{-1}.

        shape_by_dir : list
            Number of control points in each parametric direction.

        nb_dofs_per_node : int
            Number of degrees of freedom associated with each control point.

        nb_nodes : int
            Total number of control points of the patch.

        is_initialized : bool
            Indicates whether the preconditioner has been constructed by
            calling ``compute(model)``.
            """ 

        self.mass_matrices_1d = None
        self.diag_parametric_1d = None
        self.diag_true_mass = None

        self.A_factors = None
        self.A_solvers = None
        self.shape_by_dir = None

        self.nb_dofs_per_node = None
        self.nb_nodes = None
        self.is_initialized = False

    def compute(self, model):
        """
        Builds the scaled mass preconditioner for the given IGA model.

        This function initializes all the quantities required by the
        preconditioner. It computes:
            - the one-dimensional parametric mass matrices,
            - the diagonal of the physical mass matrix,
            - the normalized one-dimensional matrices used in the
              Kronecker-product formulation of the preconditioner.

        Once this function has been executed, the preconditioner is fully
        initialized and can be applied through
        ``apply_spatial_preconditioner()``.

        Parameters
        ----------
        model : ExplicitDynamicsModel
            IGA model containing the geometry, quadrature rules,
            material properties and matrix-free operators required to
            construct the preconditioner.
        """
        
        self.nb_dofs_per_node = model.nbvars
        self.nb_nodes = model.part.nbctrlpts_total

        self._compute_parametric_mass_matrices(model)
        self._compute_true_mass_diagonal(model)
        self._compute_scaled_1d_matrices()

        self.is_initialized = True
          
    def _compute_parametric_mass_matrices(self, model):
        """
        Computes the one-dimensional parametric mass matrices.

        For each parametric direction, the corresponding one-dimensional
        mass matrix M_k is extracted from the weighted quadrature operators.
        These matrices are the building blocks of the global parametric
        mass matrix

            M_hat = M_d ⊗ ... ⊗ M_1,

        which is never assembled explicitly.

        The number of control points in each parametric direction is also
        stored for the subsequent tensor-product operations.

        Parameters
        ----------
        model : ExplicitDynamicsModel
            IGA model providing the quadrature rules associated with the
            spatial discretization.
        """
        
        quadrule_list = model.part.quadrule_list

        self.mass_matrices_1d = [
            q.weights[0] @ q.basis[0]
            for q in quadrule_list
        ]

        self.shape_by_dir = [
            M.shape[0]
            for M in self.mass_matrices_1d
        ]
        
    def _compute_true_mass_diagonal(self, model):
        """
        Computes the diagonal of the physical mass matrix without assembling
        the global mass matrix.

        The diagonal entries

            D = diag(M)

        are obtained directly from the matrix-free operators by integrating
        only the terms associated with the diagonal of the mass matrix.
        The resulting scalar diagonal is then duplicated for each degree of
        freedom of the mechanical problem.

        Parameters
        ----------
        model : ExplicitDynamicsModel
            IGA model providing the quadrature rules, mass properties and
            matrix-free operators required to evaluate the physical mass
            diagonal.
        """

        if model._mass_property is None:
            model.compute_mass_property()

        quadrule_list = model.part.quadrule_list
        coefficients = model.mass_property

        diagonal_operators = []

        for q in quadrule_list:
            basis = q.basis[0]      
            weights = q.weights[0] 

            # opérateur diagonal 1D :
            # diag_i = sum_q weights[i,q] * basis[q,i] * coeff[q]
            diag_op = weights.multiply(basis.T)

            diagonal_operators.append(diag_op)

        diag_scalar = MatrixFree.apply(
            diagonal_operators,
            coefficients,
            is_transpose=False,
        )

        diag = np.zeros(model.get_size_of_arrays())
        diag = np.reshape(diag, (self.nb_dofs_per_node, self.nb_nodes))

        for dof in range(self.nb_dofs_per_node):
            diag[dof, :] = diag_scalar

        self.diag_true_mass = np.ravel(diag)

        if np.any(self.diag_true_mass <= 0.0):
            raise ValueError("True mass diagonal must be positive.")
            
    def _compute_scaled_1d_matrices(self):
        """
        Builds the normalized one-dimensional matrices used by the
        scaled mass preconditioner.

        For each parametric direction, the diagonal of the one-dimensional
        mass matrix is extracted to construct

            A_k = Dhat_k^{-1/2} M_k Dhat_k^{-1/2},

        where Dhat_k = diag(M_k).

        The corresponding sparse factorizations are also computed in order
        to efficiently apply A_k^{-1} during the application of the
        preconditioner.
        """

        self.A_factors = []
        self.A_solvers = []
        self.diag_parametric_1d = []

        for M in self.mass_matrices_1d:
            M = sp.csr_array(M)

            diag_M = M.diagonal()

            if np.any(diag_M <= 0.0):
                raise ValueError("Diagonal of parametric mass matrix must be positive.")

            self.diag_parametric_1d.append(diag_M)

            inv_sqrt_diag = 1.0 / np.sqrt(diag_M)
            D_inv_sqrt = sp.diags(inv_sqrt_diag)

            A = D_inv_sqrt @ M @ D_inv_sqrt
            A = sp.csc_array(A)

            solver_A = scsplin.factorized(A)

            self.A_factors.append(A)
            self.A_solvers.append(solver_A)
            
    def _apply_kron_inverse(self, x):
        """
        Applies the inverse of the Kronecker product of the normalized
        one-dimensional matrices.

        This method computes

            (A_d^{-1} ⊗ ... ⊗ A_1^{-1}) x

        without assembling the global Kronecker matrix. The input vector is
        reshaped as a tensor and each inverse A_k^{-1} is applied along the
        corresponding parametric direction.

        Parameters
        ----------
        x : ndarray
            Input vector associated with one scalar field.

        Returns
        -------
        ndarray
            Result of the Kronecker inverse application.
        """
        
        tensor = np.reshape(x, self.shape_by_dir[::-1])

        for axis, solver_A in enumerate(self.A_solvers[::-1]):
            tensor = np.moveaxis(tensor, axis, 0)
            old_shape = tensor.shape

            tensor = tensor.reshape(old_shape[0], -1)

            for j in range(tensor.shape[1]):
                tensor[:, j] = solver_A(tensor[:, j])

            tensor = tensor.reshape(old_shape)
            tensor = np.moveaxis(tensor, 0, axis)

        return np.ravel(tensor)

    def apply_spatial_preconditioner(self, x):
        """
        Applies the inverse of the scaled mass preconditioner.

        This method computes

            M_tilde^{-1} x
            = D^{-1/2} (A_d^{-1} ⊗ ... ⊗ A_1^{-1}) D^{-1/2} x,

        where D is the diagonal of the physical mass matrix and A_k are the
        normalized one-dimensional parametric mass matrices.

        Parameters
        ----------
        x : ndarray
            Input vector to be preconditioned.

        Returns
        -------
        ndarray
            Preconditioned vector.
        """
        
        if not self.is_initialized:
            raise RuntimeError("ScaledMassPreconditioner must be computed before use.")

        x = np.asarray(x)

        y = x / np.sqrt(self.diag_true_mass)

        y = np.reshape(y, (self.nb_dofs_per_node, self.nb_nodes))

        out = np.zeros_like(y)

        for dof in range(self.nb_dofs_per_node):
            out[dof, :] = self._apply_kron_inverse(y[dof, :])

        out = np.ravel(out)
        out = out / np.sqrt(self.diag_true_mass)

        return out