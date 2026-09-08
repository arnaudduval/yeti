from .operations import Operations
from typing import List, Optional, Any
from abc import ABC, abstractmethod
from scipy import sparse as sp
from scipy import linalg as sclin
import numpy as np


class IGAQuadratureRule(ABC):
    def __init__(self, degree: int, knotvector: np.ndarray, is_periodic: bool = False):
        # Public variables that need to be checked
        self._set_degree(degree)
        self._set_knotvector(knotvector)
        self._set_is_periodic(is_periodic)

        # Public variables that are fixed
        self._nbctrlpts = int(len(self.knotvector) - self.degree - 1) - (
            self.degree if self.is_periodic else 0
        )
        self._unique_kv = np.unique(
            self.knotvector[(self.knotvector >= 0.0) & (self.knotvector <= 1.0)]
        )
        self._nbelem = int(len(self.unique_kv) - 1)

        # Public variables that may vary
        self._quadrature_class = ""
        self._quadrature_type = ""
        self._quadpts: np.ndarray = np.array([])
        self._basis: List[sp.csr_array] = []
        self._weights: List[sp.csr_array] = []
        self._knots_to_sample = None
        self._basis_to_sample = None

        # Private variables
        self._lu_gram_matrix = None
        self._coo_indices: list = []
        self._coo_basis: np.ndarray = np.array([])
        self._coo_weigts: np.ndarray = np.array([])

    @property
    def degree(self):
        return self._degree

    def _set_degree(self, value):
        assert isinstance(value, (int, np.integer))
        self._degree = int(value)

    @property
    def knotvector(self):
        return self._knotvector

    def _set_knotvector(self, value):
        assert isinstance(value, (np.ndarray))
        self._knotvector = np.asarray(value)

    @property
    def is_periodic(self):
        return self._is_periodic

    def _set_is_periodic(self, value):
        assert isinstance(value, bool)
        self._is_periodic = value

    @property
    def quadrature_class(self) -> str:
        "Quadrature class"
        return self._quadrature_class

    @property
    def quadrature_type(self) -> str:
        "Quadrature type"
        return self._quadrature_type

    @property
    def quadpts(self):
        "Quadrature points"
        return self._quadpts

    @quadpts.setter
    def quadpts(self, value):
        assert isinstance(value, np.ndarray)
        self._quadpts = value

    @property
    def basis(self):
        "Quadrature basis"
        return self._basis

    @basis.setter
    def basis(self, value):
        assert isinstance(value, list)
        assert all(isinstance(x, sp.csr_array) for x in value)
        self._basis = value

    @property
    def weights(self):
        "Quadrature weights"
        return self._weights

    @weights.setter
    def weights(self, value):
        assert isinstance(value, list)
        assert all(isinstance(x, sp.csr_array) for x in value)
        self._weights = value

    @property
    def nbctrlpts(self):
        "Number of control points"
        return self._nbctrlpts

    @property
    def unique_kv(self):
        "Unique knots of knotvector"
        return self._unique_kv

    @property
    def nbelem(self):
        "Number of elements"
        return self._nbelem

    @property
    def nbquadpts(self):
        "Number of quadrature points"
        return len(self.quadpts)

    @property
    def basis_to_sample(self):
        if self._basis_to_sample is not None:
            return self._basis_to_sample
        return self._basis

    @property
    def knots_to_sample(self):
        if self._knots_to_sample is not None:
            return self._knots_to_sample
        return self._quadpts

    @knots_to_sample.setter
    def knots_to_sample(self, value):
        assert isinstance(value, np.ndarray)
        self._knots_to_sample = value
        self._basis_to_sample = Operations.eval_ders_basis_sparse(
            self.degree,
            self.knotvector,
            self._knots_to_sample,
            nders=1,
            is_periodic=self.is_periodic,
        )

    @property
    def max_h_size(self) -> float:
        "Maximal size among knot-spans"
        return float(np.max(np.abs(np.diff(self.unique_kv))))

    def clear_sample(self):
        "Erases the information (knots and basis) of the sample"
        self._knots_to_sample = None
        self._basis_to_sample = None

    def _set_coo_basis_weights(
        self, basis: np.ndarray, weights: np.ndarray, indices: list
    ):
        self._coo_indices = indices
        self._coo_basis = basis
        self._coo_weights = weights

    def _assemble_csr_basis_weights(self):
        basis, weights = [], []
        indi, indj = self._coo_indices
        indi = indi % self.nbctrlpts

        for i in range(np.size(self._coo_basis, axis=1)):
            b = sp.coo_array(
                (self._coo_basis[:, i], (indj, indi)),
                shape=(self.nbquadpts, self.nbctrlpts),
            ).tocsr()
            b.eliminate_zeros()
            basis.append(b)

        for i in range(np.size(self._coo_weights, axis=1)):
            w = sp.coo_array(
                (self._coo_weights[:, i], (indi, indj)),
                shape=(self.nbctrlpts, self.nbquadpts),
            ).tocsr()
            w.eliminate_zeros()
            weights.append(w)
        self.basis, self.weights = basis, weights

    def eval_basis(self, knots_to_interp: np.ndarray, nders: int = 1):
        """
        Evaluates basis functions and its first derivative at given knots.
        Args:
            knots_to_interp (np.ndarray): the list of knots where
                the spline is evaluated
        """
        return Operations.eval_ders_basis_sparse(
            self.degree,
            self.knotvector,
            knots_to_interp,
            nders=nders,
            is_periodic=self.is_periodic,
        )

    def compute_dual_basis(
        self, knots_to_interp: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Evaluates dual basis functions at given knots.
        Note: dual basis is the inverse of Gram matrix in parametric space.
        Args:
            knots_to_interp (np.ndarray): the list of knots where
                the spline is evaluated. By default it is
                the quadrature points of the quadrature rule
        """

        # Evaluate B-spline basis at all points
        if knots_to_interp is not None:
            assert isinstance(knots_to_interp, np.ndarray)
            basis = self.eval_basis(knots_to_interp, nders=1)[0]
        else:
            basis = self.basis_to_sample[0].copy()

        assert sp.issparse(basis)

        # Compute dual basis
        if self._lu_gram_matrix is None:
            # Compute Gram matrix and inverse
            assert len(self.weights) > 0 and len(self.basis) > 0
            G = sp.csr_array(self.weights[0] @ self.basis[0]).toarray()  # Mass matrix
            self._lu_gram_matrix = sclin.lu_factor(G)

        luG = self._lu_gram_matrix
        dualbasis: np.ndarray = sclin.lu_solve(
            luG, basis.toarray().T
        )  # each row is lambda_i(xi)

        return dualbasis.T

    @abstractmethod
    def export_quadrature_rules(self) -> Any:
        "Initializes the quadrature class"
        raise NotImplementedError("Subclasses must implement this method")

    def __repr__(self) -> str:
        message = f""""
            \n {self.quadrature_class} QUADRATURE:
            type: {self.quadrature_type}
            degree: {self.degree}
            nb of functions: {self.nbctrlpts}
            nb of elements: {self.nbelem}
            max knot-span: {self.max_h_size:.2e}
            nb of quadrature points: {self.nbquadpts}
        """
        return message
