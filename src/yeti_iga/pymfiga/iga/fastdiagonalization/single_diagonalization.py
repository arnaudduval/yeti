from yeti_iga.pymfiga.common.base.enum import Constants
from yeti_iga.pymfiga.common.base.cls import BasePreconditioner
from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from yeti_iga.pymfiga.common.numerics.operations import MatrixFree
from .single.cspace import SpaceFD

# from .single.ctime_legacy import SpTimeFD
from .single.ctime_arrowhead import SpTimeFD
from typing import List, Sequence, Optional
import numpy as np
import logging

logger = logging.getLogger("SRC.FASTDIAG")


class SingleFastDiagonalization(BasePreconditioner):
    """
    This class implements the fast-diagonalization preconditioner for a single patch.
    """

    def __init__(self):
        self._space_fd: Optional[SpaceFD] = None
        self._sptime_fd: Optional[SpTimeFD] = None

    @property
    def space_preconditioner(self):
        assert self._space_fd is not None
        return self._space_fd

    @property
    def sptm_preconditioner(self):
        assert self._sptime_fd is not None
        return self._sptime_fd

    @property
    def global_space_eigenvalues_inverse(self) -> np.ndarray:
        if self._space_fd is None:
            raise ValueError("Space preconditioner not initialized")
        # TODO: maybe save in order to avoid constructing every time
        # However, space eigenvalues may vary
        nnz = self._space_fd.sp_nnz
        local_eig = []
        for free_nodes, eigenvalues in zip(
            self._space_fd.space_free_nodes, self._space_fd.space_eigenvalues
        ):
            if np.any(eigenvalues <= Constants.TINY):
                # Apply regularization
                logger.warning("Apply regularization")
                eigenvalues += Constants.TINY * np.max(eigenvalues)
            reference = np.zeros(nnz)
            reference[free_nodes] = 1.0 / eigenvalues
            local_eig.append(reference)
        eigenvalues = np.hstack(local_eig)
        return eigenvalues

    def matvec_space_eigenvectors(self, array_in: np.ndarray, is_transpose: bool):
        if self._space_fd is None:
            raise ValueError("Space preconditioner not initialized")
        array_in = np.reshape(array_in, (self._space_fd.nbDoFsPerNode, -1))
        array_out = np.zeros_like(array_in)
        for i in range(self._space_fd.nbDoFsPerNode):
            array_out[i, self._space_fd.space_free_nodes[i]] = MatrixFree.apply(
                self._space_fd.eigenvec_by_dir_space[i],
                array_in[i, self._space_fd.space_free_nodes[i]],
                is_transpose=is_transpose,
            )
        return np.ravel(array_out)

    def compute_space_eigendecomposition(
        self,
        space_quadrule_list: List[IGAQuadratureRule],
        space_table_dirichlet: np.ndarray,
    ):
        """
        Compute the eigendecomposition of the pencils (stiffness, mass)
        for each direction in space and for each DoFs per node.

        Args:
            space_quadrule_list (List[IGAQuadratureRule]): List of quadrature rules
                for each spatial direction.
            space_table_dirichlet (np.ndarray): A boolean array indicating the
                Dirichlet boundary conditions for each DoFs per node and direction.
        """
        self._space_fd = SpaceFD()
        self._space_fd.compute_space_eigendecomposition(
            space_quadrule_list, space_table_dirichlet
        )

    def compute_time_schurdecomposition(self, time_quadrule: IGAQuadratureRule):
        """
        Compute the Schur decomposition of the pencil (advection, mass) for the time direction.

        Args:
            time_quadrule (IGAQuadratureRule): Quadrature rule for the time direction.
        """
        if self._space_fd is None:
            raise ValueError("Space preconditioner not initialized")
        self._sptime_fd = SpTimeFD(self._space_fd)
        self._sptime_fd.compute_time_schurdecomposition(time_quadrule)

    def update_space_eigenvalues(self, scalar_coefs: Sequence[float]):
        """
        Let say that the matrix is a linear combination of
        the mass and stiffness matrices (in the physical space), i.e.
        A = scalar_coefs[0] * M + scalar_coefs[1] * K,
        Then the preconditioner is built following the
        same linear combination but with the eigenvalues, i.e.
            P = scalar_coefs[0] * I + scalar_coefs[1] * Lambda,
        where Lambda is the combination of the eigenvalues of
        the pencils (stiffness, mass) for each direction in space.

        Args:
            scalar_coefs (Sequence[float]):
                Coefficients for the linear combination
                of mass and stiffness in the physical space.
                The first entry corresponds to the mass,
                and the second entry corresponds to the stiffness.
        """
        if self._space_fd is None:
            raise ValueError("Space preconditioner not initialized")
        self._space_fd.update_space_eigenvalues(scalar_coefs)

    def add_scalar_space_time_correctors(
        self,
        mass_corrector: Optional[List[float]] = [],
        stiffness_corrector: Optional[List[np.ndarray]] = [],
        advection_corrector: Optional[List[float]] = [],
    ):
        """
        It is important to remember that the (vanilla) preconditioner
        does not consider any information about geometry or material properties,
        then in order to improve the performance of the preconditioner, we can add some
        scalar correctors to the eigenvalues of the pencils (stiffness, mass) for each direction in space
        and to the Schur decomposition of the pencil (advection, mass) for the time direction.

        Args:
            mass_corrector (List[float]): List of scalar correctors for the mass matrix.
                The length of the list should be equal to the number of DoFs per node.
            stiffness_corrector (List[np.ndarray]): List of scalar correctors for the stiffness matrix.
                Each entry should be an array of shape (nbdirs,) corresponding to the number of spatial directions.
                The length of the list should be equal to the number of DoFs per node.
            advection_corrector (List[float]): List of scalar correctors for the advection matrix.
                The length of the list should be equal to the number of DoFs per node.
        """
        if self._space_fd is not None:
            self._space_fd.add_scalar_space_correctors(
                mass_corrector=mass_corrector, stiffness_corrector=stiffness_corrector
            )
        if self._sptime_fd is not None:
            self._sptime_fd.add_scalar_time_correctors(
                advection_corrector=advection_corrector
            )

    def apply_spatial_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        """
        Apply the spatial part of the fast-diagonalization preconditioner to the input array.

        Args:
            array_in (np.ndarray): Input array to which the spatial preconditioner will be applied
        """
        if self._space_fd is None:
            raise ValueError("Space preconditioner not initialized")

        return self._space_fd.apply_spatial_preconditioner(array_in)

    def apply_spacetime_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        """
        Apply the full space-time fast-diagonalization preconditioner to the input array.

        Args:
            array_in (np.ndarray): Input array to which the full space-time preconditioner will be applied
        """
        if self._sptime_fd is None:
            raise ValueError("Space-Time preconditioner not initialized")
        return self._sptime_fd.apply_spacetime_preconditioner(array_in)

    def __repr__(self) -> str:
        if self._space_fd is None:
            raise ValueError("Space preconditioner not initialized")

        message = f""""
            Fast diagonalization with
            {self._space_fd.nbDoFsPerNode} DoFs per node
            Total number of nodes in space: {self._space_fd.sp_nnz}
            Number of nodes free per DoFs (space): {[len(nodes) for nodes in self._space_fd.space_free_nodes]}
        """

        if self._sptime_fd is not None:
            message += f"""
            Total number of nodes in time: {self._sptime_fd.tm_nnz}
            Number of nodes free per DoFs (space-time): {[len(nodes) for nodes in self._sptime_fd.sptm_free_nodes]}
            """
        return message
