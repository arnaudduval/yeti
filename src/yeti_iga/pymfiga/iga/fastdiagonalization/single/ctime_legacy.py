from yeti_iga.pymfiga.common.base.enum import ParametricDirection
from yeti_iga.pymfiga.common.base.math import kron_nonzero_indices
from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from yeti_iga.pymfiga.common.numerics.operations import MatrixFree
from .core import Template
from .cspace import SpaceFD
from typing import Sequence, List, Optional
from scipy import linalg as sclin, sparse as sp
from time import time
import numpy as np
import logging

logger = logging.getLogger("SRC.FASTDIAG")


class SpTimeFD(Template):
    def __init__(self, spacefd: SpaceFD):
        self._space_preconditioner = spacefd
        self._sptm_free_nodes: list = []
        self._advection_time_corrector: Sequence[float] = []
        self._schur_adv: Optional[np.ndarray] = None
        self._schur_mass: Optional[np.ndarray] = None
        self._schur_VSL: Optional[np.ndarray] = None
        self._schur_VSR: Optional[np.ndarray] = None

    @property
    def nbDoFsPerNode(self):
        return self._space_preconditioner.nbDoFsPerNode

    @property
    def nonzeros_by_dir(self):
        return self._space_preconditioner.nonzeros_by_dir

    @property
    def indices_by_dir(self):
        return self._space_preconditioner.indices_by_dir

    @property
    def tm_nnz(self):
        return self.nonzeros_by_dir.get(ParametricDirection.TAU, 0)

    @property
    def space_preconditioner(self):
        return self._space_preconditioner

    @property
    def sptm_free_nodes(self):
        return self._sptm_free_nodes

    @property
    def advection_time_corrector(self):
        return self._advection_time_corrector

    @property
    def schur_adv(self):
        assert self._schur_adv is not None
        return self._schur_adv

    @property
    def schur_mass(self):
        assert self._schur_mass is not None
        return self._schur_mass

    @property
    def schur_VSL(self):
        assert self._schur_VSL is not None
        return self._schur_VSL

    @property
    def schur_VSR(self):
        assert self._schur_VSR is not None
        return self._schur_VSR

    def _propagate_nodes_in_time(self):
        indx_time = self.indices_by_dir[(0, ParametricDirection.TAU)]
        nnz_time = self.nonzeros_by_dir[ParametricDirection.TAU]
        space_nodes = self.space_preconditioner.space_free_nodes.copy()
        nnz_list = [nnz_time, self.space_preconditioner.sp_nnz]
        free_nodes = [[] for _ in range(self.nbDoFsPerNode)]
        for ii in range(self.nbDoFsPerNode):
            indices_ii_list = [indx_time, space_nodes[ii]]
            global_indices = kron_nonzero_indices(indices_ii_list, nnz_list)
            free_nodes[ii] = global_indices
        return free_nodes

    def compute_time_schurdecomposition(self, time_quadrule: IGAQuadratureRule):
        """
        Compute the Schur decomposition of the pencil (advection, mass) for the time direction.

        Args:
            time_quadrule (IGAQuadratureRule): Quadrature rule for the time direction.
        """
        q = Template.rewrite_quadrature([time_quadrule])[0]
        matrices: List[sp.csr_array] = [
            q.weights[0] @ q.basis[0],
            q.weights[1] @ q.basis[1],
        ]
        # We always assume that time is constraint at the begining
        nnz = q.nbctrlpts
        indices = np.arange(1, nnz, dtype=int).tolist()
        self.nonzeros_by_dir[ParametricDirection.TAU] = nnz
        self.indices_by_dir[(0, ParametricDirection.TAU)] = indices
        A = matrices[1].toarray()[np.ix_(indices, indices)]
        B = matrices[0].toarray()[np.ix_(indices, indices)]
        S, T, Q, Z = sclin.qz(A, B, output="complex")
        self._schur_adv = S
        self._schur_mass = T
        self._schur_VSL = Q
        self._schur_VSR = Z
        self._advection_time_corrector = [1.0] * self.nbDoFsPerNode
        self._sptm_free_nodes = self._propagate_nodes_in_time()

    def add_scalar_time_correctors(
        self,
        advection_corrector: List[float] = [],
    ):
        """
        It is important to remember that the (vanilla) preconditioner
        does not consider any information about geometry or material properties,
        then in order to improve the performance of the preconditioner, we can add some
        scalar correctors to the Schur decomposition of the pencil (advection, mass) for the time direction.

        Args:
            advection_corrector (List[float]): List of scalar correctors for the advection matrix.
                The length of the list should be equal to the number of DoFs per node.
        """

        def is_verified(entry):
            return isinstance(entry, (list, tuple)) and len(entry) == self.nbDoFsPerNode

        if is_verified(advection_corrector):
            self._advection_time_corrector = advection_corrector

    def apply_spacetime_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        """
        Apply the full space-time fast-diagonalization preconditioner to the input array.

        Args:
            array_in (np.ndarray): Input array to which the full space-time preconditioner will be applied
        """
        start = time()
        array_in = np.reshape(array_in, (self.nbDoFsPerNode, -1))
        adv_corrector = self.advection_time_corrector
        array_out = np.zeros_like(array_in)

        for i in range(self.nbDoFsPerNode):
            array1 = MatrixFree.apply(
                self.space_preconditioner.eigenvec_by_dir_space[i]
                + [self.schur_VSL.conj()],
                array_in[i, self.sptm_free_nodes[i]],
                is_transpose=True,
            )
            # NOTE: In the reshape we consider nnz_time - 1 because we assume that
            # the first node in time is constrained, then the Schur decomposition
            # is computed only for the free nodes in time, which are nnz_time - 1.
            # NOTE: Here, the order "F" has the meaning of slicing in time
            # so maybe it's not worth spending time in changing it to "C"
            nnz_time = self.nonzeros_by_dir[ParametricDirection.TAU]
            array1_reshape = np.reshape(array1, (-1, nnz_time - 1), order="F")
            array2_reshape = np.zeros_like(array1_reshape)

            # TODO: Could we apply multi-threading here since we have
            # to solve many independent triangular systems?
            for idx, row in enumerate(array1_reshape):
                mat = (
                    adv_corrector[i] * self.schur_adv
                    + self.space_preconditioner.space_eigenvalues[i][idx]
                    * self.schur_mass
                )
                array2_reshape[idx] = sclin.blas.ztrsv(mat, row, lower=False)

            array_out[i, self.sptm_free_nodes[i]] = np.real(
                MatrixFree.apply(
                    self.space_preconditioner.eigenvec_by_dir_space[i]
                    + [self.schur_VSR],
                    np.ravel(array2_reshape, order="F"),
                    is_transpose=False,
                )
            )

        logger.debug(
            f"Single patch fast-diagonalization in {time() - start:.2e} seconds"
        )
        return np.ravel(array_out)
