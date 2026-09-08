from yeti_iga.pymfiga.common.base.enum import ParametricDirection
from yeti_iga.pymfiga.common.base.math import kron_nonzero_indices
from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from yeti_iga.pymfiga.common.numerics.operations import MatrixFree
from .helpers import solve_special_arrowhead, compute_special_schur
from .core import Template
from .cspace import SpaceFD
from typing import Sequence, List, Optional, Tuple
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
        self._eigenvec_time: np.ndarray = np.array([])
        self._eigenblocks_time: Sequence[np.ndarray] = []
        self._blocks_time = None

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
    def eigenvec_time(self):
        return self._eigenvec_time

    @property
    def eigenblocks_time(self):
        return self._eigenblocks_time

    @property
    def blocks_time(self):
        if self._blocks_time is None:
            self._blocks_time = self._compute_blocks()
        return self._blocks_time

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
        # q = self._rewrite_quadrature([time_quadrule])[0]
        q = time_quadrule
        matrices: List[sp.csr_array] = [
            q.weights[0] @ q.basis[0],
            q.weights[1] @ q.basis[1],
        ]
        # We always assume that time is constraint at the begining
        nnz = q.nbctrlpts
        indices = np.arange(1, nnz, dtype=int).tolist()
        self.nonzeros_by_dir[ParametricDirection.TAU] = nnz
        self.indices_by_dir[(0, ParametricDirection.TAU)] = indices
        Wt = matrices[1].toarray()[np.ix_(indices, indices)]
        Mt = matrices[0].toarray()[np.ix_(indices, indices)]

        # Extract blocks
        Wp = Wt[:-1, :-1]
        Mp = Mt[:-1, :-1]
        m = Mt[:-1, -1]
        w = np.atleast_2d(Wt[:-1, -1])

        # Compute eigen decomposition
        Dtp, Utp = sclin.eig(Wp, Mp)
        for i in range(Utp.shape[1]):
            u_i = Utp[:, i].copy()
            Utp[:, i] /= np.sqrt(np.abs(u_i.conj().T @ Mp @ u_i))

        # Compute [v, 1] and normalized it -> [k rho]
        v = np.linalg.solve(Mp, -m)
        v_1 = np.hstack([v, [1]])
        k_rho = v_1 / np.sqrt(v_1 @ (Mt @ v_1))
        k = np.atleast_2d(k_rho[:-1])

        # Assemble Ut: Ut^H @ Mt @ Ut = It
        Ut = np.block([[Utp, k.T], [np.zeros_like(k), k_rho[-1]]])

        # Compute blocks of Dt
        g = np.atleast_2d(Utp.conj().T @ (np.block([Wp, w.T]) @ k_rho))
        sigma = k_rho @ (Wt @ k_rho)

        # Assemble Dt: Ut^H @ Wt @ Ut = Dt
        Dt = np.block([[np.diag(Dtp), g.T], [-g.conj(), sigma]])

        # We only need to save Ut and the blocks of Dt (not necessary Dt)
        self._eigenvec_time = Ut
        self._eigenblocks_time = [np.diag(Dt), g.ravel()]

        self._advection_time_corrector = [1.0] * self.nbDoFsPerNode
        self._sptm_free_nodes = self._propagate_nodes_in_time()

    def _compute_blocks(self):
        adv_corr = self.advection_time_corrector
        space_eigvals = self.space_preconditioner.space_eigenvalues
        Dt, g = self.eigenblocks_time
        bfac_list: List[np.ndarray] = []
        H_list: List[np.ndarray] = []
        schur_list: List[np.ndarray] = []
        for i in range(self.nbDoFsPerNode):
            bfac = adv_corr[i] * g
            H = np.add.outer(adv_corr[i] * Dt, space_eigvals[i])
            S = compute_special_schur(H, bfac)
            bfac_list.append(bfac)
            H_list.append(H)
            schur_list.append(S)
        return H_list, bfac_list, schur_list

    def add_scalar_time_correctors(
        self,
        advection_corrector: Optional[List[float]] = [],
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

        if advection_corrector is not None and is_verified(advection_corrector):
            self._advection_time_corrector = advection_corrector
            self._blocks_time = None

    def apply_spacetime_preconditioner(self, array_in: np.ndarray) -> np.ndarray:
        """
        Apply the full space-time fast-diagonalization preconditioner to the input array.

        Args:
            array_in (np.ndarray): Input array to which the full space-time preconditioner will be applied
        """
        start = time()
        array_in = np.reshape(array_in, (self.nbDoFsPerNode, -1))
        array_out = np.zeros_like(array_in)
        H, bfac, S = self.blocks_time

        for i in range(self.nbDoFsPerNode):
            array1 = MatrixFree.apply(
                self.space_preconditioner.eigenvec_by_dir_space[i]
                + [self.eigenvec_time.conj()],
                array_in[i, self.sptm_free_nodes[i]],
                is_transpose=True,
            )

            # Solve arrowhead
            array2 = solve_special_arrowhead(
                H=H[i], bfac=bfac[i], rhs=np.atleast_2d(array1), schur=S[i]
            )

            array_out[i, self.sptm_free_nodes[i]] = np.real(
                MatrixFree.apply(
                    self.space_preconditioner.eigenvec_by_dir_space[i]
                    + [self.eigenvec_time],
                    np.ravel(array2),
                    is_transpose=False,
                )
            )

        logger.debug(
            f"Single patch fast-diagonalization in {time() - start:.2e} seconds"
        )
        return np.ravel(array_out)
