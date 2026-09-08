from yeti_iga.pymfiga.common.base.enum import BoundarySide, ParametricDirection, Constants
from yeti_iga.pymfiga.common.base.cls import MortarInterface
from yeti_iga.pymfiga.common.numerics.quadrature_rules import IGAQuadratureRule
from yeti_iga.pymfiga.common.numerics.operations import BsplineOperations, NurbsOperations
from yeti_iga.pymfiga.common.numerics.solvers import LinearSolver
from yeti_iga.pymfiga.iga.single_model.cls.cspace import SingleSpatialModel
from yeti_iga.pymfiga.iga.fastdiagonalization import SingleFD
from yeti_iga.pymfiga.iga.model_manager import ModelManager
from .helpers import kron_axis0
from typing import Dict, List, Union, Any, Tuple, Callable
from time import time
import scipy.sparse as sp
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.MULTI_MODEL")


class MortarCoupling:
    """
    Physically consistent dual mortar coupling.
    Dual basis constructed from physical slave Gram matrix.
    """

    def __init__(self, multipatch: ModelManager):
        self._multipatch = multipatch
        self._interfaces = self._translate_interfaces(multipatch.interfaces)
        self._M = self._assemble_global_coupling_matrix()

    # ------------------------------------------------------------------
    # Interface translation
    # ------------------------------------------------------------------

    def _translate_interfaces(
        self, interfaces: List[MortarInterface]
    ) -> Dict[tuple, Dict[str, Any]]:
        translated = {}
        for interface in interfaces:
            translated[(interface.master_patch, interface.slave_patch)] = {
                "master_boundary_id": interface.master_side,
                "slave_boundary_id": interface.slave_side,
            }
        return translated

    # ------------------------------------------------------------------
    # Assembly
    # ------------------------------------------------------------------

    def _assemble_global_coupling_matrix(self) -> sp.csr_array:

        rows = []
        cols = []
        data = []

        current_row = 0

        for (p_master, p_slave), info in self._interfaces.items():

            master = self._multipatch.patch_id_model_dict[p_master]
            slave = self._multipatch.patch_id_model_dict[p_slave]
            assert isinstance(master, SingleSpatialModel) and isinstance(
                slave, SingleSpatialModel
            )

            master_boundary_id: Tuple[ParametricDirection, BoundarySide, List[int]] = (
                info["master_boundary_id"]
            )
            slave_boundary_id: Tuple[ParametricDirection, BoundarySide, List[int]] = (
                info["slave_boundary_id"]
            )

            # Build physical Gram and dual
            M_operator, P_operator, slave_data = self._assemble_slave_gram(
                slave,
                slave_boundary_id,
                master.part.dir_degrees,
            )

            # Compute local master block
            M_master_local = self._assemble_master_block(
                p_master,
                master_boundary_id,
                M_operator,
                P_operator,
                slave_data,
            )

            master_offset = self._multipatch.patch_offsets[p_master]
            slave_offset = self._multipatch.patch_offsets[p_slave]

            # Insert master block
            M_coo = M_master_local.tocoo()
            for i, j, v in zip(M_coo.row, M_coo.col, M_coo.data):
                rows.append(current_row + i)
                cols.append(master_offset + j)
                data.append(v)

            # Insert slave block: Identity block due to biorthogonality
            slave_boundary_nodes = slave_boundary_id[2]

            for comp in range(self._multipatch.nbvars):
                for local_idx, node in enumerate(slave_boundary_nodes):

                    global_dof = (
                        slave_offset
                        + node
                        + comp * self._multipatch.num_nodes_patch[p_slave]
                    )

                    row_id = current_row + local_idx + comp * len(slave_boundary_nodes)

                    rows.append(row_id)
                    cols.append(global_dof)
                    data.append(-1.0)

            Mshape = M_master_local.shape
            assert Mshape is not None
            current_row += Mshape[0]

        M = sp.coo_array(
            (data, (rows, cols)),
            shape=(current_row, self._multipatch.get_size_of_arrays()),
        ).tocsr()
        M.eliminate_zeros()
        return sp.csr_array(M)

    # ------------------------------------------------------------------
    # Physical slave Gram matrix
    # ------------------------------------------------------------------

    def _assemble_slave_gram(
        self,
        slave: SingleSpatialModel,
        slave_boundary_id: Tuple[ParametricDirection, BoundarySide, List[int]],
        master_degrees: Dict[ParametricDirection, int],
    ) -> Tuple[Callable, Callable, tuple]:

        # Compute data on slave boundary
        slave_indices = slave_boundary_id[2]
        data_on_slave_boundary = slave.part.get_data_on_surface(
            slave_indices,
            slave_boundary_id[0],
            p_degrees=master_degrees,
        )
        _, jacobian, nurbs_weights, quadrules = data_on_slave_boundary

        # M_ij = ∫ N_i N_j dΓ (in physical space)
        def matvec(x):
            v = slave.operator_engine.compute_mf_scalar_u_v(
                quadrules,
                jacobian,
                x,
                allow_lumping=False,
                nurbs_weights=nurbs_weights,
            )
            return v

        # P_ij = ∫ hatN_i hatN_j dxi (in parametric space)
        # TODO: investigate if this preconditioner works
        # compare it with the preconditioner in eigs (Physics)
        fd = SingleFD()
        fd.compute_space_eigendecomposition(quadrules, np.zeros((1, len(quadrules), 2)))
        fd.update_space_eigenvalues(scalar_coefs=[1.0, 0.0])
        return matvec, fd.apply_spatial_preconditioner, data_on_slave_boundary

    # ------------------------------------------------------------------
    # Master block assembly using physical dual
    # ------------------------------------------------------------------

    def _assemble_master_block(
        self,
        p_master: Union[int, str],
        master_boundary_id: Tuple[ParametricDirection, BoundarySide, List[int]],
        mass_slave: Callable,
        preconditioner_slave: Callable,
        data_on_slave_boundary: Tuple[
            np.ndarray,
            np.ndarray,
            np.ndarray,
            List[IGAQuadratureRule],
        ],
    ) -> sp.csr_array:
        rows = []
        cols = []
        data = []

        # Unravel slave data
        assert len(data_on_slave_boundary) == 4
        (
            slv_phypts,
            slv_jacobian,
            slv_nurbs_weights,
            slv_quadrules,
        ) = data_on_slave_boundary

        # Compute projection of quadrature points on slave to master
        master = self._multipatch.patch_id_model_dict[p_master]
        assert isinstance(master, SingleSpatialModel)
        mst_indices = master_boundary_id[2]
        mst_quadpts = master.part.inverse_projection_on_boundary(
            slv_phypts,
            (master_boundary_id[0], master_boundary_id[1]),
        )  # the shape is (d, nbquadpts)

        # Compute master basis
        mst_directions = ParametricDirection.get_integration_dirs(
            master_boundary_id[0], master.ndim
        )
        mst_basis_dir: List[np.ndarray] = []
        for direction in mst_directions:
            q = master.part.quadrule_list[direction.value]
            q.knots_to_sample = mst_quadpts[direction.value]
            mst_basis_dir.append(q.basis_to_sample[0].toarray().copy().T)
            q.clear_sample()

        # NOTE: this lines is a variant of BsplineOperations.spkron_product_quadrature
        mst_basis = kron_axis0(mst_basis_dir)  # the shape is (N_master, nbquadpts)
        if master.part.nurbs_weights.size > 0:
            mst_nurbs_weights = master.part.nurbs_weights[mst_indices]
            mst_inv_weights_proj = 1.0 / (mst_basis.T @ mst_nurbs_weights)
            mst_basis = np.einsum(
                "i,ij,j->ij", mst_nurbs_weights, mst_basis, mst_inv_weights_proj
            )

        # We should compute [M_s]^-1 @ M_p, where M_s is the slave Gram matrix
        # and M_p_AB = int N_s_A N_m_B dΓ with N_s_A the slave basis
        # and N_m_B the master basis. However we are computing them iteratively.
        # mass_slave is M_s operator and preconditioner_slave is an approximation of M_s
        # M_p_AB is actually a set of "forces"
        n_master = mst_basis.shape[0]
        n_slave = np.prod([q.nbctrlpts for q in slv_quadrules])
        integrand = np.zeros((n_slave, n_master))
        operator_engine = (
            NurbsOperations if slv_nurbs_weights.size > 0 else BsplineOperations
        )
        linear_solver = LinearSolver(
            tolerance=Constants.TINY, maxiters=100, linear_type="cg", verbose=False
        )
        start = time()
        for j in range(n_master):
            # NOTE: it is safe to use atleast_2d and then ravel since we know that is a column
            rhs = operator_engine.assemble_scalar_u_force(
                slv_quadrules,
                np.atleast_2d(slv_jacobian * mst_basis[j]),
                slv_nurbs_weights,
            )
            sol = linear_solver.solve(
                mass_slave, np.ravel(rhs), Pfun=preconditioner_slave
            )["sol"]
            integrand[:, j] = sol
        logger.info(f"Compute mortar matrix in {time()-start:.2e} seconds")

        # Assemble master block
        assert n_master == len(mst_indices)
        for comp in range(self._multipatch.nbvars):
            for i in range(n_slave):
                for j in range(n_master):
                    row_id = i + comp * n_slave
                    col_id = (
                        mst_indices[j]
                        + comp * self._multipatch.num_nodes_patch[p_master]
                    )
                    rows.append(row_id)
                    cols.append(col_id)
                    data.append(integrand[i, j])

        M = sp.coo_array(
            (data, (rows, cols)),
            shape=(
                n_slave * self._multipatch.nbvars,
                self._multipatch.num_nodes_patch[p_master] * self._multipatch.nbvars,
            ),
        ).tocsr()
        M.eliminate_zeros()
        return sp.csr_array(M)

    # ------------------------------------------------------------------
    # Public
    # ------------------------------------------------------------------

    @property
    def coupling_matrix(self):
        if self._M is None:
            self._M = self._assemble_global_coupling_matrix()
        return self._M

    def build_constraint_matrix_vector(self) -> Tuple[sp.csr_array, np.ndarray]:
        Cmat = self.coupling_matrix
        gvec = np.zeros(Cmat.shape[0])
        return Cmat, gvec
