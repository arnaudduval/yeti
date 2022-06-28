"""
.. Testing anisotry
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
import scipy 
from scipy import sparse as sp

# My libraries
from lib.physics import powden_annulus
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq
from lib.fortran_mf_wq import wq_find_basis_weights_fortran
from lib.base_functions import create_knotvector, generate_rand_positive_matrix, wq_find_basis_weights_opt


# Set global variables
DEGREE = 7
CUTS = 3
NBEL = 2**CUTS

# # Create geometry using geomdl
# modelGeo = create_geometry(DEGREE, CUTS, 'QA')

# # ===========================================
# # IGA WQ MF APPROACH
# # ===========================================

# # Creation of thermal model object
# conductivity = np.array([[1, 0],[0, 1]])
# properties = {"conductivity": conductivity}
# modelPhy = fortran_mf_wq(modelGeo, **properties)

# # Block boundaries
# dof = modelPhy._thermal_dof
# dod = modelPhy._thermal_dod

# # Assemble conductivity matrix K
# K2nn = modelPhy.eval_conductivity_matrix(dof, dof)

# # Assemble source vector F
# F2n = modelPhy.eval_source_vector(powden_annulus, dof)

# # Solve system
# Tn = scipy.sparse.linalg.spsolve(K2nn, F2n)

# # Assembly
# Tsolution = np.zeros(modelPhy._nb_ctrlpts_total)
# Tsolution[dof] = Tn

# modelPhy.export_results(u_ctrlpts= Tsolution)

# =================================================

# Define basis (the same for all directions)
_, qp_wq, dB0, dB1, dW00, dW01, \
dW10, dW11, indi, indj = wq_find_basis_weights_fortran(DEGREE, NBEL)
nb_qp_wq = len(qp_wq)
nb_ctrlpts = DEGREE + NBEL
nb_ctrlpts_total = nb_ctrlpts**3
nb_qp_wq_total = nb_qp_wq**3
indi -= 1; indj -= 1

# Create basis and weights from fortran
B0 = sp.csr_matrix((dB0, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
B1 = sp.csr_matrix((dB1, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
W00 = sp.csr_matrix((dW00, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
W01 = sp.csr_matrix((dW01, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
W10 = sp.csr_matrix((dW10, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
W11 = sp.csr_matrix((dW11, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))

# knotvector = create_knotvector(DEGREE, NBEL)
# qp_wq, B0, B1, W00, W01, W10, W11 = wq_find_basis_weights_opt(DEGREE, knotvector, 2)
# nb_qp_wq = len(qp_wq)
# nb_ctrlpts = DEGREE + NBEL
# nb_ctrlpts_total = nb_ctrlpts**3
# nb_qp_wq_total = nb_qp_wq**3

DB = [B0, B1]
DW = [[W00, W01], [W10, W11]]

conductivity = generate_rand_positive_matrix(3, nb_qp_wq_total)

# Initialize conductivity matrix 
K = sp.csr_matrix((nb_ctrlpts_total, nb_ctrlpts_total))
for j in range(3):
    beta = np.zeros(3, dtype = int); beta[j] = 1
    Bt = 1 
    for dim in range(3): 
        bt = beta[dim]
        Bt = sp.kron(DB[bt], Bt)

    for i in range(3):
        Cij = conductivity[i, j, :]
        alpha = np.zeros(3, dtype = int); alpha[i] = 1
        
        Wt = 1
        for dim in range(3): 
            at = alpha[dim]
            bt = beta[dim]
            Wt = sp.kron(DW[at][bt], Wt)  

        # Evaluates Cij * W in each point
        Wt = sp.csr_matrix.dot(Wt, sp.diags(Cij))

        # Find K = W C B
        K += sp.csr_matrix.dot(Wt.tocsr()[:,:], Bt.tocsr()[:,:].T)

# Extract diagonal
diag_K = K.diagonal()
print(np.min(diag_K), np.max(diag_K))