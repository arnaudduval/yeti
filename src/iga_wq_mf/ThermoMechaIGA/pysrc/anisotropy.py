"""
.. Testing anisotry
.. Joaquin Cornejo 
"""

# Python libraries
from cProfile import label
import os
import numpy as np
from scipy import sparse as sp
import matplotlib.pyplot as plt

# My libraries
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq
from lib.fortran_mf_wq import wq_find_basis_weights_fortran
from lib.base_functions import compute_jacobien_mean, generate_rand_positive_matrix

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
DEGREE = 5
CUTS = 6
NBEL = 2**CUTS
NB_CTRLPTS = (DEGREE+NBEL)
GEO = 'TR'

# ===========================================
# Test MF-WQ methods
# ===========================================
file_exist = True
filename = folder + 'jacobienMean' + GEO

if not file_exist:
    # Create geometry using geomdl
    modelGeo = create_geometry(DEGREE, CUTS, GEO)

    # Creation of thermal model object
    modelPhy = fortran_mf_wq(modelGeo, isThermal= False)

    # Initialize
    data_save = []
    nb_pts = 2
    step = 1
    threshold = min(NB_CTRLPTS*2, 100)
    while nb_pts < threshold:
        print(nb_pts)
        
        jacobien_PS, qp_PS, detJ, _ = modelPhy.interpolate_field(nb_pts=nb_pts)
        dimensions = compute_jacobien_mean(jacobien_PS)
        data_save.append([nb_pts, *dimensions])

        if nb_pts < 5: step = 1
        elif nb_pts < 50: step = 5
        else: step = 10
        nb_pts += step

    data_save = np.array(data_save)
    np.savetxt(filename, data_save)

else: 
    data_save = np.loadtxt(filename)
    nb_pts = data_save[:, 0]**3
    dim = np.shape(data_save)[1]-1
    plt.figure(1)
    for _ in range(dim):
        current_data = data_save[:, _+1]
        ref = current_data[-1]
        current_data = abs(current_data-ref)/ref*100
        plt.loglog(nb_pts[:-1], current_data[:-1], 'o--', label='L'+str(_+1))
    
    # Set properties
    plt.grid()
    plt.ylim([1e-10, 1e2])
    plt.xlabel("Sample size", fontsize= 16)
    plt.ylabel("Relative error (%)", fontsize= 16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='best')

    plt.tight_layout()
    plt.savefig(filename+'.png')

# # ===========================================
# # Verify diagonal of K
# # ===========================================

# # Define basis (the same for all directions)
# _, qp_wq, dB0, dB1, dW00, dW01, \
# dW10, dW11, indi, indj = wq_find_basis_weights_fortran(DEGREE, NBEL)
# nb_qp_wq = len(qp_wq)
# nb_ctrlpts = DEGREE + NBEL
# nb_ctrlpts_total = nb_ctrlpts**3
# nb_qp_wq_total = nb_qp_wq**3
# indi -= 1; indj -= 1

# # Create basis and weights from fortran
# B0 = sp.csr_matrix((dB0, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
# B1 = sp.csr_matrix((dB1, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
# W00 = sp.csr_matrix((dW00, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
# W01 = sp.csr_matrix((dW01, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
# W10 = sp.csr_matrix((dW10, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
# W11 = sp.csr_matrix((dW11, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))

# DB = [B0, B1]
# DW = [[W00, W01], [W10, W11]]

# conductivity = generate_rand_positive_matrix(3, nb_qp_wq_total)

# # Initialize conductivity matrix 
# K = sp.csr_matrix((nb_ctrlpts_total, nb_ctrlpts_total))
# for j in range(3):
#     beta = np.zeros(3, dtype = int); beta[j] = 1
#     Bt = 1 
#     for dim in range(3): 
#         bt = beta[dim]
#         Bt = sp.kron(DB[bt], Bt)

#     for i in range(3):
#         Cij = conductivity[i, j, :]
#         alpha = np.zeros(3, dtype = int); alpha[i] = 1
        
#         Wt = 1
#         for dim in range(3): 
#             at = alpha[dim]
#             bt = beta[dim]
#             Wt = sp.kron(DW[at][bt], Wt)  

#         # Evaluates Cij * W in each point
#         Wt = sp.csr_matrix.dot(Wt, sp.diags(Cij))

#         # Find K = W C B
#         K += sp.csr_matrix.dot(Wt.tocsr()[:,:], Bt.tocsr()[:,:].T)

# # Extract diagonal
# diag_K = K.diagonal()
# min_val = np.min(diag_K)
# min_pos = np.where(diag_K == min_val)[0]
# min_row = K.getrow(min_pos)
# print(min_row.todense())

# max_val = np.max(diag_K)
# max_pos = np.where(diag_K == max_val)[0]
# max_row = K.getrow(max_pos)
# print(max_row.todense())

# plt.figure(1)
# plt.imshow(K.todense(), interpolation='none',cmap='binary')
# plt.colorbar()
# plt.savefig('sparsity.png')