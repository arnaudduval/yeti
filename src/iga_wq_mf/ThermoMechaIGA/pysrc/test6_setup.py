"""
.. Test of setup time
.. We test how much time it takes to compute K and C matrices
.. Joaquin Cornejo 
"""

# Python libraries
import os, time
import numpy as np
from scipy import sparse as sp
import pandas as pd
import matplotlib.pyplot as plt

# My libraries
from lib.fortran_mf_wq import fortran_mf_wq, wq_find_basis_weights_fortran
from iga_wq_mf import assembly

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test6/'
if not os.path.isdir(folder): os.mkdir(folder)

dataExist = True
nbel = 64

# Set filename
filename_WQ = 'setup_time_WQ_'  + 'nbel_' + str(nbel) 
filename_WQ = folder + filename_WQ

# Define inputs
degree_list = np.arange(2, 7)

# Output
time_matrix = np.zeros((len(degree_list), 5))
time_matrix[:, 0] = degree_list

if not dataExist:
    for i, degree in enumerate(degree_list):
        # Define basis (the same for all directions)
        nnz_I, qp_wq, dB0, dB1, dW00, dW01, \
        dW10, dW11, indi, indj = wq_find_basis_weights_fortran(degree, nbel)

        # Some other inputs
        nb_rows = nb_cols = (degree+nbel)**3
        capacity_coefs = np.ones(len(qp_wq)**3)
        conductivity_coefs = np.ones((3, 3, len(qp_wq)**3))

        # Initialize
        shape_matrices, indices, data_C, data_K, size_I = [], [], [], [], []
        for dim in range(3):
            shape_matrices.append(len(qp_wq))
            indices.append(indi); indices.append(indj) 
            data_C.append(dB0); data_C.append(dW00)
            data_K.append(dB0); data_K.append(dB1)
            data_K.append(dW00); data_K.append(dW01)
            data_K.append(dW10); data_K.append(dW11)
            size_I.append(nnz_I)

        inputs_C = [capacity_coefs, *shape_matrices, *indices, *data_C, *size_I]
        inputs_K = [conductivity_coefs, *shape_matrices, *indices, *data_K, *size_I]
        
        start = time.time()
        val_C, indi_C, indj_C = assembly.wq_get_capacity_3d(*inputs_C)
        stop1 = time.time()
        CC = sp.csr_matrix((val_C, indj_C, indi_C), shape=(nb_rows, nb_cols))
        stop2 = time.time()
        time_matrix[i, 1] = stop1 - start 
        time_matrix[i, 2] = stop2 - start 

        start = time.time()
        val_K, indi_K, indj_K = assembly.wq_get_conductivity_3d(*inputs_K)
        stop1 = time.time()
        KK = sp.csr_matrix((val_K, indj_K, indi_K), shape=(nb_rows, nb_cols))
        stop2 = time.time()
        stop = time.time()
        time_matrix[i, 3] = stop1 - start 
        time_matrix[i, 4] = stop2 - start 

        del CC, KK

        print(time_matrix[i, :])
        np.savetxt(filename_WQ, time_matrix)

else :

    if nbel == 40:
        # Load data
        setup_WQ = np.loadtxt(filename_WQ)
        litterature = pd.read_table(folder + 'lit_WQ.dat', sep='\t', names=['degree', 'time'])

        # plot
        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
        ax1.loglog(litterature.degree, litterature.time, 'o--', label='Litterature')
        ax1.loglog(setup_WQ[:,0], setup_WQ[:,1], 'x--', label='My algo')

        # Plot reference
        c1, c2 = 0.125, 0.125
        xx = litterature.degree
        yy = c1*xx**3
        ax1.loglog(xx, yy, '-', label=r'$O(C p^3)$')

        # Properties
        for ax in [ax1]:
            ax.grid()
            ax.legend()
            ax.set_xlim([2, 10])
            ax.set_ylim([1e-1, 1e3])
            ax.set_ylabel('CPU setup time (s)', fontsize=16)
            ax.set_xlabel('Polynomial degree p', fontsize=16)
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)

        plt.tight_layout()
        plt.savefig(filename_WQ)

    if nbel == 64:
        # Load data
        setup_WQ = np.loadtxt(filename_WQ)
        litterature = pd.read_table(folder + 'lit_MF.dat', sep='\t', names=['degree', 'C', 'K'])

        # plot
        fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12,4))
        ax1.loglog(litterature.degree, litterature.C, 'o--', label='Litterature')
        ax1.loglog(setup_WQ[:,0], setup_WQ[:,1], 'x--', label='My algo')
        ax2.loglog(litterature.degree, litterature.K, 'o--', label='Litterature')
        ax2.loglog(setup_WQ[:,0], setup_WQ[:,3], 'x--', label='My algo')

        # Plot reference
        c1, c2 = 0.125, 0.125
        xx = litterature.degree
        yy = c1*xx**3
        yy2 = c2*xx**3
        ax1.loglog(xx, yy, '-', label=r'$O(C p^3)$')
        ax2.loglog(xx, yy2, '-', label=r'$O(C p^3)$')

        # Properties
        for ax in [ax1, ax2]:
            ax.grid()
            ax.legend()
            ax.set_xlim([2, 10])
            ax.set_ylim([1, 1e4])
            ax.set_ylabel('CPU setup time (s)', fontsize=16)
            ax.set_xlabel('Polynomial degree p', fontsize=16)
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)

        plt.tight_layout()
        plt.savefig(filename_WQ)

        


