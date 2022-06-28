"""
.. Test of setup time
.. We test how much time it takes to compute K and C matrices
.. Joaquin Cornejo 
"""

# Python libraries
import os, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpltools import annotation

# My libraries
from lib.fortran_mf_wq import wq_find_basis_weights_fortran
from iga_wq_mf import assembly

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test6/'
if not os.path.isdir(folder): os.mkdir(folder)

dataExist = True
nbel = 40

# Set filename
filename_WQ = 'setup_time_WQ_'  + 'nbel_' + str(nbel) 
filename_WQ = folder + filename_WQ

# Define inputs
degree_list = np.arange(9, 11)

# Output
time_matrix = np.zeros((len(degree_list), 3))
time_matrix[:, 0] = degree_list

if not dataExist:
    for i, degree in enumerate(degree_list):
        # Define basis (the same for all directions)
        nnz_I, qp_wq, dB0, dB1, dW00, dW01, \
        dW10, dW11, indi, indj = wq_find_basis_weights_fortran(degree, nbel)

        # Some other inputs
        nb_rows = nb_cols = (degree+nbel)**3
        
        # ------------
        capacity_coefs = np.ones(len(qp_wq)**3)
        shape_matrices, indices, data_C, size_I = [], [], [], []
        for dim in range(3):
            shape_matrices.append(len(qp_wq))
            indices.append(indi); indices.append(indj) 
            data_C.append(dB0); data_C.append(dW00)
            size_I.append(nnz_I)

        inputs_C = [capacity_coefs, *shape_matrices, *indices, *data_C, *size_I]
        
        start = time.time()
        assembly.wq_get_capacity_3d(*inputs_C)
        stop = time.time()
        time_matrix[i, 1] = stop - start 
        del inputs_C, capacity_coefs

        # ------------
        conductivity_coefs = np.ones((3, 3, len(qp_wq)**3))
        data_K = []
        for dim in range(3):
            data_K.append(dW00); data_K.append(dW01)
            data_K.append(dW10); data_K.append(dW11)
            data_K.append(dB0); data_K.append(dB1)

        inputs_K = [conductivity_coefs, *shape_matrices, *indices, *data_K, *size_I]

        start = time.time()
        assembly.wq_get_conductivity_3d(*inputs_K)
        stop = time.time()
        time_matrix[i, 2] = stop - start 
        del inputs_K, conductivity_coefs

        print(time_matrix[i, :])
        # np.savetxt(filename_WQ, time_matrix)

else :

    if nbel == 40:
        # Load data
        setup_WQ = np.loadtxt(filename_WQ)
        litterature = pd.read_table(folder + 'lit_WQ.dat', sep='\t', names=['degree', 'time'])

        # Plot
        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
        ax1.loglog(litterature.degree, litterature.time, 'o--', label='Litterature')
        ax1.loglog(setup_WQ[:,0], setup_WQ[:,1], 'x--', label='My algo')

        # Slope
        slope, _ = np.polyfit(np.log10(setup_WQ[:,0]),np.log10(setup_WQ[:,1]), 1)
        slope = round(slope, 3)
        annotation.slope_marker((setup_WQ[2,0], setup_WQ[2,1]), slope, 
                                text_kwargs={'fontsize': 14},
                                poly_kwargs={'facecolor': (0.73, 0.8, 1)})

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

        # Plot
        fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12,4))
        ax1.loglog(litterature.degree, litterature.C, 'o--', label='Litterature')
        ax1.loglog(setup_WQ[:,0], setup_WQ[:,1], 'x--', label='My algo')

        ax2.loglog(litterature.degree, litterature.K, 'o--', label='Litterature')
        ax2.loglog(setup_WQ[:,0], setup_WQ[:,2], 'x--', label='My algo')


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

        


