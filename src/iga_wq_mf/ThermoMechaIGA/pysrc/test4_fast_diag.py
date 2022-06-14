"""
.. Test of fast diagonalization
.. We verify speed of fast diagonalization as direct method to solve Sylvester system
.. Since Create Geomdl is slow, we avoid it and we define all necessary object to do Fast Diag.
.. Joaquin Cornejo 
"""

# Python libraries
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from mpltools import annotation
import time

# My libraries
from lib import blockPrint, enablePrint
from lib.base_functions import erase_rows_csr
from lib.fortran_mf_wq import wq_find_basis_weights_fortran
from iga_wq_mf import solver

dataExist = False

if not dataExist:
    # Set global variables
    DEGREE = 5
    for NBEL in range(300, 601, 100): 
        NB_CTRLPTS = DEGREE + NBEL
        # print('********')
        # print(NB_CTRLPTS)

        start = time.time()
        # Define basis (the same for all directions)
        _, qp_wq, dB0, dB1, dW00, dW01, \
        dW10, dW11, indi, indj = wq_find_basis_weights_fortran(DEGREE, NBEL)
        stop = time.time()
        # print('Time computing basis: %.3e s' %(stop-start,))

        # Erase data
        rows2erase = [0, -1]
        indi_t, indj_t, data_t = erase_rows_csr(rows2erase, indi, indj, 
                                [dB0, dB1, dW00, dW01, dW10, dW11])
        [dB0_t, dB1_t, dW00_t, _, _, dW11_t] = data_t

        # Initialize
        shape_matrices, indices, data = [], [], []
        for dim in range(3):
            shape_matrices.append(len(qp_wq))
            indices.append(indi_t); indices.append(indj_t) 
            data.append(dB0_t); data.append(dB1_t)
            data.append(dW00_t); data.append(dW11_t)

        # Solve sylvester equation P s = r
        start = time.time()
        inputs = [*shape_matrices, *indices, *data]
        solver.test_precondfd(*inputs)
        stop = time.time()
        # print('Time computing Fast Diag: %.3e s' %(stop-start,))
        # time.sleep(1)

else:
    # PLOT 
    # Chosee folder
    full_path = os.path.realpath(__file__)
    folder_file = os.path.dirname(full_path) + '/data/'
    folder_figure = os.path.dirname(full_path) + '/results/'

    # Import data
    file1 = pd.read_table(folder_file + 'Algo1.dat', sep='\t', names=['nbel', 'time'])
    file2 = pd.read_table(folder_file + 'Algo2.dat', sep='\t', names=['nbel', 'time'])
    file3 = pd.read_table(folder_file + 'Algo3.dat', sep='\t', names=['nbel', 'time'])

    # Post treatment
    files = [file1, file2, file3]
    labels = ['New algorithm', 'Former algorithm', 'Litterature']
    for file, label in zip(files, labels):

        # Extract data
        nbel = file.nbel**3
        times = file.time
        nbdata = len(nbel)

        # Plot data
        plt.loglog(nbel, times, 'o--', label=label)

        # Get slope
        slope, _ = np.polyfit(np.log10(nbel),np.log10(times), 1)
        slope = round(slope, 3)
        annotation.slope_marker((nbel[round((nbdata-1)/2)], times[round((nbdata-1)/2)]), slope, 
                                text_kwargs={'fontsize': 14},
                                poly_kwargs={'facecolor': (0.73, 0.8, 1)})

    # Set properties
    plt.grid()
    plt.xlabel("Total DOF", fontsize= 16)
    plt.ylabel("CPU time (s)", fontsize= 16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='best')

    plt.tight_layout()
    plt.savefig(folder_figure + 'ComplexityAlgo' + '.png')
