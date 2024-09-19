"""
.. Test of basis and weights 
.. We test if functions done in python and fortran for WQ approach 
.. gives the expected results.
.. Joaquin Cornejo 
"""

# Python libraries
import os
import numpy as np
from scipy import sparse as sp
import matplotlib.pyplot as plt

# My libraries
from lib.base_functions import (create_knotvector, 
                                eval_basis_python,
                                iga_find_positions_weights,
                                wq_find_basis_weights_opt
)
from lib.fortran_mf_wq import wq_find_basis_weights_fortran

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set number of elements
CUTS = np.arange(1, 6)
NBEL = [2**_ for _ in CUTS]

for varName in ['I00', 'I01', 'I10', 'I11']:
    plt.figure(1)
    ax = plt.gca()
    for degree in range(2, 6):

        norm_fortran = []; norm_python = []; ddl =[]
        color = next(ax._get_lines.prop_cycler)['color']

        for nb_el in NBEL: 

            # ========================================
            # FORTRAN
            # ======================================== 
            _, _, dB0, dB1, dW00, dW01, \
            dW10, dW11, indi, indj = wq_find_basis_weights_fortran(degree, nb_el)
            nb_ctrlpts = degree + nb_el
            nb_qp_wq = np.max(indj)
            indi -= 1; indj -= 1

            # Create basis and weights from fortran
            B0f = sp.csr_matrix((dB0, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            B1f = sp.csr_matrix((dB1, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            W00f = sp.csr_matrix((dW00, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            W01f = sp.csr_matrix((dW01, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            W10f = sp.csr_matrix((dW10, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))
            W11f = sp.csr_matrix((dW11, indj, indi), shape=(nb_ctrlpts, nb_qp_wq))

            # Calculate I
            I00f = W00f @ B0f.T
            I01f = W01f @ B1f.T
            I10f = W10f @ B0f.T
            I11f = W11f @ B1f.T

            # ========================================
            # PYTHON
            # ========================================
            knotvector = create_knotvector(degree, nb_el)
            _, B0p, B1p, W00p, W01p, W10p, W11p = wq_find_basis_weights_opt(degree, knotvector, 2)

            # Calculate I
            I00p = W00p @ B0p.T
            I01p = W01p @ B1p.T
            I10p = W10p @ B0p.T
            I11p = W11p @ B1p.T

            # ========================================
            # REFERENCE
            # ========================================
            qp_cgg, Wcgg = iga_find_positions_weights(degree, knotvector)
            B0, B1 = eval_basis_python(degree, knotvector, qp_cgg)

            # Calculate I
            I00 = B0 @ np.diag(Wcgg) @ B0.T
            I01 = B0 @ np.diag(Wcgg) @ B1.T
            I10 = B1 @ np.diag(Wcgg) @ B0.T
            I11 = B1 @ np.diag(Wcgg) @ B1.T

            # To choose variables
            if varName == 'I00': var1 = I00; var2 = I00f; var3 = I00p
            elif varName == 'I01': var1 = I01; var2 = I01f; var3 = I01p
            elif varName == 'I10': var1 = I10; var2 = I10f; var3 = I10p
            elif varName == 'I11': var1 = I11; var2 = I11f; var3 = I11p

            # Compare results 
            error_fortran = var1 - var2
            norm_temp = np.linalg.norm(error_fortran, np.inf)/np.linalg.norm(var1, np.inf)
            if norm_temp > 1e-5:
                raise Warning("Something happend. Fortran basis are wrong")
            norm_fortran.append(norm_temp)

            error_python = var1 - var3
            norm_temp = np.linalg.norm(error_python, np.inf)/np.linalg.norm(var1, np.inf)
            if norm_temp > 1e-5:
                raise Warning("Something happend. Python basis are wrong")
            norm_python.append(norm_temp)

        # Change type 
        norm_fortran = np.asarray(norm_fortran)
        norm_python = np.asarray(norm_python)
        ddl = np.asarray(NBEL)

        # Figure 
        plt.figure(1)
        plt.plot(ddl, norm_fortran*100, '-o', label='F p = ' + str(degree), color=color)
        plt.plot(ddl, norm_python*100, '--P', color=color)

    # Plot configurations
    plt.grid()
    plt.xscale("log")
    plt.xlabel("Number of elements $nb_{el}$", fontsize= 16)
    plt.yscale("log")
    plt.ylabel("Relative error (%)", fontsize= 16)
    plt.xlim([1, 100])
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(bbox_to_anchor= (1.05, 1.0), loc= 'upper left')
    plt.tight_layout()

    plt.savefig(folder + 'Error_basisweights_' + varName +'.png')
    plt.figure(1).clear()
