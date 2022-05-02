""" 
.. Calculating the time to assembly stiffness matrix
.. Joaquin Cornejo
"""

# Python libraries
import numpy as np
import time
import os, sys
from datetime import datetime
from scipy import sparse as sp
from matplotlib import pyplot as plt

# Yeti libraries 
from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix

# My libraries
from pysrc.lib.fortran_mf_wq import fortran_mf_wq

full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) 

# Enable and disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

# Yeti algorithms
# ---------------
FILENAME='volumetricBeam'

tas_IGA = []
tas_WQ = []
for DEGREE in range(2, 3):
    tas_IGA_t = []
    tas_WQ_t = []
    for CUTS in range(3, 4):
        time_IGA = 0.0
        time_WQ = 0.0

        # Print current time
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Current Time =", current_time)
        print([DEGREE, 2**CUTS])

        blockPrint()
        modeleIGA = IGAparametrization(filename=FILENAME)
        modeleIGA.refine(nb_degreeElevationByDirection=(DEGREE-1)*np.array([1, 1, 1]),
                            nb_refinementByDirection=CUTS*np.array([1, 1, 1]))

        start = time.time()
        data, row, col, Fb = build_stiffmatrix(*modeleIGA.get_inputs4system_elemStorage())
        Kside = sp.coo_matrix((data, (row,col)), shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                        dtype='float64').tocsc()
        S_python = Kside+Kside.transpose()
        plt.figure(1)
        plt.spy(S_python)
        plt.savefig(folder + '/StiffPython'+ '.png')
        plt.close(1)
        stop = time.time()
        time_IGA = stop - start

        # My algorithms 
        Model1 = fortran_mf_wq(modeleIGA, isThermal= False, isMechanical= True)
        E = Model1._youngModule
        nu = Model1._poissonCoef
        start = time.time()
        S_fortran = Model1.eval_stiffness_matrix()
        S_fortran *= E/((1 + nu)*(1 - 2*nu))
        plt.figure(1)
        plt.spy(S_fortran)
        plt.savefig(folder + '/StiffFortran'+ '.png')
        stop = time.time()
        time_WQ = stop - start
        
        del Model1, modeleIGA

        enablePrint()
        # Append information per nbel
        tas_IGA_t.append(time_IGA)
        tas_WQ_t.append(time_WQ)

#     # Append information per degree
#     tas_IGA.append(tas_IGA_t)
#     tas_WQ.append(tas_WQ_t)

# # np.savetxt('IGA_assembly.txt', tas_IGA, delimiter=',')
# # np.savetxt('WQ_assembly.txt', tas_WQ, delimiter=',')

# # ======================================
# WQ_time = np.loadtxt(folder + '/WQ_assembly.txt', delimiter=',')
# IGA_time = np.loadtxt(folder + '/IGA_assembly.txt', delimiter=',')

# m, n = np.shape(WQ_time)
# IGA_WQ = np.zeros((m, n))
# for i in range(m):
#     for j in range(n):
#         IGA_WQ[i, j] = IGA_time[i,j]/WQ_time[i, j]

# CUTS = np.arange(n)
# NBEL = [2**(cut+1) for cut in CUTS]
# plt.figure()
# for i in range(m):
#     plt.plot(NBEL, IGA_WQ[i, :], label='p = '+str(i+2))
# # Properties
# plt.grid()
# plt.yscale("log")
# plt.xlabel('Number of elements', fontsize=16)
# plt.ylabel('Time IGA / time WQ', fontsize=16)
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.legend(prop={'size': 14})
# plt.tight_layout()
# plt.savefig(folder + '/Time_acceleration'+ '.png')