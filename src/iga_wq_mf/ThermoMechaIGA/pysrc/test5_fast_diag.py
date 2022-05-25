"""
.. Test of fast diagonalization
.. We verify speed of fast diagonalization as direct method to solve Sylvester system
.. Since Create Geomdl is slow, we avoid it and we define all necessary object to do Fast Diag.
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
import time

# My libraries
from lib.base_functions import erase_rows_csr
from lib.fortran_mf_wq import wq_find_basis_weights_fortran
from iga_wq_mf import solver

# Set global variables
filename = 'parallelepiped'
DEGREE = 3
CUTS = 7
NBEL = 2**CUTS
NB_CTRLPTS = DEGREE + NBEL

start = time.time()
# Define basis (the same for all directions)
_, qp_wq, dB0, dB1, dW00, dW01, \
dW10, dW11, indi, indj = wq_find_basis_weights_fortran(DEGREE, NBEL)
stop = time.time()
print('Time computing basis: %.3e s' %(stop-start,))

# Erase data
rows2erase = [0, -1]
indi_t, indj_t, data_t = erase_rows_csr(rows2erase, indi, indj, 
                        [dB0, dB1, dW00, dW01, dW10, dW11])
[dB0_t, dB1_t, dW00_t, _, _, dW11_t] = data_t

# Initialize
shape_matrices, indexes, data = [], [], []
for dim in range(3):
    shape_matrices.append(len(qp_wq))
    indexes.append(indi_t); indexes.append(indj_t) 
    data.append(dB0_t); data.append(dB1_t)
    data.append(dW00_t); data.append(dW11_t)

# Assemble r vector
start = time.time()
len_dof = (NB_CTRLPTS-2)**3
r = np.random.uniform(low=0, high=0.1, size=(len_dof,))
stop = time.time()
print('Time computing random vector: %.3e s' %(stop-start,))

# Solve sylvester equation P s = r
start = time.time()
inputs = [*shape_matrices, *indexes, *data, r]
s = solver.preconditioner_fd(*inputs)
stop = time.time()
print('Time computing Fast Diag: %.3e s' %(stop-start,))