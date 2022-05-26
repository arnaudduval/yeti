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
from lib import blockPrint, enablePrint
from lib.base_functions import erase_rows_csr
from lib.fortran_mf_wq import wq_find_basis_weights_fortran
from iga_wq_mf import solver
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.methods_mf import MF
from lib.physics import powden_thickring

# # Set global variables
# DEGREE = 5
# NBEL = 128
# NB_CTRLPTS = DEGREE + NBEL

# start = time.time()
# # Define basis (the same for all directions)
# _, qp_wq, dB0, dB1, dW00, dW01, \
# dW10, dW11, indi, indj = wq_find_basis_weights_fortran(DEGREE, NBEL)
# stop = time.time()
# print('Time computing basis: %.3e s' %(stop-start,))

# # Erase data
# rows2erase = [0, -1]
# indi_t, indj_t, data_t = erase_rows_csr(rows2erase, indi, indj, 
#                         [dB0, dB1, dW00, dW01, dW10, dW11])
# [dB0_t, dB1_t, dW00_t, _, _, dW11_t] = data_t

# # Initialize
# shape_matrices, indexes, data = [], [], []
# for dim in range(3):
#     shape_matrices.append(len(qp_wq))
#     indexes.append(indi_t); indexes.append(indj_t) 
#     data.append(dB0_t); data.append(dB1_t)
#     data.append(dW00_t); data.append(dW11_t)

# # Assemble r vector
# start = time.time()
# len_dof = (NB_CTRLPTS-2)**3
# r = np.random.uniform(low=0, high=0.1, size=(len_dof,))
# stop = time.time()
# print('Time computing random vector: %.3e s' %(stop-start,))

# # Solve sylvester equation P s = r
# start = time.time()
# inputs = [*shape_matrices, *indexes, *data, r]
# s = solver.preconditioner_fd(*inputs)
# stop = time.time()
# print('Time computing Fast Diag: %.3e s' %(stop-start,))

# # =============================
# blockPrint()
# # Create geometry
# filename = 'thick_ring'

# DEGREE = 3
# CUTS = 4
# NBEL = 2**CUTS
# NB_CTRLPTS = DEGREE + NBEL

# # Create vector
# len_dof = (NB_CTRLPTS-2)**3
# r = np.ones(len_dof)
# rext = np.zeros(NB_CTRLPTS**3)

# geometry = {'degree': [DEGREE, DEGREE, DEGREE]}
# modelGeo = geomdlModel(filename=filename, **geometry)
# modelGeo.knot_refinement(nb_refinementByDirection= CUTS*np.array([1, 1, 1]))

# # Creation of thermal model object
# Model1 = fortran_mf_wq(modelGeo)
# dof = Model1._thermal_dof
# rext[dof] = r

# Ku1 = Model1.eval_Ku(r, table=Model1._thermalblockedboundaries)
# Cu1 = Model1.eval_Cu(r, table=Model1._thermalblockedboundaries)

# Model2 = MF(modelGeo)
# Ku2 = Model2.eval_conductivity_matrix() @ rext
# Ku2 = Ku2[dof]
# Ku3 = Model2.eval_Ku(rext)[dof]

# Cu2 = Model2.eval_capacity_matrix() @ rext
# Cu2 = Cu2[dof]
# Cu3 = Model2.eval_Cu(rext)[dof]

# # -------------
# enablePrint()
# error = np.abs(Cu3-Cu2)
# print(np.min(error), np.max(error))

# error = np.abs(Cu1-Cu2)
# print(np.min(error), np.max(error))

# print('-------------------')
# error = np.abs(Ku3-Ku2)
# print(np.min(error), np.max(error))

# error = np.abs(Ku1-Ku2)
# print(np.min(error), np.max(error))

# ==========================================
# Create geometry
filename = 'thick_ring'

DEGREE = 4
CUTS = 3
NBEL = 2**CUTS
NB_CTRLPTS = DEGREE + NBEL

# Create vector
len_dof = (NB_CTRLPTS-2)**3
r = np.ones(len_dof)
rext = np.zeros(NB_CTRLPTS**3)

geometry = {'degree': [DEGREE, DEGREE, DEGREE]}
modelGeo = geomdlModel(filename=filename, **geometry)
modelGeo.knot_refinement(nb_refinementByDirection= CUTS*np.array([1, 1, 1]))

# Creation of thermal model object
Model1 = fortran_mf_wq(modelGeo)
Model2 = MF(modelGeo)

# Block boundaries
dof = Model1._thermal_dof
sol_direct = None

# Recursive solver 
# ----------------
# Assemble source vector F
F2solve = Model1.eval_source_vector(powden_thickring, indi=dof)

if sol_direct is None: 
    sol_direct = np.ones(len(F2solve))
    print("Direct solution unknown. Default: ones chosen. Be aware of residue results")

# With and without preconditioner
epsilon  = 1e-10
iterations = 5

for name in ['TDS']:
    # With fortran
    _, residue_f, error_f = Model1.mf_conj_grad(F2solve, iterations, epsilon, name, sol_direct, True)
    print(residue_f)

    # With python
    _, residue_p = Model2.mf_wq_conj_grad(F2solve, dof, iterations, epsilon)
    print(residue_p)