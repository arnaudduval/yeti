"""
.. Test of matrix free
.. We test is fortran subroutines work correctly
.. Joaquin Cornejo 
"""
# Python libraries
import numpy as np
from scipy import sparse as sp

# My libraries
from lib.geomdl_geometry import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq
from lib.physics import (powden_cube, 
                        powden_prism,
                        powden_thickring
)
from lib.methods_mf import MF

# Set global variables
DEGREE = 3
CUTS = 3
GEOMETRY_CASE = 0
IS_CG = True

if GEOMETRY_CASE == 0: txtname = 'C'; funpower = powden_cube
elif GEOMETRY_CASE == 1: txtname = 'VB'; funpower = powden_prism
elif GEOMETRY_CASE == 2: txtname = 'TR'; funpower = powden_thickring

# Create geometry using geomdl
modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)

# Creation of thermal model object
Model1 = fortran_mf_wq(modelGeo)

# Block boundaries
dof = Model1._thermal_dof

# Compute diagonal of K
# ---------------------
Kwq = Model1.eval_conductivity_matrix()
Kwq_diag = sp.csr_matrix.diagonal(Kwq)
Kmf_diag = Model1.mf_find_K_diagonal()

# Compare results 
error = Kwq_diag - Kmf_diag
try: norm = sp.linalg.norm(error, np.inf)/sp.linalg.norm(Kwq_diag, np.inf)
except: norm = np.linalg.norm(error, np.inf)/np.linalg.norm(Kwq_diag, np.inf)
print('Error calculating diagonal: %.05f' %(norm,))

# =======================================================
# Compute Ku 
# -------------
# Assemble conductivity matrix K
K2solve = Kwq.tocsc()[dof, :][:, dof]

# Assemble source vector F
F2solve = np.asarray(Model1.eval_source_vector(funpower))[dof]

# Solve system
LU = sp.linalg.splu(K2solve)
sol_d  = LU.solve(F2solve)

# Set parameters
epsilon  = 1e-14 
iterations = 4
inputs = [F2solve, dof, iterations, epsilon, "WP", sol_d, IS_CG]   
sol_t1, residue_t1, error_t4 = Model1.mf_conj_grad(*inputs)

# Compare with MF solver
Model1 = MF(modelGeo)
sol_t2, residue_t2 = Model1.mf_wq_conj_grad(F2solve, dof, iterations, epsilon)
error = sol_t1 - sol_t2 
norm = np.linalg.norm(error, np.inf)/np.linalg.norm(sol_t2, np.inf)
print('Error calculating CG: %.05f' %(norm,))








