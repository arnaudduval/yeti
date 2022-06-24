"""
.. Testing anisotry
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
import scipy 

# My libraries
from lib.physics import powden_annulus
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq

# Set global variables
DEGREE = 4
CUTS = 4 

# Create geometry using geomdl
modelGeo = create_geometry(DEGREE, CUTS, 'QA')

# ===========================================
# IGA WQ MF APPROACH
# ===========================================

# Creation of thermal model object
# conductivity = np.eye(2)
conductivity = np.array([[1, 0],[0, 0.1]])
properties = {"conductivity": conductivity}
modelPhy = fortran_mf_wq(modelGeo, **properties)

# Block boundaries
dof = modelPhy._thermal_dof
dod = modelPhy._thermal_dod

# Assemble conductivity matrix K
K2nn = modelPhy.eval_conductivity_matrix(dof, dof)

# Assemble source vector F
F2n = modelPhy.eval_source_vector(powden_annulus, dof)

# Solve system
Tn = scipy.sparse.linalg.spsolve(K2nn, F2n)

# Assembly
Tsolution = np.zeros(modelPhy._nb_ctrlpts_total)
Tsolution[dof] = Tn

modelPhy.export_results(u_ctrlpts= Tsolution)