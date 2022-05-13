"""
.. Test of matrix free "reduced"
.. We test is Matrix-Free fortran subroutines work correctly
.. Joaquin Cornejo 
"""

# Python libraries
import os, sys
import numpy as np
from scipy import sparse as sp

# My libraries
from lib.base_functions import erase_rows_csr
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq

# Enable and disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

# Set global variables
DEGREE = 4
CUTS = 3
NBEL = 2**CUTS

# Create geometry using geomdl
modelGeo = create_geometry(DEGREE, CUTS, 3)

# ===========================================
# IGA WQ MF APPROACH
# ===========================================
# Creation of thermal model object
Model1 = fortran_mf_wq(modelGeo, thermalblockedboundaries=[[1,1], [1,1], [1,1]])
indi, indj, data = erase_rows_csr([0, -1], *Model1._indexes[0], [*Model1._DB[0], *Model1._DW[0][0], *Model1._DW[0][1]])
[dB0, dB1, dW00, dW01, dW10, dW11] = data
nb_rows = np.size(indi)-1
nb_cols = Model1._nb_qp_wq[0][0]
B0 = sp.csr_matrix((dB0, indj, indi), shape=(nb_rows, nb_cols))
print(B0.toarray())