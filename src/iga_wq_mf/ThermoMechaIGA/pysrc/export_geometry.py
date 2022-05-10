"""
.. Test geometry, plot determinant and export results
.. We test if geomdl works and if we can exprot the results in VTK format
.. Joaquin Cornejo 
"""

# Python libraries
import os

# My libraries
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'

DEGREE = 4
CUTS = 3

for GEOMETRY_CASE in range(4):

    if GEOMETRY_CASE == 0: txtname = 'CB' 
    elif GEOMETRY_CASE == 1: txtname = 'VB' 
    elif GEOMETRY_CASE == 2: txtname = 'TR' 
    elif GEOMETRY_CASE == 3: txtname = 'RQA' 

    # Create geometry using geomdl
    modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)

    # Creation of thermal model object
    Model1 = fortran_mf_wq(modelGeo)

    # Different type of plots and data
    Model1.export_results(filename= folder + txtname)
