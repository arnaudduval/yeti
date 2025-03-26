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

DEGREE = 6
CUTS = 5

for GEOMETRY_CASE in ['CB', 'VB', 'TR', 'RQA']:

    # Create geometry using geomdl
    modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)

    # Creation of thermal model object
    Model1 = fortran_mf_wq(modelGeo, isThermal=False)

    # Different type of plots and data
    Model1.export_results(filename= folder + GEOMETRY_CASE)
    