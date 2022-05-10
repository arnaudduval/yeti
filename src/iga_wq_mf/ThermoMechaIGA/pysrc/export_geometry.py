"""
.. Test geometry, plot determinant and export results
.. We test if geomdl works and if we can exprot the results in VTK format
.. Joaquin Cornejo 
"""

# Python libraries
import os
import numpy as np

# My libraries
from lib.geomdl_geometry import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.geomdl_geometry import create_geometry

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'

DEGREE = 4
CUTS = 3
# GEOMETRY_CASE = 3

# if GEOMETRY_CASE == 0: txtname = 'CB' 
# elif GEOMETRY_CASE == 1: txtname = 'VB' 
# elif GEOMETRY_CASE == 2: txtname = 'TR' 
# elif GEOMETRY_CASE == 3: txtname = 'RQA' 

# # Create geometry using geomdl
# modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)

filename = 'quarter_annulus'
txtname = 'QA2D'
geometry = {'degree': [DEGREE, DEGREE]}
modelGeo = geomdlModel(filename=filename, **geometry)
modelGeo.knot_refinement(nb_refinementByDirection= CUTS*np.array([1, 1, 1]))

# Creation of thermal model object
Model1 = fortran_mf_wq(modelGeo)

# Different type of plots and data
Model1.export_results(filename= folder + txtname)
