"""
.. Test geometry, plot determinant and export results
.. We test if geomdl works and if we can exprot the results in VTK format
.. Joaquin Cornejo 
"""

# Python libraries
import os

# My libraries
from preprocessing.igaparametrization import IGAparametrization
from lib.create_geomdl import create_geometry, geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'

DEGREE = 3
CUTS = 5

# for GEOMETRY_CASE in ['CB', 'VB', 'TR', 'RQA']:
for GEOMETRY_CASE in ['TR']:

    if GEOMETRY_CASE == 'CB': filename = 'parallelepiped'
    elif GEOMETRY_CASE == 'VB': filename = 'prism'
    elif GEOMETRY_CASE == 'TR': filename = 'thick_ring'
    elif GEOMETRY_CASE == 'RQA': filename = 'rotated_quarter_annulus'
    elif GEOMETRY_CASE == 'SQ': filename = 'quadrilateral'

    # Create and refine model
    geometry = {'degree': [DEGREE, DEGREE, DEGREE]}
    modelGeo = geomdlModel(filename=filename, **geometry)
    modelGeo.write_YETI_inputfile(filename= folder + GEOMETRY_CASE)
    modelIGA = IGAparametrization(filename= folder + GEOMETRY_CASE)

    # # ==============================
    # # Create geometry using geomdl
    # modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)

    # # Creation of thermal model object
    # Model1 = fortran_mf_wq(modelGeo, isThermal=False)

    # # Different type of plots and data
    # Model1.export_results(filename= folder + GEOMETRY_CASE)
    