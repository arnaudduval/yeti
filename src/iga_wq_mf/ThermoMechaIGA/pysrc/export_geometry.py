"""
.. Test geometry, plot determinant and export results
.. We test if geomdl works and if we can exprot the results in VTK format
.. Joaquin Cornejo 
"""

# My libraries
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq

degree = 3
cuts = 4

for geometry_case in ['SQ', 'VB', 'TR', 'RQA', 'CB']:

    # Create geometry using geomdl
    modelGeo = create_geometry(degree, cuts, geometry_case)

    # Creation of thermal model object
    modelPhy = fortran_mf_wq(modelGeo, isThermal=False)

    # Different type of plots and data
    modelPhy.export_results(filename= geometry_case)
    