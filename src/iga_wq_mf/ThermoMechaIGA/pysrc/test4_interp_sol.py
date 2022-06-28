"""
.. Test of matrix free / interpolation
.. We test is Matrix-Free fortran subroutines work correctly
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
import scipy 

# My libraries
from lib.physics import powden_rotring, temperature_rotring
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq
from lib.fortran_mf_iga import fortran_mf_iga

# Set global variables
DEGREE = 4
CUTS = 4 # It can not exceed 5 

# Create geometry using geomdl
modelGeo = create_geometry(DEGREE, CUTS, 'RQA')

# ===========================================
# IGA WQ MF APPROACH
# ===========================================

for fortran_model in [fortran_mf_wq, fortran_mf_iga]:

    # Creation of thermal model object
    modelPhy = fortran_model(modelGeo)
    Temp_CP, Temp_Surface = modelPhy.interpolate_ControlPoints(temperature_rotring)
    _, qp_sample, _, Temp_Sample_interp = modelPhy.interpolate_field(u_ctrlpts=Temp_CP)
    Temp_Sample_exact = temperature_rotring(qp_sample)
    error_1 = np.linalg.norm(Temp_Sample_exact-Temp_Sample_interp, np.inf)/np.linalg.norm(Temp_Sample_exact, np.inf)*100

    # Block boundaries
    dof = modelPhy._thermal_dof
    dod = modelPhy._thermal_dod

    # Assemble conductivity matrix K
    K2nn = modelPhy.eval_conductivity_matrix(dof, dof)

    # Assemble source vector F
    F2n = modelPhy.eval_source_vector(powden_rotring, dof, indj=dod, Td=Temp_Surface)

    # Solve system
    Tn = scipy.linalg.solve(K2nn.todense(), F2n)

    # Assembly
    Tsolution = np.zeros(modelPhy._nb_ctrlpts_total)
    Tsolution[dof] = Tn
    Tsolution[dod] = Temp_Surface

    # Compare solutions 
    _, qp_sample, _, Temp_Sample_interp2 = modelPhy.interpolate_field(u_ctrlpts=Tsolution)
    error_2 = np.linalg.norm(Temp_Sample_interp-Temp_Sample_interp2, np.inf)/np.linalg.norm(Temp_Sample_interp, np.inf)*100

    print("Error using interpolation : %.3e %%" %(error_1,))
    print("Error interpolation/direct solution : %.3e %%" %(error_2,))

    # Solution using conjugate gradient
    iterations = 80
    epsilon = 1e-15

    # With preconditioner
    method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]
    for name in method_list: 
        inputs = [F2n, iterations, epsilon, name, Tn, True]   
        Tn_t, residue_t, error_t = modelPhy.MFsolver(*inputs)
        Tsolution_t = np.zeros(modelPhy._nb_ctrlpts_total)
        Tsolution_t[dof] = Tn_t
        Tsolution_t[dod] = Temp_Surface
        error_3 = np.linalg.norm(Tsolution-Tsolution_t, np.inf)/np.linalg.norm(Tsolution, np.inf)*100
        print("Error direct/iterative solution : %.3e %%" %(error_3,))
    
    print("------------------")