"""
.. Test of matrix free / interpolation
.. We test is Matrix-Free fortran subroutines work correctly
.. Joaquin Cornejo 
"""

# Python libraries
import os, sys
import numpy as np
import scipy 

# My libraries
from lib.physics import powden_rotring, temperature
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq
from lib.fortran_mf_iga import fortran_mf_iga

# Enable and disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

# Set global variables
DEGREE = 4
CUTS = 3

# Create geometry using geomdl
modelGeo = create_geometry(DEGREE, CUTS, 3)

# ===========================================
# IGA WQ MF APPROACH
# ===========================================

for fortran_model in [fortran_mf_iga, fortran_mf_wq]:

    blockPrint()
    # Creation of thermal model object
    Model1 = fortran_model(modelGeo, thermalblockedboundaries=[[1,1], [1,1], [1,1]])
    Temp_CP, Td = Model1.MSE_ControlPoints(temperature)
    _, qp_PS_1, _, Temp_qp_PS_1 = Model1.interpolate_results(u_ctrlpts=Temp_CP)
    Treal_sample = np.asarray([temperature(3, qp_PS_1[:, :, _][0]) for _ in range(np.shape(qp_PS_1)[2])])
    error_Temp_1 = np.linalg.norm(Treal_sample-Temp_qp_PS_1, np.inf)/np.linalg.norm(Treal_sample, np.inf)

    # Block boundaries
    dof = Model1._thermal_dof
    dod = Model1._thermal_dod

    # Assemble conductivity matrix K
    K2nn = Model1.eval_conductivity_matrix(dof, dof)

    # Assemble source vector F
    F2n = Model1.eval_source_vector(powden_rotring, dof, indj=dod, Td=Td)

    # Solve system
    Tn = scipy.linalg.solve(K2nn.todense(), F2n)

    # Assembly
    Tsolution = np.zeros(Model1._nb_ctrlpts_total)
    Tsolution[dof] = Tn
    Tsolution[dod] = Td

    # Compare solutions 
    _, qp_PS_2, _, Temp_qp_PS_2 = Model1.interpolate_results(u_ctrlpts=Tsolution)
    error_Temp_2 = np.linalg.norm(Temp_qp_PS_1 - Temp_qp_PS_2, np.inf)/np.linalg.norm(Temp_qp_PS_1, np.inf)
    enablePrint()

    print("Error using interpolation : %.5f" %(error_Temp_1,))
    print("Error interpolation/direct solution : %.5f" %(error_Temp_2,))

    # Solution using conjugate gradient
    iterations = 80
    epsilon = 1e-15

    # With preconditioner
    method_list = ["C", "TDS", "JM"]
    for name in method_list:
        inputs = [F2n, dof, iterations, epsilon, name, Tn, True]   
        Tn_t, residue_t, error_t = Model1.mf_conj_grad(*inputs)
        Tsolution_t = np.zeros(Model1._nb_ctrlpts_total)
        Tsolution_t[dof] = Tn_t
        Tsolution_t[dod] = Td
        error_precond = np.linalg.norm(Tsolution-Tsolution_t, np.inf)/np.linalg.norm(Tsolution, np.inf)
        print("Error direct/iterative solution : %.5f" %(error_precond,))
    
    print("------------------")