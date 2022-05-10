# Python libraries
import re
import numpy as np
import scipy 

# My libraries
from lib.geomdl_geometry import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq
from lib.fortran_mf_iga import fortran_mf_iga

def temperature(dim, P):
    "T = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z)"
    x = P[0]
    y = P[1]
    z = P[2]
    u = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*np.sin(np.pi*z)
    return u

def powerdensity(dim, P):
    """
    u = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z)
    f = -(d2u/dx2 + d2u/dy2 + d2u/dz2)
    """
    x = P[0]
    y = P[1]
    z = P[2]
    f = (8*x*y**4*np.sin(np.pi*z) + 8*x**3*y**2*np.sin(np.pi*z) 
        + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 1) + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
        + 2*x*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
        - x*y**2*np.pi**2*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
    )
    return f

# Set global variables
DEGREE = 4
CUTS = 3

# Create geometry using geomdl
modelGeo = create_geometry(DEGREE, CUTS, 3)

# Creation of thermal model object
Model1 = fortran_mf_iga(modelGeo, thermalblockedboundaries=[[1,1], [1,1], [1,1]])
Temp_CP, Td = Model1.interpolate_MSE_CP(temperature, dirichlet0= [[1,0], [0,0], [1,0]])
_, qp_PS_1, _, Temp_qp_PS_1 = Model1.interpolate_results(u_ctrlpts=Temp_CP)
Treal_sample = np.asarray([temperature(3, qp_PS_1[:, :, _][0]) for _ in range(np.shape(qp_PS_1)[2])])
error_Temp = np.linalg.norm(Treal_sample-Temp_qp_PS_1, np.inf)/np.linalg.norm(Treal_sample, np.inf)
print(error_Temp)

# # ===========================================
# # WQ MF APPROACH
# # ===========================================
# Model1 = fortran_mf_wq(modelGeo, thermalblockedboundaries=[[1,1], [1,1], [1,1]])

# # Block boundaries
# dof = Model1._thermal_dof
# dod = Model1._thermal_dod

# # Assemble conductivity matrix K
# K2nn = Model1.eval_conductivity_matrix().tocsc()[dof, :][:, dof]
# K2nd = Model1.eval_conductivity_matrix().tocsc()[dof, :][:, dod]

# # Assemble source vector F
# F2n_v1 = np.asarray(Model1.eval_source_vector(powerdensity))[dof] - K2nd @ Td
# F2n = F2n_v2 = np.asarray(Model1.eval_source_vector(powerdensity, dod=dod, Td=Td))[dof]
# error_F = np.linalg.norm(F2n_v1-F2n_v2, np.inf)/np.linalg.norm(F2n_v2, np.inf)

# # Solve system
# Tn = scipy.linalg.solve(K2nn.todense(), F2n)

# # Assembly
# Tsolution = np.zeros(Model1._nb_ctrlpts_total)
# Tsolution[dof] = Tn
# Tsolution[dod] = Td

# # Compare solutions 
# _, qp_PS_2, _, Temp_qp_PS_2 = Model1.interpolate_results(u_ctrlpts=Tsolution)
# error_Temp = np.linalg.norm(Temp_qp_PS_2-Temp_qp_PS_1, np.inf)/np.linalg.norm(Temp_qp_PS_1, np.inf)
# print(error_Temp)

# # Solution using conjugate gradient
# iterations = 100
# epsilon = 1e-15

# # # With preconditioner
# # method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]
# # for name in method_list:
# #     inputs = [F2n, dof, iterations, epsilon, name, Tn, True]   
# #     Tn_t, residue_t, error_t = Model1.mf_conj_grad(*inputs)
 