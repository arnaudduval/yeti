"""
.. Test of matrix free / interpolation
.. We test is Matrix-Free fortran subroutines work correctly
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np, scipy

# My libraries
from lib.physics import powden_rotring, temperature_rotring
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.fortran_mf_iga import fortran_mf_iga

# Set global variables
isIGA = False
degree = 4
cuts = 3 # It can not exceed 5 

if isIGA: classfortran = fortran_mf_iga
else: classfortran = fortran_mf_wq

# Create geometry using geomdl
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('RQA', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# ===========================================
# IGA WQ MF APPROACH
# ===========================================
# Interpolation of u
modelPhy = classfortran(modelIGA)
u_interp = modelPhy.interpolate_ControlPoints(temperature_rotring)
_, qp_sample, _, u_interp_sample = modelPhy.interpolate_field(u_ctrlpts=u_interp)
u_exact_sample = temperature_rotring(qp_sample)
error_1 = np.linalg.norm(u_exact_sample-u_interp_sample, np.inf)/np.linalg.norm(u_exact_sample, np.inf)*100

# Add material 
material = {'capacity':1, 'conductivity':np.eye(3)}
modelPhy._set_material(material)

# Block boundaries
Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
modelPhy._set_dirichlet_boundaries(Dirichlet)
dof = modelPhy._thermal_dof
dod = modelPhy._thermal_dod 

# Solve system
ud = u_interp[dod]
K = modelPhy.eval_conductivity_matrix()
Knn = K[dof, :][:, dof]
Knd = K[dof, :][:, dod]
Fn = modelPhy.eval_source_vector(powden_rotring, dof) - Knd @ ud
un = scipy.linalg.solve(Knn.todense(), Fn)

usolved = np.zeros(modelPhy._nb_ctrlpts_total)
usolved[dof] = un
usolved[dod] = ud

# Compare solutions 
_, _, _, u_interp_sample2 = modelPhy.interpolate_field(u_ctrlpts=usolved)
error_2 = np.linalg.norm(u_interp_sample-u_interp_sample2, np.inf)/np.linalg.norm(u_interp_sample, np.inf)*100

print("Error using interpolation : %.3e %%" %(error_1,))
print("Error interpolation/direct solution : %.3e %%" %(error_2,))

# Solution using conjugate gradient
iterations = 80
epsilon = 1e-15

# With preconditioner
method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]
for method_name in method_list: 
    inputs = [Fn, iterations, epsilon, method_name, un]   
    un_t, _, _ = modelPhy.MFsteadyHeat(*inputs)
    usolved_t = np.zeros(modelPhy._nb_ctrlpts_total)
    usolved_t[dof] = un_t
    usolved_t[dod] = ud
    error_3 = np.linalg.norm(usolved-usolved_t, np.inf)/np.linalg.norm(usolved, np.inf)*100
    print("Error direct/iterative solution : %.3e %%" %(error_3,))
