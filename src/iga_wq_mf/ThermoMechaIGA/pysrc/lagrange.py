"""
.. Test of Lagrange solver
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np, scipy
from copy import deepcopy

# My libraries
# from lib.physics import powden_rotring, temperature_rotring
from lib.physics import powden_cube
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Set global variables
degree = 4
cuts = 4

# Create geometry using geomdl
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('CB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=np.array([cuts, cuts, cuts]))

# ===========================================
# IGA WQ MF APPROACH
# ===========================================
# Interpolation of u
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'capacity':1, 'conductivity':np.eye(3)}
modelPhy._set_material(material)

# Block boundaries
Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
modelPhy._set_dirichlet_boundaries(Dirichlet)
dof = modelPhy._thermal_dof
dod = deepcopy(modelPhy._thermal_dod)

# Compute source vector
F = modelPhy.eval_source_vector(powden_cube)
F[dod] = 0

# Solution using conjugate gradient with preconditioner
iterations = 100; epsilon = 1e-11
inputs = [F, iterations, epsilon]   
usol, relres = modelPhy.MFsteadyHeat_Lagrange(*inputs, indi=dod)
print(relres)

from matplotlib import pyplot as plt
plt.semilogy(relres*100)
plt.xlabel('Number of iterations')
plt.ylabel("Relative error (%)")
plt.savefig('Residue.png')

# modelPhy.export_results(u_ctrlpts=usol)



# ==============================================================================




# # Set global variables
# degree = 4
# cuts = 3 

# # Create geometry using geomdl
# geometry = {'degree':[degree, degree, degree]}
# modelGeo = geomdlModel('RQA', **geometry)
# modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
#                                             np.array([cuts, cuts, cuts]))

# # ===========================================
# # IGA WQ MF APPROACH
# # ===========================================
# # Interpolation of u
# modelPhy = fortran_mf_wq(modelIGA)
# u_interp = modelPhy.interpolate_ControlPoints(temperature_rotring)

# # Add material 
# material = {'capacity':1, 'conductivity':np.eye(3)}
# modelPhy._set_material(material)

# # Block boundaries
# Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
# modelPhy._set_dirichlet_boundaries(Dirichlet)
# dof = modelPhy._thermal_dof
# dod = modelPhy._thermal_dod 

# # Solve system
# ud = u_interp[dod]
# F = modelPhy.eval_source_vector(powden_rotring)

# # Solution using conjugate gradient
# iterations = 20
# epsilon = 1e-15

# # With preconditioner
# inputs = [F, iterations, epsilon]   
# usol = modelPhy.MFsteadyHeat_Lagrange(*inputs, g=ud, indi=dod)[0]
# print(u_interp)
# print(usol)