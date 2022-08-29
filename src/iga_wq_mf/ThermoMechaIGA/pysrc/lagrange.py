"""
.. Test of Lagrange solver
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np, scipy
from copy import deepcopy

# My libraries
from lib.physics import powden_rotring, temperature_rotring
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
def temperature_cube(P: list):
    " T = sin(pi*x)*sin(pi*y)*sin(pi*z)"
    x = P[0, :]
    y = P[1, :]
    z = P[2, :]
    f = np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)

    return f

# Interpolation of u
modelPhy = fortran_mf_wq(modelIGA)
u_interp = modelPhy.interpolate_ControlPoints(temperature_cube)

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
Fn = F[dof]

# Solution using conjugate gradient with preconditioner (Lagrange or penalty)
iterations = 5; epsilon = 1e-3
inputs = [F, iterations, epsilon]   
usol_0, energy = modelPhy.MFsteadyHeat_PLS(*inputs, indi=dod, method_pls='S', method_precond='WP')
error = abs(u_interp - usol_0)/abs(u_interp).max()

from matplotlib import pyplot as plt
plt.title('Energy last value: %.5f' %energy[-1])
plt.plot(np.arange(1, len(energy)+1), energy)
plt.xlabel('Number of iterations')
plt.ylabel("Energy")
plt.tight_layout()
plt.savefig('Energy_Substitution.png')

# # Define residue R = F - KT (all variables)
# KU0 = modelPhy.eval_Ku(usol_0)
# R0 = abs(F - KU0)
# E0 = 0.5*usol_0 @ KU0 - usol_0 @ F

# # Solution solving Knn un = Fn, in this example ud = 0
# inputs = [Fn, iterations, epsilon, "JM"]   
# un_1, _, _ = modelPhy.MFsteadyHeat(*inputs)
# usol_1 = np.zeros(modelPhy._nb_ctrlpts_total)
# usol_1[dof] = un_1
# error = abs(u_interp - usol_1)/abs(u_interp).max()

# # Define residue R = F - KT (all variables)
# KU1 = modelPhy.eval_Ku(usol_1)
# R1 = abs(F - KU1)
# E1 = 0.5* usol_1 @ KU1 - usol_1 @ F

# print('Energy : %.5f, %.5f'  %(E0, E1))
# print("Residue : %.5f, %.5f" %(R0.max(), R1.max()))

# # modelPhy.export_results(u_ctrlpts=usol_0)

# # ==============================================================================

# from lib.__init__ import blockPrint

# blockPrint()
# # Set global variables
# degree = 4
# cuts = 4

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
# iterations = 200
# epsilon = 1e-15

# # With preconditioner
# inputs = [F, iterations, epsilon]   
# usol, energy = modelPhy.MFsteadyHeat_PLS(*inputs, ud=ud, indi=dod, method_pls=True)
# error = abs(u_interp - usol)/abs(u_interp).max()

# from matplotlib import pyplot as plt
# plt.semilogy(np.arange(1, len(energy)+1), abs(energy))
# plt.title('Energy last value: %.5f' %energy[-1])
# plt.xlabel('Number of iterations')
# plt.ylabel("Energy")
# plt.grid()
# plt.tight_layout()
# plt.savefig('Energy_Lagrange.png')