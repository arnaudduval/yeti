"""
.. Test of plasticity
.. We test how plasticity module works
.. Unities : GPa, mm2, kT
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np
from scipy import sparse as sp

# My libraries
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

def bodyforce(P: list):
    " Gravity "
    # Isotropy
    f = np.zeros(np.shape(P))
    f[-1, :] = 9800
    return f

# Set global variables
degree = 4
cuts = 3

# Create geometry using geomdl
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('CB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create physical model
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'density': 7.8e-12, 'young': 210, 'poisson': 0.3, 'sigmaY': 0.2}
modelPhy._set_material(material)

# Set Dirichlet and Neumann boundaries
table = np.zeros((3, 2, 3))
table[0, 0, :] = 1
Dirichlet = {'mechanical':table}
_, _, _, _, Mdof, Mdod = modelPhy._set_dirichlet_boundaries(Dirichlet)

forces = [[0 for i in range(3)] for j in range(6)]
forces[1] = [0.1, 0.0, 0.0]
Neumann = {'mechanical': forces}
modelPhy._set_neumann_boundaries(Neumann)

# Set external forces
Fvol = modelPhy.eval_force_body(bodyforce)
Fsurf = modelPhy.eval_force_surf()

# # ==========================================
# # Compute Stiffness matrix
# ones = Fsurf
# S = modelPhy.eval_stiffness_matrix()
# ones_row = np.reshape(ones, (-1, 1))
# R1 = S @ ones_row
# R2 = modelPhy.eval_Su(ones)
# R2 = np.reshape(R2, (-1, 1))
# diff = R1 - R2
# error = np.linalg.norm(diff)
# print(error)
# # ==========================================

# # Solver elastic
# Fext = Fsurf
# for i in range(3):
#     Fext[i, Mdod[i]] = 0.0
# itersol = modelPhy.MFelasticity(dod=Mdod, u=Fext)

# dof_extended = []
# for i in range(3):
#     dof = np.array(Mdof[i]) + i*modelPhy._nb_ctrlpts_total
#     dof_extended.extend(list(dof))

# # Compute Stiffness matrix
# S = modelPhy.eval_stiffness_matrix()[dof_extended, :][:, dof_extended]
# F2solve = np.reshape(Fext, (-1, 1))[dof_extended]
# dirsol = sp.linalg.spsolve(S, F2solve)
# dirsol = np.atleast_2d(dirsol)

# itersol2 = np.reshape(itersol, (-1, 1))[dof_extended]
# error = dirsol.T - itersol2
# norm_error = np.linalg.norm(error)/np.linalg.norm(dirsol)
# print(norm_error)

# # Export results
# modelPhy.export_results(result)

# ============================================
# Do ramp function (Fvol is constant, but Fsurf increase linearly)
nbStep = 5
dt = 1/nbStep
Fext = np.zeros((*np.shape(Fvol), nbStep+1))
for i in range(1, nbStep+1): 
    Fext[:, :, i] = i*dt*Fsurf

for i in range(3):
    Fext[i, Mdod[i], :] = 0.0 # Is zero but to be sure

# Solve system
_, delta_F, coef_S= modelPhy.MFplasticity(Mdod, Fext[:,:,:2])

# Solve with iterative method
itersol = modelPhy.MFelasticity(coefs = coef_S, dod= Mdod, u = delta_F)
print(np.max(itersol), np.min(itersol))

# -----------------
print('=========================')
dof_extended = []
for i in range(3):
    dof = np.array(Mdof[i]) + i*modelPhy._nb_ctrlpts_total
    dof_extended.extend(list(dof))

# Compute Stiffness matrix
S = modelPhy.eval_stiffness_matrix(coefs = coef_S)[dof_extended, :][:, dof_extended]
F2solve = np.reshape(delta_F, (-1, 1))[dof_extended]
dirsol = sp.linalg.spsolve(S, F2solve)
dirsol = np.atleast_2d(dirsol)
print(np.max(dirsol), np.min(dirsol))

itersol2 = np.reshape(itersol, (-1, 1))[dof_extended]
error = dirsol.T - itersol2
norm_error = np.linalg.norm(error)/np.linalg.norm(dirsol)
print(norm_error)
