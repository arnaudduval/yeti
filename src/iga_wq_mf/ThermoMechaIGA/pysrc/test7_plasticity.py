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
material = {'density': 7.8e-12, 'young': 210, 'poisson': 0.3, 'sigmaY': 0.11}
modelPhy._set_material(material)

# Set Dirichlet and Neumann boundaries
table = np.zeros((3, 2, 3))
table[0, 0, :] = 1
Dirichlet = {'mechanical':table}
_, _, _, _, Mdof, Mdod = modelPhy._set_dirichlet_boundaries(Dirichlet)

forces = [[0 for i in range(3)] for j in range(6)]
forces[1] = [0.8, 0.0, 0.0]
Neumann = {'mechanical': forces}
modelPhy._set_neumann_boundaries(Neumann)

# Set external forces
Fvol = modelPhy.eval_force_body(bodyforce)
Fsurf = modelPhy.eval_force_surf()

# # ==========================================
# # Compute S.u
# u = Fsurf
# S = modelPhy.eval_stiffness_matrix()
# urow = np.reshape(u, (-1, 1))
# R1 = S @ urow
# R2 = modelPhy.eval_Su(u)
# R2 = np.reshape(R2, (-1, 1))
# diff = R1 - R2
# error = np.linalg.norm(diff)
# print(error)
# # ==========================================

# # Solver elastic
# u = Fsurf
# for i in range(3):
#     u[i, Mdod[i]] = 0.0
# itersol = modelPhy.MFelasticity(u=u, dod=Mdod)

# dof_extended = []
# for i in range(3):
#     dof = np.array(Mdof[i]) + i*modelPhy._nb_ctrlpts_total
#     dof_extended.extend(list(dof))

# # Compute Stiffness matrix
# S2solve = modelPhy.eval_stiffness_matrix()[dof_extended, :][:, dof_extended]
# u2solve = np.reshape(u, (-1, 1))[dof_extended]
# dirsol = sp.linalg.spsolve(S2solve, u2solve)
# dirsol = np.atleast_2d(dirsol)

# itersol2 = np.reshape(itersol, (-1, 1))[dof_extended]
# error = dirsol.T - itersol2
# norm_error = np.linalg.norm(error)/np.linalg.norm(dirsol)
# print(norm_error)

# ============================================
# Do ramp function (Fsurf increase linearly)
nbStep = 5
dt = 1/nbStep
u = np.zeros((*np.shape(Fvol), nbStep+1))
for i in range(1, nbStep+1): 
    u[:, :, i] = i*dt*Fsurf

for i in range(3):
    u[i, Mdod[i], :] = 0.0 # Is zero but to be sure

# Solve system in fortran
_, delta_F, coef_S = modelPhy.MFplasticity(u=u[:,:,:2], indi=Mdod)
modelPhy._set_extra_mechanical_properties()
coef_S_el = modelPhy.compute_elastic_coefficient(modelPhy._Jqp)
print(coef_S[:, :, 100])
print(coef_S_el[:, :, 100])


print('=========================')

# # Solve with iterative method in fortran
# itersol = modelPhy.MFelasticity(u=delta_F, indi=Mdod, coefs=coef_S)
# print(np.max(itersol), np.min(itersol))

# # Solve with direct method
# dof_extended = []
# for i in range(3):
#     dof = np.array(Mdof[i]) + i*modelPhy._nb_ctrlpts_total
#     dof_extended.extend(list(dof))

# # Compute Stiffness matrix
# S2solve = modelPhy.eval_stiffness_matrix(coefs = coef_S)[dof_extended, :][:, dof_extended]
# delta_F_dir = np.reshape(delta_F, (-1, 1))[dof_extended]
# dirsol = sp.linalg.spsolve(S2solve, delta_F_dir)
# dirsol = np.atleast_2d(dirsol)
# print(np.max(dirsol), np.min(dirsol))

# itersol2 = np.reshape(itersol, (-1, 1))[dof_extended]
# error = dirsol.T - itersol2
# norm_error = np.linalg.norm(error)/np.linalg.norm(dirsol)
# print(norm_error)
