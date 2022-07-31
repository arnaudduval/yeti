"""
.. Test of plasticity
.. We test how plasticity module works
.. Unities : GPa, mm2, kT
.. Joaquin Cornejo 
"""

# Python libraries
import numpy as np

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
material = {'density': 7.8e-12, 'young': 210, 'poisson': 0.3, 'sigmaY': 0.01}
modelPhy._set_material(material)

# Set Dirichlet and Neumann boundaries
table = np.zeros((3, 2, 3))
table[0, 0, :] = 1
Dirichlet = {'mechanical':table}
Mdod = modelPhy._set_dirichlet_boundaries(Dirichlet)[-1]

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
# Fext = Fvol + Fsurf
# for i in range(3):
#     Fext[i, Mdod[i]] = 0.0
# result = modelPhy.MFelasticity(Mdod, Fext)

# # Export results
# modelPhy.export_results(result)

# ============================================
# Do ramp function (Fvol is constant, but Fsurf increase linearly)
nbStep = 5
dt = 1/nbStep
Fext = np.zeros((*np.shape(Fvol), nbStep+1))
for i in range(1, nbStep+1): 
    Fext[:, :, i] = Fvol + i*dt*Fsurf

for i in range(3):
    Fext[i, Mdod[i], :] = 0.0

# Solve system
disp = modelPhy.MFplasticity(Mdod, Fext)
