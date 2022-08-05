"""
.. Test of plasticity
.. We test how plasticity module works
.. Unities : GPa, mm2, kT
.. Joaquin Cornejo 
"""

# Python libraries
import os
import numpy as np
from scipy import sparse as sp
from matplotlib import pyplot as plt

# My libraries
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

full_path = os.path.realpath(__file__)
folder_file = os.path.dirname(full_path) + '/data/'

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
material = {'density': 7.8e-12, 'young': 210, 'poisson': 0.3, 'sigmaY': 0.1}
modelPhy._set_material(material)

# Set Dirichlet and Neumann boundaries
table = np.zeros((3, 2, 3), dtype=int)
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

# # ==================================
# # ELASTICITY
# # ==================================
# # Initialize
# Fext = Fsurf
# for i in range(3):
#     Fext[i, Mdod[i]] = 0.0

# dof_extended = []
# for i in range(3):
#     dof = np.array(Mdof[i]) + i*modelPhy._nb_ctrlpts_total
#     dof_extended.extend(list(dof))

# # Compute iterative solution in python 
# itersol = modelPhy.MFWQ_solveElasticity(indi=Mdod, b=Fext, isPrecond=True)

# # Compute direct solution
# S2solve = modelPhy.eval_stiffness_matrix()[dof_extended, :][:, dof_extended]
# f2solve = np.reshape(Fext, (-1, 1))[dof_extended]
# dirsol = sp.linalg.spsolve(S2solve, f2solve)
# dirsol = np.atleast_2d(dirsol)

# # Compare results
# itersol2 = np.reshape(itersol, (-1, 1))[dof_extended]
# error = dirsol.T - itersol2
# relerror = np.linalg.norm(error)/np.linalg.norm(dirsol)
# print(relerror)

# # ==================================
# # PLASTICITY
# # ==================================
# # Do ramp function (Fsurf increase linearly)
# nbStep = 5
# dt = 1/nbStep
# Fext = np.zeros((*np.shape(Fvol), nbStep+1))
# for i in range(1, nbStep+1): 
#     Fext[:, :, i] = i*dt*Fsurf

# # Solve system in Python
# modelPhy.MFWQ_solvePlasticity(Fext=Fext[:,:,:2], indi=Mdod)

# # Solve system in fortran
# modelPhy.MFplasticity(u=Fext[:,:,:2], indi=Mdod)

# ==================================
# RANDOMNESS
# ==================================
# Initialize
Fext = Fsurf
for i in range(3):
    Fext[i, Mdod[i]] = 0.0

# Compute iterative solution in python 
coefs = np.load(folder_file+'CoefStiff.npy')

fig, axs = plt.subplots(nrows=9, ncols=9)
for j in range(9):
    for i in range(9):
        x = coefs[i, j, :]
        ax = axs[i, j]
        ax.hist(x)
        ax.set_xticks([], [])
        ax.set_yticks([], [])
    
plt.savefig(folder_file+'output.pdf')

itersol = modelPhy.MFWQ_solveElasticity(coefs=coefs, indi=Mdod, b=Fext, isPrecond=True, isnoised=True)