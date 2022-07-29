"""
.. Test of plasticity
.. We test how plasticity module works
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
cuts = 2

# Create geometry using geomdl
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('VB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create physical model
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'density': 7.8e-6, 'young': 210, 'poisson': 0.3, 'sigmaY': 80}
modelPhy._set_material(material)

# Set Dirichlet and Neumann boundaries
table = np.zeros((3, 2, 3))
table[0, :, :] = 1
Dirichlet = {'mechanical':table}
modelPhy._set_dirichlet_boundaries(Dirichlet)

forces = [[0 for i in range(3)] for j in range(6)]
forces[1] = [10, 20, 100]
Neumann = {'mechanical': forces}
modelPhy._set_neumann_boundaries(Neumann)

# Set external forces
Fvol = modelPhy.eval_force_body(bodyforce)
Fsurf = modelPhy.eval_force_surf()

# # Do ramp function (Fvol is constant, but Fsurf increase linearly)
# nbStep = 100
# dt = 1/nbStep
# Fext = np.zeros((len(Fvol), nbStep+1))
# for i in range(nbStep+1): 
#     Fext[:, i] = Fvol + dt*Fsurf

# # Solve system
# disp = modelPhy.MFplasticity(Fext)

Fext = Fvol + Fsurf
result = modelPhy.eval_Su(Fext)