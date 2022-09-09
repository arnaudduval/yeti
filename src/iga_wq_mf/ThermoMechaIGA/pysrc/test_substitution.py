"""
.. Test of substitution solver
.. Joaquin Cornejo 
"""

# Python libraries
from matplotlib import pyplot as plt
import numpy as np

# My libraries
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Set global variables
case = 1
iterations = 30
method_precond = 'WP'
degree, cuts = 4, 5
epsilon = 1e-10

if case == 1:

    # Create geometry using geomdl
    geometry = {'degree':[degree, degree, degree]}
    modelGeo = geomdlModel('VB', **geometry)
    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                np.array([cuts, cuts, cuts]))

    # Interpolation of u
    modelPhy = fortran_mf_wq(modelIGA)

    # Add material 
    conductivity = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
    material = {'capacity':1, 'conductivity':conductivity}
    modelPhy._set_material(material)

    # Block boundaries
    Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
    modelPhy._set_dirichlet_boundaries(Dirichlet)
    dod = modelPhy._thermal_dod 

    # Solve system
    F = modelPhy.eval_source_vector(powden_prism)

    # With preconditioner
    inputs = [F, iterations, epsilon]   
    usol, residue = modelPhy.MFsteadyHeat_PLS(*inputs, indi=dod, method_precond=method_precond)

elif case == 2:

    # Create geometry using geomdl
    geometry = {'degree':[degree, degree, degree]}
    modelGeo = geomdlModel('RQA', **geometry)
    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                np.array([cuts, cuts, cuts]))

    # Interpolation of u
    modelPhy = fortran_mf_wq(modelIGA)
    u_interp = modelPhy.interpolate_ControlPoints(temperature_rotring)

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
    F = modelPhy.eval_source_vector(powden_rotring)

    # With preconditioner
    inputs = [F, iterations, epsilon]   
    usol, residue = modelPhy.MFsteadyHeat_PLS(*inputs, ud=ud, indi=dod, method_precond=method_precond)

# Print results
CARO = 2
newres = residue[np.nonzero(residue)]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
ax.semilogy(np.arange(1, len(newres)+1), abs(newres))
ax.set_xlabel('Number of iterations')
ax.set_ylabel("Relative residue")

plt.grid()
plt.tight_layout()
plt.savefig('Substitution.png')