"""
.. Test of substitution solver
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
case = 1
degree, cuts = 4, 5
method_precond = 'JMC'
nbIterPCG, threshold = 30, 1e-10

if case == 1:

    # Create geometry and model
    geometry = {'degree':[degree, degree, degree]}
    modelGeo = geomdlModel('VB', **geometry)
    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                np.array([cuts, cuts, cuts]))
    modelPhy = fortran_mf_wq(modelIGA)

    # Add material 
    conductivity = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
    material = {'capacity':1, 'conductivity':conductivity}
    modelPhy._set_material(material)

    # Block boundaries
    Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
    modelPhy._set_dirichlet_boundaries(Dirichlet)

    # Compute vectors
    ud = None
    F = modelPhy.eval_source_vector(powden_prism)

elif case == 2:

    # Create geometry and model
    geometry = {'degree':[degree, degree, degree]}
    modelGeo = geomdlModel('RQA', **geometry)
    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                np.array([cuts, cuts, cuts]))
    modelPhy = fortran_mf_wq(modelIGA)
    
    # Interpolation of u
    u_interp = modelPhy.interpolate_ControlPoints(temperature_rotring)

    # Add material 
    material = {'capacity':1, 'conductivity':np.eye(3)}
    modelPhy._set_material(material)

    # Block boundaries
    Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
    modelPhy._set_dirichlet_boundaries(Dirichlet)
    dod = modelPhy._thermal_dod 

    # Compute vectors
    ud = u_interp[dod]
    F = modelPhy.eval_source_vector(powden_rotring)

# Solve system
inputs = [F, nbIterPCG, threshold]   
usol, residue = modelPhy.MFsteadyHeat_PLS(*inputs, ud=ud, methodPCG=method_precond)

# Plot
newres = residue[np.nonzero(residue)]
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5,4))
ax.semilogy(np.arange(1, len(newres)+1), abs(newres))

# Set properties
ax.set_xlabel('Number of iterations')
ax.set_ylabel("Relative residue")
fig.tight_layout()
fig.savefig(folder + 'Substitution.png')