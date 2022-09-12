"""
.. Test of transient heat solver
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
degree, cuts = 4, 4

# Create geometry 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('VB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create model
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'capacity':3, 'conductivity':np.eye(3)}
modelPhy._set_material(material)

# Block boundaries
Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
modelPhy._set_dirichlet_boundaries(Dirichlet)

# Create source
Fext = modelPhy.eval_source_vector(powden_prism)

fig, ax = plt.subplots(nrows=1, ncols=1)

for itmethod in ['FDC', 'JMC']:
    # Solve linear transient problem
    sol, residue = modelPhy.MF_THL(F=Fext, method=itmethod)
    newresidue = residue[residue>0]*100
    ax.semilogy(np.arange(len(newresidue )), newresidue, label=itmethod)

# Set properties
ax.legend()
ax.set_xlabel('Number of iterations')
ax.set_ylabel('Relative error (%)')
ax.set_ylabel('Relative residue (%)')
fig.tight_layout()
fig.savefig(folder + 'Transient3D.png')
