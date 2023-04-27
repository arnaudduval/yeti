from lib.__init__ import *
from lib.lib_material import thermomat
from lib.lib_geomdl import Geomdl
from lib.lib_model import part
from lib.lib_step import step
from lib.lib_job import heatproblem
from lib.lib_load import powden_prism

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test6/'
if not os.path.isdir(folder): os.mkdir(folder)

# Create geometry and model
degree, cuts = 4, 4
kwargs   = {'name': 'vb', 'degree': degree*np.ones(3, dtype=int),
		'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
modelgeo = Geomdl(**kwargs)
modeliga = modelgeo.getIGAParametrization()
model    = part(modeliga)

# Create material 
material = thermomat()
material.addConductivity(np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]]), 
                        isIsotropic=True, shape=(3, 3))

# Set boundary conditions
table    = np.ones((3, 2), dtype=bool)
boundary = step(model._nbctrlpts)
boundary.add_DirichletTemperature(table=table)

problem = heatproblem(material, model, boundary)

# Add heat force
Fext = problem.eval_heatForce(fun=powden_prism)
sol, residue = problem.solveSteadyHeatProblemFT(Fext)
model.exportResults(u_ctrlpts=sol, folder=folder, nbDOF=1)