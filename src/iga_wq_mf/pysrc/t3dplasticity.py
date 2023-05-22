"""
.. Test of elasticity 3D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : Pa (210e9)
..      - Length : m
..      - Force  : N
..      - Mass   : kg 
..      - Density: kg/m^3 (7.8e3)
..      - Gravity: m/s^2 (9.8)
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.lib_geomdl import Geomdl
from lib.lib_model import part
from lib.lib_material import mechamat
from lib.lib_step import step
from lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t3delastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 2, 3
name = 'CB'

# Create model 
inputs = {'name': name, 'degree':degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
modelGeo = Geomdl(**inputs)
modelIGA = modelGeo.getIGAParametrization()
model    = part(modelIGA)

# Add material 
kwargs = {'density': 7800, 'elastic_modulus': 1e9, 'poisson_ratio': 0.3, 'elastic_limit': 500e9, 
			'law':{'name': 'linear', 'Hbar':1445, 'theta':1.0}}
material = mechamat(kwargs)

# Set Dirichlet boundaries
boundary = step(model._nbctrlpts)
table = np.zeros((3, 2, 3), dtype=int)
table[0, 0, 0] = 1
table[1, 0, 1] = 1
table[2, 0, 2] = 1
boundary.add_DirichletDisplacement(table=table)

# Elasticity problem
problem = mechaproblem(material, model, boundary)

# Set Neumann boundaries
def forceFun(P:list):
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	ref  = np.array([1e8, 2e8, 0.0])
	prop = np.zeros((3, len(x)))
	for i in range(3): prop[i, :] = ref[i] 
	return prop

Fsurf = problem.eval_surfForce(forceFun, nbFacePosition=1)
nbStep = 6; dt = 1/nbStep
Fext   = np.zeros((*np.shape(Fsurf), nbStep+1))
for i in range(1, nbStep+1): Fext[:, :, i] = i*dt*Fsurf

displacement, stress_vm = problem.solvePlasticityProblemPy(Fext=Fext)