"""
.. Test of plasticity 3D
.. We test how plasticity module works
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
from lib.lib_part import part
from lib.lib_material import mechamat
from lib.lib_boundary import boundaryCondition
from lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t3delastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 2, 3
name = 'CB'

# Create model 
geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
quadArgs  = {'quadrule': 'wq', 'type': 1}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
model    = part(modelIGA, quadArgs=quadArgs)

# Add material 
matArgs  = {'density': 7800, 'elastic_modulus':200e3, 'elastic_limit':506, 'poisson_ratio': 0.3,
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}
material = mechamat(matArgs)

# Set Dirichlet boundaries
boundary = boundaryCondition(model.nbctrlpts)
table = np.zeros((3, 2, 3), dtype=int)
table[0, 0, 0] = 1
table[1, 0, 1] = 1
table[2, 0, 2] = 1
boundary.add_DirichletDisplacement(table=table)

# Elasticity problem
problem = mechaproblem(material, model, boundary)

def forceSurfFun(P:list):
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	ref  = np.array([1e8, 2e8, 0.0])
	prop = np.zeros((3, len(x)))
	for i in range(3): prop[i, :] = ref[i] 
	return prop

Fsurf = problem.eval_surfForce(forceSurfFun, nbFacePosition=1)
nbStep = 6; dt = 1/nbStep
Fext   = np.zeros((*np.shape(Fsurf), nbStep+1))
for i in range(1, nbStep+1): Fext[:, :, i] = i*dt*Fsurf

displacement = problem.solvePlasticityProblemPy(Fext=Fext)[0]