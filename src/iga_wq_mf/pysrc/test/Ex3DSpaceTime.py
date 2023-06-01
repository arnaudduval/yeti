from lib.__init__ import *
from lib.lib_base import createKnotVector
from lib.lib_geomdl import Geomdl
from lib.lib_part import part
from lib.lib_material import thermomat
from lib.lib_boundary import boundaryCondition
from lib.lib_job_sptm import heatproblemSpTm

def setKprop(P:list):
	cst = 10
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	T = P[3, :]
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((3, 3, len(x)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j] 
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j]*cst*(1.0 + 2.0/(1.0 + np.exp(-5.0*(T-1.0))))
	return Kprop 

def setCprop(P:list):
	cst = 1.0
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	T = P[3, :]
	Cprop = cst*(1 + np.exp(-2.0*abs(T)))
	return Cprop

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t3dtransient/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 6, 4 

# Create model 
geoArgs  = {'name': 'VB', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
quadArgs = {'quadrule': 'wq', 'type': 1}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
model    = part(modelIGA, quadArgs=quadArgs)

# Add material 
material = thermomat()
material.addConductivity(setKprop, isIsotropic=False) 
material.addCapacity(setCprop, isIsotropic=False) 

# Block boundaries
boundary = boundaryCondition(model.nbctrlpts)
boundary.add_DirichletTemperature(table=np.array([[1, 0], [0, 0], [0, 0]]))
boundary.add_DirichletTemperature(table=np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)

# ---------------------
# Transient model
# ---------------------
nbSteps  = 20
timeArgs = {'timespan': 0.25, 'degree': degree, 'knotvector': createKnotVector(degree, nbSteps), 'quadrule': 'wq'}
problem  = heatproblemSpTm(material, model, boundary, timeArgs=timeArgs)

