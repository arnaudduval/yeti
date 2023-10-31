from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat, mechamat
from pysrc.lib.lib_job import thermomechaproblem
from pysrc.lib.lib_base import sigmoid

# Set global variables
DEGREE, CUTS = 3, 3
TRACTION, RINT, REXT = 1.0, 1.0, 2.0
YOUNG, POISSON = 1e5, 0.3
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
GEONAME = 'QA'
GEOARGS = {'name': 'QA', 'degree': DEGREE*np.ones(3, dtype=int), 
			'nb_refinementByDirection': CUTS*np.ones(3, dtype=int), 
			'extra':{'Rin':RINT, 'Rex':REXT}
}
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1.5, 'poisson_ratio': POISSON, 
			'plasticLaw': {'Isoname':'none'}}
SOLVERARGS = {'nbIterationsPCG':150, 'PCGThreshold':1e-9, 'PCGmethod': 'TDC', 'NRThreshold': 1e-6}
QUADARGS   = {'quadrule': 'iga', 'type': 'leg'}

def conductivityProperty(P:list):
	Kref  = np.array([[1, 0.5],[0.5, 2]])
	Kprop = np.zeros((2, 2, np.size(P, axis=1)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j] 
	return Kprop 

def capacityProperty(P:list):
	cst = 1.0
	T   = P[-1, :]
	Cprop = cst*(1 + np.exp(-2.0*abs(T)))
	return Cprop

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = RINT**2/r_square
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = TRACTION/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = TRACTION/2*3*np.sin(3*theta)*(b**2 - b)
	return F

blockPrint()	
heatmaterial = heatmat()
heatmaterial.addConductivity(conductivityProperty, isIsotropic=False)	
heatmaterial.addCapacity(capacityProperty, isIsotropic=False) 
mechamaterial = mechamat(MATARGS)

modelGeo = Geomdl(GEOARGS)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=QUADARGS)

# Set Dirichlet boundaries
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.array([[1, 0], [0, 0], [0, 0]]))
boundary.add_DirichletConstTemperature(table=np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)
table = np.zeros((2, 2, 2), dtype=int); table[1, 1, 0] = 1; table[1, 0, 1] = 1
boundary.add_DirichletDisplacement(table=table)
enablePrint()

# Solve elastic problem
problem = thermomechaproblem(heatmaterial, mechamaterial, modelPhy, boundary)
problem.addSolverConstraints(solverArgs=SOLVERARGS)
problem.addDensity(1.0, isIsotropic=True)

# Set external forces
Fref = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
Fmech_list = np.zeros((2, modelPhy.nbctrlpts_total, len(TIME_LIST)))
for k in range(len(TIME_LIST)): Fmech_list[:, :, k] = np.sin(TIME_LIST[k])*Fref

Fref = np.zeros((problem.part.nbctrlpts_total, 1))
Fheat_list = np.kron(Fref, sigmoid(TIME_LIST))

# Set initial conditions
temperature = np.zeros(np.shape(Fheat_list))
for i in range(1, len(TIME_LIST)): temperature[boundary.thdod, i] = boundary.thDirichletBound[boundary.thdod]

displacement = np.zeros(np.shape(Fmech_list))
problem.solveThermoElasticityProblem(displacement, temperature, Fmech_list, Fheat_list, TIME_LIST, isLumped=True)
