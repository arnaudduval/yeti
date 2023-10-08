from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat, mechamat
from pysrc.lib.lib_job import thermomechaproblem

# Set global variables
degree, cuts = 3, 3
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-8}
quadArgs   = {'quadrule': 'iga', 'type': 'leg'}

def conductivityProperty(P:list):
	Kref  = np.array([[1, 0.5],[0.5, 2]])
	Kprop = np.zeros((2, 2, np.size(P, axis=1)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j] 
	return Kprop 

def capacityProperty(P:list):
	cst = 1.0
	T   = P[3, :]
	Cprop = cst*(1 + np.exp(-2.0*abs(T)))
	return Cprop

geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'Rin':1.0, 'Rex':2.0}
}

blockPrint()	
heatmaterial = heatmat()
heatmaterial.addConductivity(conductivityProperty, isIsotropic=False)	
heatmaterial.addCapacity(capacityProperty, isIsotropic=False) 
mechamaterial = mechamat({'elastic_modulus':1e3, 'elastic_limit':1e10, 'poisson_ratio':0.3,
		'plasticLaw': {'Isoname':'none'}})

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)

# Set Dirichlet boundaries
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))
enablePrint()

# Solve elastic problem
problem = thermomechaproblem(heatmaterial, mechamaterial, modelPhy, boundary, solverArgs)
