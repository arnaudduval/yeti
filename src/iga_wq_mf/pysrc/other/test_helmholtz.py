"""
	Helmholtz equation:
	-div(u) = lambda u, over Omega
	u = 0, overall boundary Omega
	The lowest eigenvalue of this problem is related to the Poincaré constant C
	from the inequality:
	||u||_{L2} <= C ||grad(u)||_{L2}, since the better approximation C = 1/lambda
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import problem

# Set global variables
geoName = 'QA'
degree_list = np.array([2, 3, 4, 6, 8])
cuts_list   = np.arange(4, 7)
solverArgs  = {'nbIterationsPCG':150, 'PCGThreshold':1e-8}
quadArgs    = {'quadrule': 'iga', 'type': 'leg'}

for i, degree in enumerate(degree_list):
	for j, cuts in enumerate(cuts_list):
		geoArgs = {'name': geoName, 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
					'extra':{'Rin':1.0, 'Rex':2.0}
		}

		blockPrint()			
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)

		# Set Dirichlet boundaries
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))
		enablePrint()

		# Solve elastic problem
		prob = problem(modelPhy, boundary, solverArgs)
		eigenval = prob.compute_eigs()[0]
		print(eigenval)
