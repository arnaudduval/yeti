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
from pysrc.lib.lib_material import heatmat, mechamat
from pysrc.lib.lib_job import heatproblem, mechaproblem

# Set global variables
geoName = 'QA'
degree_list = np.array([1])
cuts_list   = np.arange(2, 8)
matArgs 	= {'elastic_modulus':1e0, 'elastic_limit':1e10, 'poisson_ratio':0.3,
				'isoHardLaw': {'name':'none'}}

for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
	for i, degree in enumerate(degree_list):
		for j, cuts in enumerate(cuts_list):
			geoArgs = {'name': geoName, 'degree': degree*np.ones(3, dtype=int), 
						'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
						'extra':{'Rin':1.0, 'Rex':2.0}
			}
			quadArgs    = {'quadrule': quadrule, 'type': quadtype}

			blockPrint()			
			modelGeo = Geomdl(geoArgs)
			modelIGA = modelGeo.getIGAParametrization()
			modelPhy = part(modelIGA, quadArgs=quadArgs)

			heatmaterial = heatmat()
			heatmaterial.addCapacity(inpt=1.0, isIsotropic=True)
			heatmaterial.addConductivity(inpt=1.0, isIsotropic=True, shape=2)
			elasticmaterial = mechamat(matArgs=matArgs)

			# Set Dirichlet boundaries
			boundary = boundaryCondition(modelPhy.nbctrlpts)
			boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))
			boundary.add_DirichletDisplacement(table=np.ones((2, 2, 2), dtype=int))
			enablePrint()

			# Solve elastic problem
			heatprob = heatproblem(heatmaterial, modelPhy, boundary); heatprob._itersLin = 500
			mecaprob = mechaproblem(elasticmaterial, modelPhy, boundary); mecaprob._itersLin = 500

			# eigenval = heatprob.compute_eigs()[0]
			eigenval = mecaprob.compute_eigs()[0]
			print('eig:%.3e, degree:%d, nbel:%d' %(eigenval, degree, 2**cuts))
	print('***')