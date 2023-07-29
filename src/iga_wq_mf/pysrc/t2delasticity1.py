"""
.. Test of elastoplasticity 2D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa (200e3)
..      - Length : mm
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.lib_geomdl import Geomdl
from lib.lib_part import part
from lib.lib_material import (mechamat, block_dot_product)
from lib.lib_load import forceSurf, dispInfinitePlate
from lib.lib_boundary import boundaryCondition
from lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t2delasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
E, nu = 1e3, 0.3
trueEnergy = -135/32768*np.pi/E*(1024*nu**2 + 5*nu - 1019)
name = 'QA'

quadArgs = {'quadrule': 'iga', 'type': 'leg'}
# quadArgs = {'quadrule': 'wq', 'type': 2}
matArgs  = {'elastic_modulus':E, 'elastic_limit':1e10, 'poisson_ratio': nu}
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-16}

# Create model 
degree_list = np.array([2, 3, 4, 6, 8])
cuts_list   = np.arange(2, 9)
error_energy = np.ones(len(cuts_list))
L2error_disp = np.ones(len(cuts_list))
fig, ax  = plt.subplots()
for i, degree in enumerate(degree_list):
	for j, cuts in enumerate(cuts_list):
		geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
					'extra':{'Rin':1.0, 'Rex':4.0}
		}
		blockPrint()
		material = mechamat(matArgs)
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		model    = part(modelIGA, quadArgs=quadArgs)

		# Set Dirichlet boundaries
		boundary = boundaryCondition(model.nbctrlpts)
		table = np.zeros((2, 2, 2), dtype=int)
		table[1, 1, 0] = 1
		table[1, 0, 1] = 1
		boundary.add_DirichletDisplacement(table=table)
		enablePrint()

		# Solve elastic problem
		problem = mechaproblem(material, model, boundary)
		problem.addSolverConstraints(solverArgs=solverArgs)
		Fext = problem.eval_surfForce(forceSurf, nbFacePosition=1)
		displacement, _, stress_qp = problem.solveElasticityProblemFT(Fext=Fext)
		error_energy[j] = abs(trueEnergy -  block_dot_product(2, Fext, displacement))/trueEnergy*100
		# L2error_disp[j] = problem.L2NormOfError(dispInfinitePlate, displacement)

	elsize = 1.0/2**cuts_list
	ax.loglog(elsize, error_energy, marker=markerSet[i], label='degree p='+str(degree))
	# ax.loglog(elsize, L2error_disp, marker=markerSet[i], label='degree p='+str(degree))

	# if str(quadArgs['quadrule']) == 'wq':
	# 	pass
	# 	# slope = np.polyfit(np.log10(elsize[2:5]),np.log10(error_energy[2:5]), 1)[0]
	# 	# slope = round(slope, 1)
	# 	# annotation.slope_marker((elsize[3], error_energy[3]), slope, 
	# 	# 						poly_kwargs={'facecolor': (0.73, 0.8, 1)})
	# else: 
	# 	slope = np.polyfit(np.log10(elsize[:3]),np.log10(error_energy[:3]), 1)[0]
	# 	slope = round(slope, 1)
	# 	annotation.slope_marker((elsize[2], error_energy[2]), slope, 
	# 							poly_kwargs={'facecolor': (0.73, 0.8, 1)})
	
	ax.set_ylabel('Relative error of energy')
	ax.set_ylabel('L2 norm error')
	ax.set_xlabel('Meshsize h')
	ax.set_ylim(top=1e1, bottom=1e-14)
	ax.set_xlim(left=2e-3, right=0.8)

	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.tight_layout()
	fig.savefig(folder + 'FigInfinitePlate2_' + str(quadArgs['quadrule']) +'.png')
