"""
.. Test of elastoplasticity 2D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa (200e3)
..      - Length : mm
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat, block_dot_product
from pysrc.lib.lib_load import forceSurf_infPlate
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2delasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
E, nu = 1e3, 0.3
trueEnergy = -135/32768*np.pi/E*(1024*nu**2 + 5*nu - 1019)

# quadArgs = {'quadrule': 'iga', 'type': 'leg'}
quadArgs = {'quadrule': 'wq', 'type': 2}
matArgs  = {'elastic_modulus':E, 'elastic_limit':1e10, 'poisson_ratio': nu}
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-16}

# Create model 
degree_list = np.array([2, 3, 4, 6, 8])
cuts_list   = np.arange(2, 9)
error = np.ones(len(cuts_list))

fig, ax  = plt.subplots(figsize=(8, 4))
for i, degree in enumerate(degree_list):
	for j, cuts in enumerate(cuts_list):
		geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
					'extra':{'Rin':1.0, 'Rex':4.0}
		}
		blockPrint()
		material = mechamat(matArgs)
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)

		# Set Dirichlet boundaries
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		table = np.zeros((2, 2, 2), dtype=int)
		table[1, 1, 0] = 1
		table[1, 0, 1] = 1
		boundary.add_DirichletDisplacement(table=table)
		enablePrint()

		# Solve elastic problem
		problem = mechaproblem(material, modelPhy, boundary)
		problem.addSolverConstraints(solverArgs=solverArgs)
		Fext = problem.eval_surfForce(forceSurf_infPlate, nbFacePosition=1)
		displacement, _, stress_qp = problem.solveElasticityProblemFT(Fext=Fext)
		error[j] = abs(trueEnergy -  block_dot_product(2, Fext, displacement))/trueEnergy*100

	nbctrlpts = (2**cuts_list+degree)**2
	ax.loglog(nbctrlpts, error, marker=markerSet[i], label='degree p='+str(degree))

	if str(quadArgs['quadrule']) == 'iga':
		slope = np.polyfit(np.log10(nbctrlpts[:3]),np.log10(error[:3]), 1)[0]
		slope = round(slope, 1)
		annotation.slope_marker((nbctrlpts[2], error[2]), slope, 
								poly_kwargs={'facecolor': (0.73, 0.8, 1)})
	
	ax.set_ylabel('Relative error of energy ' + r'$\frac{|U-U^{h}|}{|U|}$')
	ax.set_xlabel('Total number of DOF')
	ax.set_ylim(top=1e1, bottom=1e-14)
	ax.set_xlim(left=10, right=1e5)

	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.tight_layout()
	fig.savefig(folder + 'FigInfinitePlate_' + str(quadArgs['quadrule']) +'.png')
