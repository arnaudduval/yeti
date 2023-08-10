"""
.. Test of elasticity 2D
.. We test how elasticity module works
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2delasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceSurf_infPlate(P:list):
	Tx, a = 1.0, 1.0
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = a**2/r_square # Already squared
	theta = np.arcsin(y/np.sqrt(r_square))

	# theta_x2 = 2.0*theta
	# PolarSt = np.zeros((3, nnz))
	# PolarSt[0, :] = Tx/2.0*(1.0 - b + np.cos(theta_x2)*(1.0 - 4.0*b + 3.0*b**2)) 
	# PolarSt[1, :] = Tx/2.0*(1.0 + b - np.cos(theta_x2)*(1.0 + 3.0*b**2))
	# PolarSt[2, :] = -Tx/2.0*(1.0 + 2.0*b - 3.0*b**2)*np.sin(theta_x2)	
	
	# CartSt = np.zeros((3, nnz))
	# CartSt[0, :] = PolarSt[0, :]*(np.cos(theta))**2 + PolarSt[1, :]*(np.sin(theta))**2 - PolarSt[2, :]*(np.sin(theta_x2)) 
	# CartSt[1, :] = PolarSt[0, :]*(np.sin(theta))**2 + PolarSt[1, :]*(np.cos(theta))**2 + PolarSt[2, :]*(np.sin(theta_x2)) 
	# CartSt[2, :] = PolarSt[2, :]*np.cos(theta_x2) + 0.5*np.sin(theta_x2)*(PolarSt[0, :] - PolarSt[1, :])

	# F = np.zeros((2, nnz))
	# F[0, :] = CartSt[0, :]*np.cos(theta) + CartSt[2, :]*np.sin(theta)
	# F[1, :] = CartSt[2, :]*np.cos(theta) + CartSt[1, :]*np.sin(theta)

	F = np.zeros((2, nnz))
	F[0, :] = Tx/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = Tx/2*3*np.sin(3*theta)*(b**2 - b)
	return F

def exactDisplacement_infPlate(P:list):
	Tx, a, E, nu = 1.0, 1.0, 1.e3, 0.3
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	theta = np.arcsin(y/np.sqrt(r_square))
	b = a**2/r_square # Already squared
	c = Tx*(1.0 + nu)*np.sqrt(r_square)/(2*E)

	# theta_x2  = 2.0*theta
	# PolarDisp = np.zeros((2, nnz))
	# PolarDisp[0, :] = c*((b + 1.0 - 2.0*nu) + np.cos(theta_x2)*(-b**2 + 4.0*(1.0 - nu)*b + 1.0))
	# PolarDisp[1, :] = -c*np.sin(theta_x2)*(b**2 + 2.0*(1.0 - 2.0*nu)*b + 1.0)

	# disp  = np.zeros((2, nnz))
	# disp[0, :] = PolarDisp[0, :]*np.cos(theta) - PolarDisp[1, :]*np.sin(theta)
	# disp[1, :] = PolarDisp[0, :]*np.sin(theta) + PolarDisp[1, :]*np.cos(theta)

	disp = np.zeros((2, nnz))
	disp[0, :] = c*(2*(1-nu)*np.cos(theta) + b*(4*(1-nu)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	disp[1, :] = c*(-2*nu*np.sin(theta) + b*(2*(-1 + 2*nu)*np.sin(theta) + np.sin(3*theta)) - b**2*np.sin(3*theta))
	
	return disp

# Set global variables
E, nu = 1e3, 0.3
matArgs    = {'elastic_modulus':E, 'elastic_limit':1e10, 'poisson_ratio': nu}
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-15, 'PCGmethod': 'TDC'}
degree_list = np.array([2, 3, 4, 6, 8])
cuts_list   = np.arange(2, 9)

for quadrule, quadtype in zip(['wq', 'iga'], [2, 'leg']):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	error_list = np.ones(len(cuts_list))
	fig, ax    = plt.subplots(figsize=(8, 4))

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
			Fext = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
			displacement = problem.solveElasticityProblemFT(Fext=Fext)[0]
			error_list[j] = problem.L2NormOfError_withExactFun(exactDisplacement_infPlate, displacement)

		nbctrlpts_list = (2**cuts_list+degree)**2
		ax.loglog(nbctrlpts_list, error_list, marker=markerSet[i], label='degree '+r'$p=\,$'+str(degree))

		if str(quadArgs['quadrule']) == 'iga':
			slope = np.polyfit(np.log10(nbctrlpts_list[:4]),np.log10(error_list[:4]), 1)[0]
			slope = round(slope, 1)
			annotation.slope_marker((nbctrlpts_list[3], error_list[3]), slope, 
									poly_kwargs={'facecolor': (0.73, 0.8, 1)})
			
		ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{L_2(\Omega)}}{||u||_{L_2(\Omega)}}$')
		ax.set_xlabel('Total number of DOF')
		ax.set_ylim(top=1e0, bottom=1e-15)
		ax.set_xlim(left=10, right=1e5)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigInfinitePlate2_' + str(quadArgs['quadrule']) +'.png')