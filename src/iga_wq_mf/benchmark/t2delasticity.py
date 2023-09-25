"""
.. Test of elasticity 2D
.. Infinite plate with a hole under uniaxial traction. 
.. The analytical solution of this problem is given by Timoshenko
.. The convergence curves are traced for IGA-Legendre and IGA-WQ
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
folder = os.path.dirname(full_path) + '/results/d2elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceSurf_infPlate(P:list):
	Tx, a = 1.0, 1.0
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = a**2/r_square
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = Tx/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = Tx/2*3*np.sin(3*theta)*(b**2 - b)
	return F

def exactDisplacement_infPlate(P:list):
	Tx, a, E, nu = 1.0, 1.0, 1.e3, 0.0
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	theta = np.arcsin(y/np.sqrt(r_square))
	b = a**2/r_square # Already squared
	c = Tx*(1.0 + nu)*np.sqrt(r_square)/(2*E)

	disp = np.zeros((2, nnz))
	disp[0, :] = c*(2*(1-nu)*np.cos(theta) + b*(4*(1-nu)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	disp[1, :] = c*(-2*nu*np.sin(theta) + b*(2*(-1 + 2*nu)*np.sin(theta) + np.sin(3*theta)) - b**2*np.sin(3*theta))
	
	return disp

# Set global variables
geoName = 'QA'
E, nu = 1e3, 0.0
matArgs    = {'elastic_modulus':E, 'elastic_limit':1e10, 'poisson_ratio':nu}
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-15, 'PCGmethod': 'TDC'}
degree_list = np.array([2, 3, 4, 6, 8])
cuts_list   = np.arange(2, 9)

# Plot results
normalPlot = {'marker': 'o', 'linestyle': '-', 'markersize': 10}
onlyMarker1 = {'marker': '.', 'linestyle': ':', 'markersize': 6}
onlyMarker2 = {'marker': 'x', 'linestyle': 'None', 'markersize': 6}

fig, ax = plt.subplots(figsize=(8, 7))
for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 1, 2], [normalPlot, onlyMarker1, onlyMarker2]):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	meshparam  = np.ones(len(cuts_list))
	error_list = np.ones(len(cuts_list))
	time_list  = np.zeros((len(degree_list), len(cuts_list)))

	for i, degree in enumerate(degree_list):
		color = colorList[i]
		for j, cuts in enumerate(cuts_list):
			geoArgs = {'name': geoName, 'degree': degree*np.ones(3, dtype=int), 
						'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
						'extra':{'Rin':1.0, 'Rex':4.0}
			}
			blockPrint()
			material = mechamat(matArgs)
			modelGeo = Geomdl(geoArgs)
			modelIGA = modelGeo.getIGAParametrization()
			modelPhy = part(modelIGA, quadArgs=quadArgs)
			meshparam[j] = modelPhy.compute_mesh_parameter()

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
			error_list[j] = problem.L2NormOfError(displacement, L2NormArgs={'exactFunction':exactDisplacement_infPlate})

		ax.loglog(meshparam, error_list, color=color, marker=plotpars['marker'], markerfacecolor='w',
					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
		ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{L_2(\Omega)}}{||u||_{L_2(\Omega)}}$')
		ax.set_xlabel('Mesh parameter ' + r'$h_{max}$')
		ax.set_ylim(top=1e-2, bottom=1e-14)
		ax.set_xlim(left=1e-2, right=2)
		fig.tight_layout()
		fig.savefig(folder + 'FigConvergenceAll2' +  geoName + '.pdf')
	