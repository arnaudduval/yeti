"""
.. Test of elastoplasticity 2D
.. We test how elastoplasticity module works
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem
import pickle

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
TRACTION, RINT = 1.0, 1.0
YOUNG, POISSON = 1e3, 0.3
NBSTEPS = 15
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
GEONAME = 'QA'
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':2, 'poisson_ratio': POISSON, 
			'plasticLaw': {'Isoname':'linear', 'Eiso':YOUNG/10}}
SOLVERARGS = {'nbIterationsPCG':150, 'PCGThreshold':1e-10, 'PCGmethod': 'TDC'}
isReference = True

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = RINT**2/r_square # Already squared
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = TRACTION/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = TRACTION/2*3*np.sin(3*theta)*(b**2 - b)
	return F

def simulate(degree, cuts, quadArgs, step=-2):
	geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
							'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
							'extra':{'Rin':1.0, 'Rex':4.0}
			}
	blockPrint()
	material = mechamat(MATARGS)
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	meshparam = modelPhy.compute_mesh_parameter()

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[1, 1, 0] = 1
	table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	problem.addSolverConstraints(solverArgs=SOLVERARGS)
	Fref = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
	Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for k in range(len(TIME_LIST)): Fext_list[:, :, k] = np.sin(TIME_LIST[k])*Fref
	displacement = problem.solvePlasticityProblemPy(Fext_list=Fext_list[:, :, :step+1])[0]
	
	return problem, displacement, meshparam

if isReference:
	degree, cuts = 3, 10 
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	problem, displacement, _ = simulate(degree, cuts, quadArgs)
	np.save(folder + 'disppl', displacement)
	with open(folder + 'refpartpl.pkl', 'wb') as outp:
		pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)

else:

	# Plot results
	normalPlot  = {'marker': 'o', 'linestyle': '-', 'markersize': 10}
	onlyMarker1 = {'marker': '.', 'linestyle': ':', 'markersize': 6}
	onlyMarker2 = {'marker': 'x', 'linestyle': 'None', 'markersize': 6}

	degree_list = np.array([2, 3, 4, 6])
	cuts_list   = np.arange(2, 7)

	with open(folder + 'refpartpl.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)
	disp_ref = np.load(folder + 'disppl.npy')

	fig, ax    = plt.subplots(figsize=(8, 7))
	for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 1, 2], [normalPlot, onlyMarker1, onlyMarker2]):
		quadArgs = {'quadrule': quadrule, 'type': quadtype}
		error_list = np.ones(len(cuts_list))
		
		for i, degree in enumerate(degree_list):
			meshparam = np.ones(len(cuts_list))
			color = COLORLIST[i]
			for j, cuts in enumerate(cuts_list):
				step = 12
				problem, displacement, meshparam[j] = simulate(degree, cuts, quadArgs, step)
				error_list[j] = problem.normOfError(displacement[:, :, step], 
								normArgs={'referencePart':part_ref, 'u_ref': disp_ref[:, :, step]})

			ax.loglog(meshparam, error_list, color=color, marker=plotpars['marker'], markerfacecolor='w',
						markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{H_1(\Omega)}}{||u||_{H_1(\Omega)}}$')
			ax.set_xlabel('Mesh parameter ' + r'$h_{max}$')
			ax.set_ylim(top=1, bottom=1e-14)
			ax.set_xlim(left=1e-2, right=2)
			fig.tight_layout()
			fig.savefig(folder + 'FigConvergencePlasticityAllH1' +  GEONAME + '.pdf')
		