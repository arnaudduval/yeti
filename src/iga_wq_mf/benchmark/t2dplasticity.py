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
TRACTION, RINT, REXT = 1.0, 1.0, 2.0
YOUNG, POISSON = 1e3, 0.3
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
GEONAME = 'QA'
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1.5, 'poisson_ratio': POISSON, 
			'isoHardLaw': {'name':'linear', 'Eiso':YOUNG/10}}
isReference = False

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
				'extra':{'Rin':RINT, 'Rex':REXT}
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
	table[1, 1, 0] = 1; table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	Fref = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
	Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for k in range(len(TIME_LIST)): Fext_list[:, :, k] = np.sin(TIME_LIST[k])*Fref
	displacement = np.zeros(np.shape(Fext_list))
	_, internalVars = problem.solvePlasticityProblem(displacement, Fext_list[:, :, :step+1])
	return problem, displacement[:, :, :step+1], meshparam, internalVars

if isReference:
	degree, cuts = 2, 8
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	problem, displacement, _, internalVars = simulate(degree, cuts, quadArgs, step=50)
	np.save(folder + 'disppl', displacement)
	with open(folder + 'refpartpl.pkl', 'wb') as outp:
		pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)
	hardening_qp = internalVars.get('hardening', None)
	for i in range(30, 50, 3):
		hardening_cp = problem.L2projectionCtrlpts(hardening_qp[:, :, i])
		problem.part.exportResultsCP(fields={'hardening': hardening_cp}, name='pls'+str(i))

else:

	degree_list = np.array([1, 2, 3])
	cuts_list   = np.arange(2, 7)
	step_max    = 50
	step_list  = range(20, step_max, 2)
	error_list = np.ones((len(step_list), len(degree_list), len(cuts_list)))

	with open(folder + 'refpartpl.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)
	disp_ref = np.load(folder + 'disppl.npy')
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	
	for i, degree in enumerate(degree_list):
		for j, cuts in enumerate(cuts_list):
			problem, displacement, _, _= simulate(degree, cuts, quadArgs, step=step_max)

			for k, step in enumerate(step_list):
				error_list[k, i, j], _ = problem.normOfError(displacement[:, :, step], 
														normArgs={'typenorm':'H1', 
														'part_ref':part_ref, 
														'u_ref': disp_ref[:, :, step]})

	np.save(folder + 'plasticity2D', error_list)
	error_list = np.load(folder + 'plasticity2D.npy')

	for k, step in enumerate(step_list):
		fig, ax = plt.subplots(figsize=(9,6))
		for i, degree in enumerate(degree_list):
			color = COLORLIST[i]
			meshparam_list = np.zeros(len(cuts_list))
			for j, cuts in enumerate(cuts_list):
				geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
					'extra':{'Rin':RINT, 'Rex':REXT}
				}
				quadArgs = {'quadrule': 'iga', 'type': 'leg'}
				blockPrint()
				modelGeo = Geomdl(geoArgs)
				modelIGA = modelGeo.getIGAParametrization()
				modelPhy = part(modelIGA, quadArgs=quadArgs)
				meshparam = modelPhy.compute_mesh_parameter(); meshparam_list[j] = meshparam
				enablePrint()

			ax.loglog(meshparam_list, error_list[k, i, :], color=color, marker='o', markerfacecolor='w',
						markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))
		ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{H_1(\Omega)}}{||u||_{H_1(\Omega)}}$')
		ax.set_xlabel('Mesh parameter ' + r'$h_{max}$')
		ax.set_ylim(top=1, bottom=1e-10)
		ax.set_xlim(left=5e-2, right=1e0)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigConvergencePlasticityAllH1' + str(step) +  GEONAME + '.pdf')
