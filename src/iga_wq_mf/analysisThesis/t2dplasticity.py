"""
.. Test of elastoplasticity 2D
.. We test how elastoplasticity module works
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import *
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem
import pickle
from pyevtk.vtk import VtkGroup

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

def run(folder=None):
	assert folder is not None, 'Folder unknown'
	print("Running group...")
	g = VtkGroup(folder)
	for i in range(20):
		g.addFile(filepath = folder + "pls"+str(i)+".vts", sim_time = i)
	g.save()

# Set global variables
GEONAME = 'SQ'
TRACTION = 150.0
YOUNG, POISSON = 2500, 0.25
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':5, 'poisson_ratio': POISSON, 
			'isoHardLaw': {'name':'linear', 'Eiso':0.0}, 'kineHardLaw':{'parameters':np.array([[500, 0]])}}
isReference = True

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	tmp = np.zeros((2, nnz)); tmp[1, :] = x**2-1/4
	F = np.zeros((2, nnz))
	F[1, :] = -TRACTION*(np.min(tmp, axis=0))**2
	return F

def simulate(degree, cuts, quadArgs, step=-2):
	geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'XY':np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])}
			}
	blockPrint()
	material = mechamat(MATARGS)
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	meshparam = modelPhy.compute_global_mesh_parameter()

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int); table[1, 0, :] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	problem._thresNL = 1e-6; problem._itersNL = 100
	Fref = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=3)[0]
	Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for k in range(len(TIME_LIST)): Fext_list[:, :, k] = np.sin(TIME_LIST[k])*Fref
	displacement = np.zeros(np.shape(Fext_list))
	_, internalVars = problem.solveElastoPlasticityProblem(displacement, Fext_list[:, :, :step+1])
	return problem, displacement[:, :, :step+1], meshparam, internalVars

if isReference:
	degree, cuts = 1, 6
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	problem, displacement, _, internalVars = simulate(degree, cuts, quadArgs)
	np.save(folder + 'disppl', displacement)
	with open(folder + 'refpartpl.pkl', 'wb') as outp:
		pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)
	stress_qp = internalVars.get('stress', None)
	for j, i in enumerate(range(0, NBSTEPS-1, 4)):
		devstress_qp = computeDeviatoric4All(stress_qp[:, :, i], dim=2)
		vonMises_qp = np.sqrt(3/2)*computeSymTensorNorm4All(devstress_qp, dim=2)
		vonMises_cp = problem.L2projectionCtrlpts(vonMises_qp)
		problem.part.exportResultsCP(fields={'stress': vonMises_cp}, name='pls'+str(j), folder=folder)
	run(folder=folder)
	
else:

	degree_list = np.array([1, 2, 3])
	cuts_list  = np.arange(2, 6)
	step_list  = range(0, NBSTEPS, 4)
	error_list = np.ones((len(step_list), len(degree_list), len(cuts_list)))
	mesh_list  = np.ones((len(degree_list), len(cuts_list)))

	with open(folder + 'refpartpl.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)
	disp_ref = np.load(folder + 'disppl.npy')
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	
	for i, degree in enumerate(degree_list):
		for j, cuts in enumerate(cuts_list):
			problem, displacement, meshparam, _= simulate(degree, cuts, quadArgs)
			mesh_list[i, j] = meshparam

			for k, step in enumerate(step_list):
				error_list[k, i, j], _ = problem.normOfError(displacement[:, :, step], 
														normArgs={'type':'H1', 
														'part_ref':part_ref, 
														'u_ref': disp_ref[:, :, step]})

	np.save(folder + 'plasticity2D', error_list)
	np.save(folder + 'meshpar2D', mesh_list)

	error_list = np.load(folder + 'plasticity2D.npy')
	mesh_list = np.load(folder+'meshpar2D.npy')

	for k, step in enumerate(step_list):
		fig, ax = plt.subplots(figsize=(9,6))
		for i, degree in enumerate(degree_list):
			color = COLORLIST[i]
			ax.loglog(mesh_list[i, :], error_list[k, i, :], color=color, marker='o', markerfacecolor='w',
						markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))
		ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{H_1(\Omega)}}{||u||_{H_1(\Omega)}}$')
		ax.set_xlabel('Mesh parameter ' + r'$h_{max}$')
		ax.set_ylim(top=1, bottom=1e-10)
		ax.set_xlim(left=5e-2, right=1e0)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigConvergencePlasticityAllH1' + str(step) +  GEONAME + '.pdf')
