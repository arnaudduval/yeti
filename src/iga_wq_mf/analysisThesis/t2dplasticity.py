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
TRACTION = 400.0
YOUNG, POISSON = 2500, 0.25
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi/2, NBSTEPS)
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':5, 'poisson_ratio': POISSON, 
			'isoHardLaw': {'name':'linear', 'Eiso':0.0}, 
			'kineHardLaw':{'parameters':np.array([[500, 0]])}
			}
isReference = False

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	tmp = np.zeros((2, nnz)); tmp[1, :] = x**2-1/4
	F = np.zeros((2, nnz))
	F[1, :] = -TRACTION*(np.min(tmp, axis=0))**2
	return F

def simulate(degree, cuts, quadArgs):
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
	_, internalVars = problem.solveElastoPlasticityProblem(displacement, Fext_list)
	return problem, displacement, meshparam, internalVars

if isReference:
	degree, cuts = 2, 8
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	problem, displacement, _, internalVars = simulate(degree, cuts, quadArgs)
	# np.save(folder + 'disppl', displacement)
	# with open(folder + 'refpartpl.pkl', 'wb') as outp:
	# 	pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)

	stress_qp = internalVars.get('stress', None)
	alpha_qp = internalVars.get('hardening', None)
	plastic_qp = np.where(np.abs(alpha_qp)<1e-6, 0.0, 1.0)
	for j, i in enumerate(range(0, NBSTEPS, 4)):
		devstress_qp = computeDeviatoric4All(stress_qp[:, :, i], dim=2)
		vonMises_qp = np.sqrt(3/2)*computeSymTensorNorm4All(devstress_qp, dim=2)
		vonMises_cp = problem.L2projectionCtrlpts(vonMises_qp)
		alpha_cp = problem.L2projectionCtrlpts(alpha_qp[0, :, i])
		plastic_cp = problem.L2projectionCtrlpts(plastic_qp[0, :, i])
		problem.part.exportResultsCP(fields={'stress': vonMises_cp, 'straineq': alpha_cp, 'plastic':plastic_cp}, 
									name='pls'+str(j), folder=folder, sampleSize=401)
	run(folder=folder)
	
else:
	
	with open(folder + 'refpartpl.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)
	disp_ref = np.load(folder + 'disppl.npy')

	normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
	onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}

	degree_list = np.array([1, 2, 3, 4])
	cuts_list  = np.arange(2, 7)
	step_list  = range(1, NBSTEPS, 3)

	# # quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	# quadArgs = {'quadrule': 'wq', 'type': 1}
	# error_list = np.ones((len(step_list), len(degree_list), len(cuts_list)))
	# mesh_list  = np.ones((len(degree_list), len(cuts_list)))
	# for i, degree in enumerate(degree_list):
	# 	for j, cuts in enumerate(cuts_list):
	# 		problem, displacement, meshparam, _= simulate(degree, cuts, quadArgs)
	# 		mesh_list[i, j] = meshparam

	# 		for k, step in enumerate(step_list):
	# 			error_list[k, i, j], _ = problem.normOfError(displacement[:, :, step], 
	# 													normArgs={'type':'H1', 
	# 													'part_ref':part_ref, 
	# 													'u_ref': disp_ref[:, :, step]})

	# quadrule = quadArgs['quadrule']; quadtype = quadArgs['type']
	# np.save(folder + 'plasticity2D'+quadrule+str(quadtype), error_list)
	# np.save(folder + 'meshpar2D'+quadrule+str(quadtype), mesh_list)

	for k, step in enumerate(step_list):
		fig, ax = plt.subplots(figsize=(6, 5))
		for quadrule, quadtype, plotpars in zip(['iga', 'wq'], ['leg', 1], [normalPlot, onlyMarker1]):
			error_list = np.ones(len(cuts_list))
			error_list = np.load(folder + 'plasticity2D'+quadrule+str(quadtype)+'.npy')
			mesh_list = np.load(folder+'meshpar2D'+quadrule+str(quadtype)+'.npy')
		
			for i, degree in enumerate(degree_list):
				color = COLORLIST[i]
				if quadrule == 'iga': 
					ax.loglog(mesh_list[i, :], error_list[k, i, :], label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
							markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
					# slope = round(np.polyfit(np.log(mesh_list[i, 2:]), np.log(error_list[k, i, 2:]), 1)[0], 1)
					# annotation.slope_marker((mesh_list[i, -1],  error_list[k, i, -1]), slope, 
					# 				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
				else: 
					ax.loglog(mesh_list[i, :], error_list[k, i, :], color=color, marker=plotpars['marker'], markerfacecolor='w',
						markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])

		ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ 4")

		ax.set_ylabel(r'$\displaystyle ||u - u^h||_{H^1(\Omega)}$')
		ax.set_xlabel('Mesh parameter ' + r'$h_{max}$')
		ax.set_ylim(top=1, bottom=1e-8)
		ax.set_xlim(left=2e-2, right=1e0)
		ax.legend(loc='lower right')
		fig.tight_layout()
		fig.savefig(folder + 'FigConvergencePlasticity' + str(step) +  GEONAME + '.png')
		plt.close(fig=fig)
