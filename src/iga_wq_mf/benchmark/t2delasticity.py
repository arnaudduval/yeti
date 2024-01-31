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
import pickle

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
TRACTION, RINT, REXT = 1.0, 1.0, 2.0
YOUNG, POISSON = 1e3, 0.0
GEONAME = 'QA'
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e10, 'poisson_ratio':POISSON,
		'isoHardLaw': {'name':'none'}}
isReference = False

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = RINT**2/r_square
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = TRACTION/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = TRACTION/2*3*np.sin(3*theta)*(b**2 - b)
	return F

def exactDisplacement_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	theta = np.arcsin(y/np.sqrt(r_square))
	b = RINT**2/r_square # Already squared
	c = TRACTION*(1.0 + POISSON)*np.sqrt(r_square)/(2*YOUNG)

	disp = np.zeros((2, nnz))
	disp[0, :] = c*(2*(1-POISSON)*np.cos(theta) + b*(4*(1-POISSON)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	disp[1, :] = c*(-2*POISSON*np.sin(theta) + b*(2*(-1 + 2*POISSON)*np.sin(theta) + np.sin(3*theta)) - b**2*np.sin(3*theta))
	
	return disp

def simulate(degree, cuts, quadArgs, useElastoAlgo=False):
	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
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
	if useElastoAlgo:
		Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, 2))
		Fext_list[:, :, 1] = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
		tmp = np.zeros(np.shape(Fext_list))
		problem.solvePlasticityProblem(tmp, Fext_list)
		displacement = tmp[:, :, -1]
	else:
		Fext = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
		displacement = problem.solveElasticityProblem(Fext)[0]
	return problem, displacement, meshparam

if isReference:
	degree, cuts = 7, 8
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	problem, displacement, _ = simulate(degree, cuts, quadArgs)
	np.save(folder + 'dispel', displacement)
	with open(folder + 'refpartel.pkl', 'wb') as outp:
		pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)

else:	

	# Plot results
	normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
	onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
	onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}

	degree_list = np.array([1, 2, 3, 4])
	cuts_list   = np.arange(1, 8)

	disp_ref = np.load(folder + 'dispel.npy')
	with open(folder + 'refpartel.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)

	fig, ax = plt.subplots(figsize=(8, 7))
	figname = folder + 'FigElasLinearConvergenceAllL2_nu0' + '.pdf'
	# for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 1, 2], [normalPlot, onlyMarker1, onlyMarker2]):
	for quadrule, quadtype, plotpars in zip(['iga', 'wq'], ['leg', 1], [normalPlot, onlyMarker1]):
		quadArgs = {'quadrule': quadrule, 'type': quadtype}
		error_list = np.ones(len(cuts_list))

		for i, degree in enumerate(degree_list):
			meshparam = np.ones(len(cuts_list))
			color = COLORLIST[i]
			for j, cuts in enumerate(cuts_list):
				problem, displacement, meshparam[j] = simulate(degree, cuts, quadArgs, useElastoAlgo=False)
				# error_list[j], _ = problem.normOfError(displacement, normArgs={'type':'H1', 
				#										'part_ref':part_ref, 'u_ref':disp_ref})
				error_list[j], _ = problem.normOfError(displacement, normArgs={'type':'L2', 
														'exactFunction':exactDisplacement_infPlate})
			if quadrule == 'iga': 
				ax.loglog(meshparam, error_list, label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
						markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			else: 
				ax.loglog(meshparam, error_list, color=color, marker=plotpars['marker'], markerfacecolor='w',
					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			
			fig.savefig(figname)

ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ 4")
# ax.loglog([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
# 		markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label="IGA-WQ 2")

ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Omega)}$')
# ax.set_ylabel(r'$\displaystyle ||u - u^h||_{H^1(\Omega)}$')
ax.set_xlabel('Mesh parameter ' + r'$h_{max}$')
ax.set_ylim(top=1e-2, bottom=1e-14)
ax.set_xlim(left=1e-2, right=5)
ax.legend()
fig.tight_layout()
fig.savefig(figname)