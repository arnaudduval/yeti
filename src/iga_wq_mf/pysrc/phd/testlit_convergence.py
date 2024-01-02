from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder2save = os.path.dirname(full_path) + '/results/biblio/'
folder2find = os.path.dirname(full_path) + '/data/'

# Set global variables
degree_list = np.array([2, 4, 6, 8])
cuts_list   = np.array([i for i in range(1, 6)])
fig, ax = plt.subplots(figsize=(8, 6))

# Set global variables
TRACTION, RINT, REXT = 10.0, 1.0, 4.0
YOUNG, POISSON = 1e5, 0.3
GEONAME = 'QA'
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e10, 'poisson_ratio':POISSON,
		'plasticLaw': {'Isoname':'none'}}
SOLVERARGS  = {'nIterKrylov':150, 'thresholdKrylov':1e-15, 'KrylovPreconditioner': 'TDC'}

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

	# disp = np.zeros((2, nnz))
	# disp[0, :] = c*(2*(1-POISSON)*np.cos(theta) + b*(4*(1-POISSON)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	# disp[1, :] = c*(-2*POISSON*np.sin(theta) + b*(2*(-1 + 2*POISSON)*np.sin(theta) + np.sin(3*theta)) - b**2*np.sin(3*theta))
	disp = c*(2*(1-POISSON)*np.cos(theta) + b*(4*(1-POISSON)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	
	return disp

def simulate(degree, cuts, quadArgs):
	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':RINT, 'Rex':REXT}
				}
	blockPrint()
	material = mechamat(MATARGS)
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[1, 1, 0] = 1; table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	problem.addSolverConstraints(solverArgs=SOLVERARGS)
	Fext = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
	displacement = problem.solveElasticityProblem(Fext=Fext)[0]
	return problem, displacement

quadArgs = {'quadrule': 'wq', 'type': 2}
error_list = np.ones(len(cuts_list))

for i, degree in enumerate(degree_list):
	color = COLORLIST[i]
	for j, cuts in enumerate(cuts_list):
		problem, displacement = simulate(degree, cuts, quadArgs)
		error_list[j] = problem.normOfError(displacement[0, :], normArgs={'type':'L2', 'exactFunction':exactDisplacement_infPlate})
	ax.loglog(2**cuts_list, error_list, color=color, marker='.', linestyle='--')

# Load data
file = pd.read_table(folder2find + 'elasticconverg'   + '.dat', sep='\t', names=['nbel', 'p2', 'p4', 'p6', 'p8', 'p10']) 
nbel = file.nbel
error = [file.p2, file.p4, file.p6, file.p8]

for i, deg in enumerate(degree_list):
	color = COLORLIST[i]
	ax.loglog(nbel, error[i], label='degree '+r'$p=\,$' + str(deg), marker='s', color=color)

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel('Number of elements')
ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{L_2(\Omega)}}{||u||_{L_2(\Omega)}}$')
ax.set_xlim(left=1e0, right=80)
ax.set_ylim(bottom=1e-12, top=1e0)
fig.tight_layout()
fig.savefig(folder2save + 'elasticerror' + '.pdf')