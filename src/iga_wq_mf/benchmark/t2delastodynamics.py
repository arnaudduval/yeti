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
YOUNG, POISSON = 1e3, 0.3
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
GEONAME = 'QA'
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e10, 'poisson_ratio':POISSON,
		'plasticLaw': {'Isoname':'none'}}
SOLVERARGS  = {'nbIterationsPCG':150, 'PCGThreshold':1e-15, 'PCGmethod': 'TDC'}

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = RINT**2/r_square
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = TRACTION/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = TRACTION/2*3*np.sin(3*theta)*(b**2 - b)
	return F

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
	problem.addSolverConstraints(solverArgs=SOLVERARGS)
	Fref = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
	Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for k in range(len(TIME_LIST)): Fext_list[:, :, k] = np.sin(TIME_LIST[k])*Fref
	displacement = np.zeros(np.shape(Fext_list))
	problem.solveDynamicsProblem(displacement, Fext_list, TIME_LIST, isLumped=False)
	return problem, displacement, meshparam

degree, cuts = 3, 4
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
problem, displacement, _ = simulate(degree, cuts, quadArgs)
