from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part1D, part
from pysrc.lib.lib_job1d import mechaproblem1D, heatproblem1D, problem1D
from pysrc.lib.lib_job3d import mechaproblem, heatproblem
from pysrc.lib.lib_material import mechamat, heatmat
from pysrc.lib.lib_boundary import boundaryCondition
import pickle
from pyevtk.vtk import VtkGroup
from pysrc.lib.lib_base import cropImage

FOLDER2SAVE = os.path.dirname(os.path.realpath(__file__)) + '/results/'
if not os.path.isdir(FOLDER2SAVE): os.mkdir(FOLDER2SAVE)

FOLDER2DATA = os.path.dirname(os.path.realpath(__file__)) + '/datafromsimu/'
if not os.path.isdir(FOLDER2DATA): os.mkdir(FOLDER2DATA)

def run(filename=None, folder=None, nbFiles=1):
	assert folder is not None, 'Folder unknown'
	if filename is None: filename = 'out'
	print("Running group...")
	g = VtkGroup(folder)
	for i in range(nbFiles):
		g.addFile(filepath = folder + filename + str(i) + '.vts', sim_time = i)
	g.save()
	return

TRACTION = 400.0
YOUNG, POISSON = 2500, 0.25
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi/2, NBSTEPS)
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':5, 'poisson_ratio': POISSON, 
			'isoHardLaw': {'name':'linear', 'Eiso':0.0}, 
			'kineHardLaw':{'parameters':np.array([[500, 0]])}
			}

degList = np.arange(1, 4)
cutList = np.arange(1, 7)
stepList = np.arange(1, NBSTEPS, 4)
RUNSIMU = False

def forceSurf(P:list):
	x = P[0, :]; nnz = np.size(P, axis=1)
	tmp = np.zeros((2, nnz)); tmp[1, :] = x**2-1/4
	F = np.zeros((2, nnz))
	F[1, :] = -TRACTION*(np.min(tmp, axis=0))**2
	return F

def simulate_2d(degree, cuts, quadArgs, precond='JMC'):
	geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			}
	blockPrint()
	material = mechamat(MATARGS)
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[0, 0, 0] = 1; table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary); problem._linPreCond = precond
	Fref = problem.compute_surfForce(forceSurf, nbFacePosition=3)[0]
	FextList = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for k in range(len(TIME_LIST)): FextList[:, :, k] = np.sin(TIME_LIST[k])*Fref
	displacement = np.zeros(np.shape(FextList))
	resLin, internalVars = problem.solveElastoPlasticityProblem(displacement, FextList)
	return problem, displacement, resLin, internalVars
