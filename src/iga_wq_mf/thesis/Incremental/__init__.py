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

FOLDER2DATA = FOLDER2SAVE + '/datafromsimu/'
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

def buildmatrix_ht(problem:heatproblem, prototype=True):
	args = {'position':problem.part.qpPhy}

	quadrules = problem.part._quadraturerules
	matrix = sp.csr_matrix((problem.part.nbctrlpts_total, problem.part.nbctrlpts_total))
	if prototype:
		submatrices = []
		for j in range(problem.part.dim):
			submatrices.append(quadrules[0]._denseWeights[0] @ quadrules[0]._denseBasis[0].T)
		if problem.part.dim == 2: matrix = sp.kron(submatrices[1], submatrices[0])
		elif problem.part.dim == 3:	matrix = sp.kron(submatrices[2], sp.kron(submatrices[1], submatrices[0]))
	else:
		prop = np.einsum('ilk,jmk,lmk,k->ijk', problem.part.invJ, problem.part.invJ,
				problem.heatmaterial.conductivity(args), problem.part.detJ)
	
		for j in range(problem.part.dim):
			beta = np.zeros(problem.part.dim, dtype=int); beta[j] = 1
			if problem.part.dim == 2: tmp1 = sp.kron(quadrules[1]._denseBasis[beta[1]], quadrules[0]._denseBasis[beta[0]]).T
			elif problem.part.dim == 3: tmp1 = sp.kron(quadrules[2]._denseBasis[beta[2]], sp.kron(quadrules[1]._denseBasis[beta[1]], quadrules[0]._denseBasis[beta[0]])).T
			for i in range(problem.part.dim):
				alpha = np.zeros(problem.part.dim, dtype=int); alpha[i] = 1
				zeta = beta + 2*alpha
				tmp2 = sp.diags(prop[i, j, :]) @ tmp1
				if problem.part.dim == 2: tmp3 = sp.kron(quadrules[1]._denseWeights[zeta[1]], quadrules[0]._denseWeights[zeta[0]]) @ tmp2
				elif problem.part.dim == 3: tmp3 = sp.kron(quadrules[2]._denseWeights[zeta[2]], sp.kron(quadrules[1]._denseWeights[zeta[1]], quadrules[0]._denseWeights[zeta[0]])) @ tmp2
				matrix += tmp3
	
	return matrix 