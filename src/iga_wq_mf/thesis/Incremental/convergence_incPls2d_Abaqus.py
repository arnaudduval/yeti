"""
.. Test of elastoplasticity 2D
.. We test how elastoplasticity module works
.. Joaquin Cornejo 
"""

from thesis.Incremental.__init__ import *
from pysrc.lib.lib_material import computeDeviatoric4All, computeSymTensorNorm4All

FOLDER2SAVE = os.path.dirname(os.path.realpath(__file__)) + '/abaqus/'
if not os.path.isdir(FOLDER2SAVE): os.mkdir(FOLDER2SAVE)

FOLDER2DATA = FOLDER2SAVE + '/datafromsimu/'
if not os.path.isdir(FOLDER2DATA): os.mkdir(FOLDER2DATA)

NBSTEPS = 101
TIME_LIST = np.linspace(0, 1, NBSTEPS)
MATARGS = {'elastic_modulus':1e5, 
			'elastic_limit':100, 
			'poisson_ratio':0.3, 
			'isoHardLaw': {'name':'linear', 'Eiso':1e4}, 
		}

degList = np.arange(1, 4)
cutList = np.arange(1, 7)
stepList = np.arange(1, NBSTEPS, 4)
RUNSIMU = False

def forceSurf_abaqus(P:list):
	x = P[0, :]; nnz = np.size(P, axis=1)
	F = np.zeros((2, nnz))
	F[0, :] = 200
	return F

def simulate_2d_abaqus(degree, cuts, quadArgs):

	blockPrint()
	material = mechamat(MATARGS)
	modelGeo = Geomdl({'name': 'abaqus'})
	modelIGA = modelGeo.getIGAParametrization()
	modelIGA.refine(nb_degreeElevationByDirection=np.array([degree-2, degree-1, 0]),
					nb_refinementByDirection=cuts*np.ones(3, dtype=int),)
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[0, 0, 0] = 1; table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	Fref = problem.compute_surfForce(forceSurf_abaqus, nbFacePosition=1)[0]
	FextList = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for i, k in enumerate(TIME_LIST): FextList[:, :, i] = k*Fref
	displacement = np.zeros(np.shape(FextList))
	resLin, internalVars = problem.solveElastoPlasticityProblem(displacement, FextList)
	return problem, displacement, resLin, internalVars

if RUNSIMU:
	degree, cuts = 3, 5
	quadArgs = {'quadrule': 'wq', 'type': 2}
	problem, displacement, _, internalVars = simulate_2d_abaqus(degree, cuts, quadArgs)
	np.save(FOLDER2DATA + 'disppl2d', displacement)
	with open(FOLDER2DATA + 'refpartpl2d.pkl', 'wb') as outp:
		pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)

	stress_qp = internalVars.get('stress', None)
	plseq_qp = internalVars.get('plseq', None)
	plastic_qp = np.where(np.abs(plseq_qp)<1e-8, 0.0, 1.0)
	filename = 'out_'
	for j, i in enumerate(stepList):
		devstress_qp = computeDeviatoric4All(stress_qp[:, :, i])
		vonMises_qp = np.sqrt(3/2)*computeSymTensorNorm4All(devstress_qp)
		problem.part.postProcessingDual(name=filename+str(j), folder=FOLDER2DATA, 
										fields={'stress': vonMises_qp, 		
												'straineq': plseq_qp[0, :, i], 
												'plastic': plastic_qp[0, :, i]})
	run(folder=FOLDER2DATA, filename=filename, nbFiles=len(stepList))