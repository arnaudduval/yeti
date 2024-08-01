"""
.. Test of elastoplasticity 2D
.. We test how elastoplasticity module works
.. Joaquin Cornejo 
"""

from thesis.Incremental.__init__ import *
from pysrc.lib.lib_material import computeDeviatoric4All, computeSymTensorNorm4All
from pysrc.lib.lib_base import vtk2png

FOLDER2SAVE = os.path.dirname(os.path.realpath(__file__)) + '/abaqus/'
if not os.path.isdir(FOLDER2SAVE): os.mkdir(FOLDER2SAVE)

FOLDER2DATA = FOLDER2SAVE + '/datafromsimu_el/'
if not os.path.isdir(FOLDER2DATA): os.mkdir(FOLDER2DATA)

modelGeo = Geomdl({'name': 'abaqus'})
modelIGA = modelGeo.getIGAParametrization()
modelIGA.refine(nb_degreeElevationByDirection=np.array([0, 0, 0]),
				nb_refinementByDirection=4*np.ones(3, dtype=int),)
quadArgs = {'quadrule': 'wq', 'type': 2}
modelPhy = part(modelIGA, quadArgs=quadArgs)
COORSIGA = modelPhy.interpolateMeshgridField(sampleSize=201)[0]

filename = FOLDER2SAVE + 'job.dat'
COORDSABAQUS = np.loadtxt(filename, delimiter=',')

x = np.linspace(0,50,201)
y = x**2/250+10

fig, ax = plt.subplots()
ax.scatter(COORDSABAQUS[:, 1], COORDSABAQUS[:, 2], s=0.2, label='FEA')
ax.scatter(COORSIGA[0, :], COORSIGA[1, :], s=0.2, label='IGA')
ax.plot(x, y, color='k', label='Parabola')
ax.legend(loc='upper left')
ax.grid(False)
fig.savefig(FOLDER2SAVE+'Abaqus')

NBSTEPS = 2
TIME_LIST = np.linspace(0, 1, NBSTEPS)
MATARGS = {'elastic_modulus':2e5, 
			'elastic_limit':1e6, 
			'poisson_ratio':0.3, 
			'isoHardLaw': {'name':'linear', 'Eiso':1e5}, 
		}

stepList = np.arange(1, NBSTEPS, 4)
RUNSIMU = True

def forceSurf_abaqus(P:list):
	x = P[0, :]; nnz = np.size(P, axis=1)
	F = np.zeros((2, nnz))
	F[0, :] = 200
	return F

def simulate_2d_abaqus(degree, cuts, quadArgs):

	material = mechamat(MATARGS)
	modelGeo = Geomdl({'name': 'abaqus'})
	modelIGA = modelGeo.getIGAParametrization()
	modelIGA.refine(nb_degreeElevationByDirection=np.array([degree-2, degree-1, 0]),
					nb_refinementByDirection=np.array([cuts,cuts-1,0]))
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[0, 0, 0] = 1; table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	Fref = problem.compute_surfForce(forceSurf_abaqus, nbFacePosition=1)[0]
	FextList = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for i, k in enumerate(TIME_LIST): FextList[:, :, i] = k*Fref
	displacement = np.zeros(np.shape(FextList)); resLin = None; internalVars = None
	# resLin, internalVars = problem.solveElastoPlasticityProblem(displacement, FextList)
	displacement = problem._solveLinearizedElasticityProblem(FextList[:,:,-1])[0]

	return problem, displacement, resLin, internalVars

if RUNSIMU:
	quadArgs = {'quadrule': 'wq', 'type': 2}
	for degree in range(2, 6):
		for cuts in range(2, 8):
			blockPrint()
			problem, displacement, _, _ = simulate_2d_abaqus(degree, cuts, quadArgs)
			disp_interp = problem.part.interpolateMeshgridField(displacement, sampleSize=101)[-1]
			enablePrint()
			print('degree: %d, nbDOF:%d, max:%.10e, min:%.10e' %(degree, problem.part.nbctrlpts_total, np.max(disp_interp), np.min(disp_interp)))

# 	degree, cuts = 2, 4
# 	quadArgs = {'quadrule': 'wq', 'type': 2}
# 	problem, displacement, _, internalVars = simulate_2d_abaqus(degree, cuts, quadArgs)
# 	np.save(FOLDER2DATA + 'disppl2d', displacement)
# 	with open(FOLDER2DATA + 'refpartpl2d.pkl', 'wb') as outp:
# 		pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)

# 	filename = 'out_'
# 	problem.part.postProcessingPrimal(name=filename+'last', 
# 										folder=FOLDER2DATA,
# 										sampleSize=201,
# 										fields={'disp1': displacement[0, :, -1], 		
# 												'disp2': displacement[1, :, -1], 
# 												})

# 	stress_qp = internalVars.get('stress', None)
# 	plseq_qp = internalVars.get('plseq', None)
# 	plastic_qp = np.where(np.abs(plseq_qp)<1e-8, 0.0, 1.0)
# 	filename = 'out_'
# 	for j, i in enumerate(stepList):
# 		devstress_qp = computeDeviatoric4All(stress_qp[:, :, i])
# 		vonMises_qp = np.sqrt(3/2)*computeSymTensorNorm4All(devstress_qp)
# 		problem.part.postProcessingDual(name=filename+str(j), folder=FOLDER2DATA, 
# 										fields={'stress': vonMises_qp, 		
# 												'straineq': plseq_qp[0, :, i], 
# 												'plastic': plastic_qp[0, :, i]})
# 	run(folder=FOLDER2DATA, filename=filename, nbFiles=len(stepList))

# vtk2png(FOLDER2DATA, filename='out_last', fieldname='disp1', title='Displacement X', position_y=0.1, n_colors=120, fmt='%.2e')
# vtk2png(FOLDER2DATA, filename='out_last', fieldname='disp2', title='Displacement Y', position_y=0.1, n_colors=120, fmt='%.2e')
# vtk2png(FOLDER2DATA, filename='out_24', fieldname='stress', title='Von Mises stress', position_y=0.1, n_colors=120)
# vtk2png(FOLDER2DATA, filename='out_24', fieldname='straineq', title='Equivalent plastic strain', position_y=0.1, n_colors=120, fmt='%.2e')
# vtk2png(FOLDER2DATA, filename='out_24', fieldname='plastic', title='Plastic zone', position_y=0.1, n_colors=2, n_labels=2)