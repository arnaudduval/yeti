"""
.. Test of elastoplasticity 2D
.. We test how elastoplasticity module works
.. Joaquin Cornejo 
"""

from thesis.Incremental.__init__ import *
from pysrc.lib.lib_material import computeDeviatoric4All, computeSymTensorNorm4All

# Set global variables
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

def simulate(degree, cuts, quadArgs):
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
	problem = mechaproblem(material, modelPhy, boundary)
	Fref = problem.compute_surfForce(forceSurf, nbFacePosition=3)[0]
	FextList = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for k in range(len(TIME_LIST)): FextList[:, :, k] = np.sin(TIME_LIST[k])*Fref
	displacement = np.zeros(np.shape(FextList))
	_, internalVars = problem.solveElastoPlasticityProblem(displacement, FextList)
	return problem, displacement, internalVars

if RUNSIMU:
	degree, cuts = 6, 8
	quadArgs = {'quadrule': 'wq', 'type': 2}
	problem, displacement, internalVars = simulate(degree, cuts, quadArgs)
	np.save(FOLDER2SAVE + 'disppl2d', displacement)
	with open(FOLDER2SAVE + 'refpartpl2d.pkl', 'wb') as outp:
		pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)

	stress_qp = internalVars.get('stress', None)
	plseq_qp = internalVars.get('plseq', None)
	plastic_qp = np.where(np.abs(plseq_qp)<1e-8, 0.0, 1.0)
	filename = 'out_'
	for j, i in enumerate(stepList):
		devstress_qp = computeDeviatoric4All(stress_qp[:, :, i])
		vonMises_qp = np.sqrt(3/2)*computeSymTensorNorm4All(devstress_qp)
		problem.part.postProcessingDual(name=filename+str(j), folder=FOLDER2SAVE, 
										fields={'stress': vonMises_qp, 		
												'straineq': plseq_qp[0, :, i], 
												'plastic': plastic_qp[0, :, i]})
	run(folder=FOLDER2SAVE, filename=filename, nbFiles=len(stepList))

	with open(FOLDER2SAVE + 'refpartpl2d.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)
	disp_ref = np.load(FOLDER2SAVE + 'disppl2d.npy')

	for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
		quadArgs = {'quadrule':quadrule, 'type': quadtype}
		errorL2_list = np.ones((len(stepList), len(degList), len(cutList)))
		errorH1_list = np.ones((len(stepList), len(degList), len(cutList)))

		for i, degree in enumerate(degList):
			for j, cuts in enumerate(cutList):
				problem, displacement, _= simulate(degree, cuts, quadArgs)

				for k, step in enumerate(stepList):
					errorL2_list[k, i, j], _ = problem.normOfError(displacement[:, :, step], 
															normArgs={'type':'L2', 
															'part_ref':part_ref, 
															'u_ref': disp_ref[:, :, step]})
					
					errorH1_list[k, i, j], _ = problem.normOfError(displacement[:, :, step], 
															normArgs={'type':'H1', 
															'part_ref':part_ref, 
															'u_ref': disp_ref[:, :, step]})

		quadrule = quadArgs['quadrule']; quadtype = quadArgs['type']
		np.save(FOLDER2SAVE + 'Abserror_pls2d_L2_'+quadrule+str(quadtype), errorL2_list)
		np.save(FOLDER2SAVE + 'Abserror_pls2d_H1_'+quadrule+str(quadtype), errorH1_list)
	
else:

	nbelList = 2**cutList

	for k, step in enumerate(stepList):

		fig, axs = plt.subplots(ncols=2, figsize=(9, 6))

		for error_name, ax in zip(['H1', 'L2'], axs):

			for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 1, 2], [normalPlot, onlyMarker1, onlyMarker2]):

				error_list = np.load(FOLDER2SAVE + 'Abserror_pls2d_' + error_name + '_' + quadrule + str(quadtype) + '.npy')

				for i, degree in enumerate(degList):
					color = COLORLIST[i]
					if quadrule == 'iga': 
						ax.loglog(nbelList, error_list[k, i, :], label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
							markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
						slope = round(np.polyfit(np.log(nbelList[2:]), np.log(error_list[k, i, 2:]), 1)[0], 1)
						annotation.slope_marker((nbelList[-2],  error_list[k, i, -2]), slope, 
										poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
					else: 
						ax.loglog(nbelList[:], error_list[k, i, :], color=color, marker=plotpars['marker'], markerfacecolor='w',
						markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])

			if error_name == 'H1':
				ax.set_ylabel(r'$H^1$' + ' error')
				ax.set_ylim(bottom=1e-9, top=1e-1)
			if error_name == 'L2':
				ax.set_ylabel(r'$L^2$' + ' error')
				ax.set_ylim(bottom=1e-10, top=1e-2)
		
			ax.set_xlabel('Number of elements')
			ax.set_xlim(left=1, right=10**2)

		ax.semilogy([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label='IGA-WQ 1')
		ax.semilogy([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
				markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label='IGA-WQ 2')

		ax.legend()
		fig.tight_layout()
		fig.savefig(FOLDER2SAVE + 'ConvergencePls2d_' + str(k) +'.pdf')
		plt.close(fig)