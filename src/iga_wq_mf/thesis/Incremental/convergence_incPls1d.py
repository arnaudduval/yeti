"""
.. Test of elastoplasticity 1D
.. Joaquin Cornejo 
"""

from thesis.Incremental.__init__ import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pysrc.lib.lib_job1d import problem1D

NBSTEPS = 201
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
YOUNG, CST, LENGTH  = 2e11, 4.e7, 1
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e7, 'poisson_ratio':0.3,
		'isoHardLaw': {'name':'linear', 'Eiso':YOUNG/10}}
MECHAMATERIAL = mechamat(MATARGS)

def forceVol(args:list):
	P = args['position']
	force = CST*(P - 1/10*P**2)
	return force

def simulate_1d(degree, nbel, quadArgs={}):
	modelPhy = part1D(createUniformOpenCurve(degree, nbel, LENGTH), quadArgs)
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	boundary.add_DirichletConstTemperature(table=np.array([[1, 1]]))
	problem = mechaproblem1D(mechanical_material=MECHAMATERIAL, part=modelPhy, boundary=boundary)
	problem2return = problem1D(part=modelPhy, boundary=boundary)
	Fref = np.atleast_2d(problem.compute_volForce(forceVol)).transpose()
	FextList = np.kron(Fref, np.sin(TIME_LIST))
	displacement = np.zeros(np.shape(FextList))
	_, stress, plseq, _ = problem.solvePlasticityProblem(displacement, FextList)
	return problem2return, displacement, stress, plseq

RUNSIMU = False
degList = np.arange(1, 4)
cutList = np.arange(1, 10)
stepMax = np.max([201, NBSTEPS])
stepList = range(1, stepMax, 10)

if RUNSIMU:

	degree, nbel = 2, 2048
	quadArgs = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
	problem, displacement_cp, stress_qp, plseq_qp = simulate_1d(degree, nbel, quadArgs)
	np.save(FOLDER2DATA + 'disppl1d', displacement_cp)
	with open(FOLDER2DATA + 'refpartpl1d.pkl', 'wb') as outp:
		pickle.dump(problem, outp, pickle.HIGHEST_PROTOCOL)

	disp_ref = np.load(FOLDER2DATA + 'disppl1d.npy')
	with open(FOLDER2DATA + 'refpartpl1d.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)

	errorL2_list = np.ones((len(stepList), len(degList), len(cutList)))
	errorH1_list = np.ones((len(stepList), len(degList), len(cutList)))

	for i, degree in enumerate(degList):
		for j, cuts in enumerate(cutList):
			nbel = 2**cuts
			quadArgs = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
			problem, displacement_cp, stress_qp, plseq_qp = simulate_1d(degree, nbel, quadArgs)

			for k, step in enumerate(stepList):
				errorL2_list[k, i, j], _ = problem.normOfError(displacement_cp[:, step], 
																normArgs={'type':'L2', 
																		'part_ref':part_ref, 
																		'u_ref': disp_ref[:, step]})	
				errorH1_list[k, i, j], _ = problem.normOfError(displacement_cp[:, step], 
																normArgs={'type':'H1', 
																		'part_ref':part_ref, 
																		'u_ref': disp_ref[:, step]})	

	np.save(FOLDER2DATA + 'Abserror_pls1d_L2', errorL2_list)
	np.save(FOLDER2DATA + 'Abserror_pls1d_H1', errorH1_list)

	plasticeq_cp = problem.L2projectionCtrlpts(plseq_qp)
	stress_cp = problem.L2projectionCtrlpts(stress_qp)
	basis = problem.part.quadRule.getSampleBasis(sampleSize=101)[0]
	displacement_interp = basis[0].T @ displacement_cp
	plseq_interp = basis[0].T @ plasticeq_cp
	stress_interp = basis[0].T @ stress_cp
	qpPhy_interp = basis[0].T @ problem.part.ctrlpts

	# Plot fields
	XX, STEPS = np.meshgrid(qpPhy_interp, np.arange(np.size(displacement_cp, axis=1)))
	names = ['Displacement field', 'Plastic strain field', 'Stress field']
	fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
	for ax, variable, name in zip([ax1, ax2, ax3], [displacement_interp, plseq_interp, stress_interp], names):
		im = ax.pcolormesh(XX, STEPS, variable.T, cmap='PuBu_r', shading='linear')
		ax.set_title(name)
		ax.set_ylabel('Step')
		ax.set_xlabel('Position')
		ax.grid(False)
		divider = make_axes_locatable(ax)
		cax  = divider.append_axes('right', size='5%', pad=0.05)
		cbar = fig.colorbar(im, cax=cax)
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE + 'ElastoPlasticity1D' + '.png')

else: 
	
	
	for k, step in enumerate(stepList):

		fig, axs = plt.subplots(ncols=2, figsize=(9, 6))

		for error_name, ax in zip(['H1', 'L2'], axs):

			error_list = np.load(FOLDER2DATA + 'Abserror_pls1d_' + error_name + '.npy')

			for i, degree in enumerate(degList):
				color = COLORLIST[i]
				ax.loglog(2**cutList, error_list[k, i, :], color=color, marker='s', markerfacecolor='w',
							markersize=10, linestyle='-', label='IGA-GL deg. ' + str(degree))
				
				if degree < 4 and step < 55:
					slope = round(np.polyfit(np.log(2**cutList[:-3]), np.log(error_list[k, i, :-3]), 1)[0], 1)
					annotation.slope_marker((2**cutList[1],  error_list[k, i, 1]), slope, 
									poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

			if error_name == 'H1':
				ax.set_ylabel(r'$H^1$' + ' error')
				ax.set_ylim(bottom=1e-12, top=1e-4)
			if error_name == 'L2':
				ax.set_ylabel(r'$L^2$' + ' error')
				ax.set_ylim(bottom=1e-13, top=1e-5)
		
			ax.set_xlabel('Number of elements')
			ax.set_xlim(left=1, right=10**3)

		ax.legend()
		fig.tight_layout()
		fig.savefig(FOLDER2SAVE + 'ConvergencePls1d_' + str(k) +'.pdf')
		plt.close(fig)