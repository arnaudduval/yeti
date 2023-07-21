"""
.. Test of elastoplasticity 1D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa
..      - Length : mm
..      - Force  : N
..      - Mass   : metric ton 
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.lib_geomdl import Geomdl
from lib.lib_part import part
from lib.lib_material import mechamat, computeMultiVMStressVgt
from lib.lib_boundary import boundaryCondition
from lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t2delastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def run_simulation(degree, cuts, matArgs, nbSteps, quadrule='iga', extra=None):

	# Create geometry
	geoName = 'QA'
	if   quadrule == 'iga': quadType = 'leg'
	elif quadrule == 'wq' : quadType = 1
	else: raise Warning('Not possible')
	if quadrule == 'iga' and extra is not None: quadType = 'legextra'
	if extra is None: extra = {}
	geoArgs = {'name': geoName, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'Rin':5.e2, 'Rex':1.e3}}
	quadArgs  = {'quadrule': quadrule, 'type': quadType, 'extra':extra}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	model    = part(modelIGA, quadArgs=quadArgs)

	# Add material
	material = mechamat(matArgs)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(model.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[1, 1, 0] = 1
	table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)

	# Elasticity problem
	problem = mechaproblem(material, model, boundary)

	def forceSurfFun(P:list):
		ref  = np.array([4e1, 0.0])
		prop = np.zeros((2, np.size(P, axis=1)))
		for i in range(2): prop[i, :] = ref[i] 
		return prop
	
	Fextref = problem.eval_surfForce(forceSurfFun, nbFacePosition=1)
	Fext = np.zeros((*np.shape(Fextref), nbSteps))
	for i in range(0, nbSteps): Fext[:, :, i] = i/(nbSteps-1)*Fextref

	# Solve
	disp_cp, _, stress_qp = problem.solvePlasticityProblemPy(Fext=Fext)

	return problem, disp_cp, stress_qp


# Set global variables
quadrule = 'wq'
matArgs  = {'elastic_modulus':2e5, 'elastic_limit':100, 'poisson_ratio': 0.3,
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}
nbSteps  = 41
sampleSize = 2500
FigPlot    = 1

ref = []
ref.append(np.load(folder + 'disp_interp_ref2D.npy'))
ref.append(np.load(folder + 'stress_interp_ref2D.npy'))

if FigPlot == 0:
	degree_list = np.arange(3, 7)
	cuts_list   = np.arange(3, 7)
	CPUtime  = np.ones((len(cuts_list), 3))
	error_str   = np.ones((len(cuts_list), 4))

	# First curve
	fig, axs  = plt.subplots(nrows=1, ncols=3, figsize=(14, 4))
	fig_str, axs_str  = plt.subplots(nrows=2, ncols=2, figsize=(8, 6))
	for degree in degree_list:
		for j, cuts in enumerate(cuts_list):
			interp = []
			problem, disp_cp, stress_qp = run_simulation(degree, cuts, matArgs, nbSteps, quadrule=quadrule)
			
			disp_interp = problem.part.interpolateField(u_ctrlpts=disp_cp[:,:,-2], nbDOF=2, sampleSize=sampleSize)[-1]
			disp_norm = np.sqrt(disp_interp[0, :]**2+disp_interp[1, :]**2)
			disp_interp = np.vstack([disp_interp, disp_norm])
			interp.append(disp_interp)

			stress_cp = problem.solveInterpolationProblemFT(datafield=stress_qp[:,:,-2])
			stress_interp = problem.part.interpolateField(u_ctrlpts=stress_cp, nbDOF=3, sampleSize=sampleSize)[-1]
			stress_vm = computeMultiVMStressVgt(stress_interp, dim=2)
			stress_interp = np.vstack([stress_interp, stress_vm])
			interp.append(stress_interp)
			
			for i in range(3):
				norm_ref = np.linalg.norm(ref[0][i])
				tmp		 = ref[0][i] - interp[0][i]
				CPUtime[j, i] = np.linalg.norm(tmp)/norm_ref*100

			for i in range(4):
				norm_ref = np.linalg.norm(ref[1][i])
				tmp		 = ref[1][i] - interp[1][i]
				error_str[j, i] = np.linalg.norm(tmp)/norm_ref*100

		for [i, ax], title in zip(enumerate(axs), ['Disp x', 'Disp y', 'Disp Norm']):
			ax.semilogy(2**cuts_list, CPUtime[:, i], marker=markerSet[i], label='degree p='+str(degree))
			ax.set_title(title)

		for [i, ax], title in zip(enumerate(np.ravel(axs_str)), ['Stress xx', 'Stress yy', 'Stress xy', 'Stress VM']):
			ax.semilogy(2**cuts_list, error_str[:, i], marker=markerSet[i], label='degree p='+str(degree))
			ax.set_title(title)

		for ax in axs:
			ax.set_ylabel('L2 Relative error (\%)')
			ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
			ax.set_ylim(bottom=1e-4, top=1e1)
			ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
			ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
			ax.set_xticks([8, 32, 64])

		for ax in np.ravel(axs_str):
			ax.set_ylabel('L2 Relative error (\%)')
			ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
			ax.set_ylim(bottom=1e-2, top=1e1)
			ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
			ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
			ax.set_xticks([8, 32, 64])

		axs[-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

		fig.tight_layout()
		fig.savefig(folder + 'Fig1DispErrorNBel_' + str(quadrule) +'.png')

		fig_str.tight_layout()
		fig_str.savefig(folder + 'Fig1StrErrorNBel_' + str(quadrule) +'.png')

elif FigPlot == 1:
	degree_list = np.arange(3, 6)
	cuts_list   = np.arange(3, 6)
	CPUtime     = np.ones(len(cuts_list))

	# First curve
	fig, axs  = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
	for i, quadrule in enumerate(['wq', 'iga']):
		for degree in degree_list:
			for j, cuts in enumerate(cuts_list):
				start = time.process_time()
				problem, disp_cp, stress_qp = run_simulation(degree, cuts, matArgs, nbSteps, quadrule=quadrule)
				stop = time.process_time()
				CPUtime[j] = stop - start

			axs[i].semilogy(2**cuts_list, CPUtime, label='degree p='+str(degree))

		for ax in axs:
			ax.set_ylabel('Time (s)')
			ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
			ax.set_ylim(bottom=1e0, top=1e3)
			ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
			ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
			ax.set_xticks([8, 16, 32])

		axs[-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'Fig2Time' +'.png')


elif FigPlot == 2:

	degree = 6
	cuts_list = np.arange(3, 7)
	TYPEOFFIG = 0

	if TYPEOFFIG == 0: TypeList = range(1, 5)
	if TYPEOFFIG == 1: TypeList = range(1, 5)
	if TYPEOFFIG == 2: TypeList = range(2, 6)

	CPUtime  = np.ones((len(cuts_list), 3))
	error_str   = np.ones((len(cuts_list), 4))
	# First curve
	fig, axs  = plt.subplots(nrows=1, ncols=3, figsize=(14, 4))
	fig_str, axs_str  = plt.subplots(nrows=2, ncols=2, figsize=(8, 6))
	for typeQuad in TypeList:
		for j, cuts in enumerate(cuts_list):
			if TYPEOFFIG == 0: 
				quadrule = 'iga'; extra={'nbQPEL':typeQuad}
			if TYPEOFFIG == 1: 
				quadrule = 'wq'; extra={'s':typeQuad}
			if TYPEOFFIG == 2: 
				quadrule = 'wq'; extra={'r':typeQuad}
			interp = []
			problem, disp_cp, stress_qp = run_simulation(degree, cuts, matArgs, nbSteps, quadrule=quadrule, extra=extra)
			
			disp_interp = problem.part.interpolateField(u_ctrlpts=disp_cp[:,:,-2], nbDOF=2, sampleSize=sampleSize)[-1]
			disp_norm = np.sqrt(disp_interp[0, :]**2+disp_interp[1, :]**2)
			disp_interp = np.vstack([disp_interp, disp_norm])
			interp.append(disp_interp)

			stress_cp = problem.solveInterpolationProblemFT(datafield=stress_qp[:,:,-2])
			stress_interp = problem.part.interpolateField(u_ctrlpts=stress_cp, nbDOF=3, sampleSize=sampleSize)[-1]
			stress_vm = computeMultiVMStressVgt(stress_interp, dim=2)
			stress_interp = np.vstack([stress_interp, stress_vm])
			interp.append(stress_interp)

			for i in range(3):
				norm_ref = np.linalg.norm(ref[0][i])
				tmp		 = ref[0][i] - interp[0][i]
				CPUtime[j, i] = np.linalg.norm(tmp)/norm_ref*100

			for i in range(4):
				norm_ref = np.linalg.norm(ref[1][i])
				tmp		 = ref[1][i] - interp[1][i]
				error_str[j, i] = np.linalg.norm(tmp)/norm_ref*100

		for [i, ax], title in zip(enumerate(axs), ['Disp x', 'Disp y', 'Disp Norm']):
			if TYPEOFFIG == 0: ax.semilogy(2**cuts_list, CPUtime[:, i], marker=markerSet[i], label=r'$p+$' + str(typeQuad) + ' quadPts/el')
			if TYPEOFFIG == 1: ax.semilogy(2**cuts_list, CPUtime[:, i], marker=markerSet[i], label=str(2+typeQuad) + ' int. points')
			if TYPEOFFIG == 2: ax.semilogy(2**cuts_list, CPUtime[:, i], marker=markerSet[i], label=r'$p+$' + str(typeQuad) + ' ext. points')
			ax.set_title(title)

		for [i, ax], title in zip(enumerate(np.ravel(axs_str)), ['Stress xx', 'Stress yy', 'Stress xy', 'Stress VM']):
			if TYPEOFFIG == 0: ax.semilogy(2**cuts_list, error_str[:, i], marker=markerSet[i], label=r'$p+$' + str(typeQuad) + ' quadPts/el')
			if TYPEOFFIG == 1: ax.semilogy(2**cuts_list, error_str[:, i], marker=markerSet[i], label=str(2+typeQuad) + ' int. points')
			if TYPEOFFIG == 2: ax.semilogy(2**cuts_list, error_str[:, i], marker=markerSet[i], label=r'$p+$' + str(typeQuad) + ' ext. points')
			ax.set_title(title)

		for ax in axs:
			ax.set_ylabel('L2 Relative error (\%)')
			ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
			ax.set_ylim(bottom=1e-4, top=1e1)
			ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
			ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
			ax.set_xticks([8, 32, 64])

		for ax in np.ravel(axs_str):
			ax.set_ylabel('L2 Relative error (\%)')
			ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
			ax.set_ylim(bottom=1e-2, top=1e1)
			ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
			ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
			ax.set_xticks([8, 32, 64])

		axs[-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

		fig.tight_layout()
		fig.savefig(folder + 'Fig3DispErrorNBel_' + str(quadrule) + str(TYPEOFFIG) +'.png')

		fig_str.tight_layout()
		fig_str.savefig(folder + 'Fig3StrErrorNBel_' + str(quadrule) + str(TYPEOFFIG) +'.png')

