"""
.. Test of elastoplasticity 1D
.. Joaquin Cornejo 
"""

import pickle
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve, createAsymmetricalCurve
from pysrc.lib.lib_1d import mechaproblem1D
from pysrc.lib.lib_material import mechamat

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Global variables
YOUNG, CST, LENGTH  = 2e11, 4.e8, 1
NBSTEPS = 251
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e6, 'poisson_ratio':0.3,
		'isoHardLaw': {'Isoname':'linear', 'Eiso':YOUNG/10}}
MECHAMATERIAL = mechamat(MATARGS)
isReference = False

def forceVol(P:list):
	force = CST*(P - 1/10*P**2)
	return force

def simulate(degree, nbel, args, step=-2):
	crv = createUniformCurve(degree, nbel, LENGTH)
	# crv = createAsymmetricalCurve(degree, nbel, LENGTH, xasym=0.75)
	modelPhy = mechaproblem1D(crv, args)
	model2return = deepcopy(modelPhy)
	modelPhy.activate_mechanical(MECHAMATERIAL)
	modelPhy.add_DirichletCondition(table=[1, 1])
	Fref = np.atleast_2d(modelPhy.compute_volForce(forceVol)).transpose()
	Fext_list = np.kron(Fref, np.sin(TIME_LIST))
	displacement = np.zeros(np.shape(Fext_list))
	blockPrint()
	strain, stress, plasticeq, Cep = modelPhy.solvePlasticityProblem(displacement, Fext_list[:, :step+1])
	enablePrint()
	return model2return, displacement[:, :step+1], stress, plasticeq


if isReference:

	degree, nbel = 2, 4096
	args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
	modelPhy, displacement, stress, plasticeq = simulate(degree, nbel, args)
	np.save(folder + 'disppl', displacement)
	with open(folder + 'refpartpl.pkl', 'wb') as outp:
		pickle.dump(modelPhy, outp, pickle.HIGHEST_PROTOCOL)

	# plasticeq_cp = modelPhy.L2projectionCtrlpts(plasticeq)
	# stress_cp  = modelPhy.L2projectionCtrlpts(stress)
	# from mpl_toolkits.axes_grid1 import make_axes_locatable
	# basis = modelPhy.quadRule.getSampleBasis(sampleSize=101)[0]
	# displacement_interp = basis[0].T @ displacement
	# plasticeq_interp = basis[0].T @ plasticeq_cp
	# stress_interp  = basis[0].T @ stress_cp
	# qpPhy_interp   = basis[0].T @ modelPhy.ctrlpts

	# # Plot fields
	# XX, STEPS = np.meshgrid(qpPhy_interp, np.arange(1, NBSTEPS))
	# names = ['Displacement field', 'Plastic strain field', 'Stress field']
	# fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
	# for ax, variable, name in zip([ax1, ax2, ax3], [displacement_interp, plasticeq_interp, stress_interp], names):
	# 	im = ax.pcolormesh(XX, STEPS, variable.T, cmap='PuBu_r', shading='linear')
	# 	ax.set_title(name)
	# 	ax.set_ylabel('Step')
	# 	ax.set_xlabel('Position')
	# 	ax.grid(False)
	# 	divider = make_axes_locatable(ax)
	# 	cax  = divider.append_axes('right', size='5%', pad=0.05)
	# 	cbar = fig.colorbar(im, cax=cax)

	# fig.tight_layout()
	# fig.savefig(folder + 'ElastoPlasticity1D' + '.png')

else: 

	disp_ref = np.load(folder + 'disppl.npy')
	with open(folder + 'refpartpl.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)

	degree_list = np.arange(1, 4)
	cuts_list   = np.arange(2, 9)
	step_max    = 130
	step_list  = range(50, step_max, 5)
	error_list = np.ones((len(step_list), len(degree_list), len(cuts_list)))

	for i, degree in enumerate(degree_list):
		for j, cuts in enumerate(cuts_list):
			nbel = 2**cuts
			args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
			modelPhy, displacement, _, _ = simulate(degree, nbel, args, step=step_max)

			for k, step in enumerate(step_list):
				error_list[k, i, j] = modelPhy.normOfError(displacement[:, step], normArgs={'type':'H1', 'part_ref':part_ref, 
																			'u_ref': disp_ref[:, step]}, isRelative=False)	

	np.save(folder + 'plasticity1D', error_list)
	error_list = np.load(folder + 'plasticity1D.npy')
	

	for k, step in enumerate(step_list):
		fig, ax = plt.subplots(figsize=(9, 6))
		for i, degree in enumerate(degree_list):
			color = COLORLIST[i]
			ax.loglog(2**cuts_list, error_list[k, i, :], color=color, marker='o', markerfacecolor='w',
						markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))

		ax.set_ylabel(r'$H^1$'+ ' error')
		ax.set_xlabel('Number of elements')
		ax.set_ylim(bottom=1e-12, top=1e-3)
		ax.set_xlim(left=1, right=10**3)

		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigPlasticity1_' + str(k) +'.pdf')