from lib.__init__ import *
from lib.lib_base import (createKnotVector)
from lib.thermomecha1D import mechamat1D, plot_results
from lib.lib_load import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1plasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
matArgs   = {'elastic_modulus':200e3, 'elastic_limit':506, 'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}
length    = 1.0
nbSteps   = 101
geoArgs   = {'length': length}
sampleSize =  2500

degree_list = np.arange(2, 10)
cuts_list   = np.arange(3, 9)

ref = []
ref.append(np.load(folder + 'disp_interp_ref.npy'))
ref.append(np.load(folder + 'stress_interp_ref.npy'))
error = np.zeros((len(cuts_list), 2))

# First curve
fig, axs  = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
for degree in degree_list:
	for j, cuts in enumerate(cuts_list):
		nbel = 2**cuts
		knotvector = createKnotVector(degree, nbel)
		quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'wq', 'type': 1}
		# quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
		args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
		model = mechamat1D(args)

		# Add material
		model.activate_mechanical(matArgs)

		# Add boundary condition
		model.add_DirichletCondition(table=[1, 0])
		nbqpiga = model.nbqp

		# Define boundaries conditions
		Fext        = np.zeros((model.nbctrlpts, nbSteps))
		Fext[:, -1] = model.compute_volForce(forceVol(model.qpPhy))
		for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

		blockPrint()
		# Solve 
		lastStep = -1
		disp_cp, _, stress, _, _ = model.solve(Fext=Fext[:, :lastStep])
		stress_cp = model.interpolate_CntrlPtsField(stress)

		interp = []
		interp.append(model.interpolate_sampleField(disp_cp, sampleSize=sampleSize)[0])
		interp.append(model.interpolate_sampleField(stress_cp, sampleSize=sampleSize)[0])
		enablePrint()

		for i in range(2):
			norm_ref = np.linalg.norm(ref[i][:, :lastStep], axis=0)
			tmp		 = ref[i][:, :lastStep] - interp[i]
			relerror = np.linalg.norm(tmp, axis=0)
			relerror = np.divide(relerror, norm_ref, out=np.zeros_like(relerror), where=np.abs(norm_ref)>1.e-12)*100
			error[j, i] = relerror[-1] # Last step 			

	for [i, ax], title in zip(enumerate(axs), ['Displacement', 'Stress']):
		ax.semilogy(2**cuts_list, error[:, i], label='degree p='+str(degree))
		ax.set_title(title)

for ax in axs:
	ax.set_ylim(bottom=1e-16, top=1e1)
	ax.set_ylabel('L2 Relative error (\%)')
	ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
	ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
	ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
	ax.set_xticks([8, 32, 128, 256])

axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig.tight_layout()
fig.savefig(folder + 'Fig1ErrorNBel_' + str(quadArgs['quadrule']) +'.png')
