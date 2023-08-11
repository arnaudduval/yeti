from lib.__init__ import *
from lib.lib_base import createUniformMaxregularKnotvector
from lib.thermomecha1D import mechamat1D

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceVol(P:list):
	force = 0.4*np.sin(P/1e3)
	return force

# Set global variables
matArgs   = {'elastic_modulus':2e5, 'elastic_limit':100, 'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}
nbSteps   = 50
geoArgs   = {'length': 1.e3}
sampleSize = 2500
FigPlot   = 2

ref = []
ref.append(np.load(folder + 'disp_interp_ref.npy'))
ref.append(np.load(folder + 'stress_interp_ref.npy'))

if FigPlot == 0:
	degree_list = np.arange(2, 6)
	cuts_list   = np.arange(3, 10)
	error = np.zeros((len(cuts_list), 2))

	# First curve
	fig, axs  = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
	for degree in degree_list:
		for j, cuts in enumerate(cuts_list):
			nbel = 2**cuts
			knotvector = createUniformMaxregularKnotvector(degree, nbel)
			quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'wq', 'type': 1}
			# quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
			args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
			model = mechamat1D(args)

			model.activate_mechanical(matArgs)
			model.add_DirichletCondition(table=[1, 0])
			Fext    = np.zeros((model.nbctrlpts, 2*nbSteps + 1))
			Fextref = model.compute_volForce(forceVol(model.qpPhy))
			for i in range(0, nbSteps+1): Fext[:, i] = i/nbSteps*Fextref
			for i in range(nbSteps+1, 2*nbSteps+1): Fext[:, i] = (2*nbSteps - i)/nbSteps*Fextref
			
			blockPrint()
			# Solve 
			lastStep = nbSteps+1 # -2 or nbSteps+1
			disp_cp, _, stress, _, _ = model.solve(Fext=Fext[:, :lastStep])
			stress_cp = model.L2projectionCtrlpts(stress)

			interp = []
			interp.append(model.interpolateMeshgridField(disp_cp, sampleSize=sampleSize)[0])
			interp.append(model.interpolateMeshgridField(stress_cp, sampleSize=sampleSize)[0])
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
		ax.set_ylabel('L2 Relative error (\%)')
		ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
		ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
		ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		ax.set_xticks([8, 32, 128, 256, 512])

	axs[0].set_ylim(bottom=1e-9, top=1e1)
	axs[1].set_ylim(bottom=1e-6, top=1e1)

	axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.tight_layout()
	fig.savefig(folder + 'Fig1ErrorNBel_' + str(quadArgs['quadrule']) +'.png')


elif FigPlot == 1:

	degree = 5
	cuts_list = np.arange(3, 10)
	TYPEOFFIG = 2

	if TYPEOFFIG == 0: TypeList = range(1, 5)
	if TYPEOFFIG == 1: TypeList = range(1, 5)
	if TYPEOFFIG == 2: TypeList = range(2, 6)

	error = np.zeros((len(cuts_list), len(TypeList), 2))
	fig, axs  = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
	for j, typeQuad in enumerate(TypeList):
		for k, cuts in enumerate(cuts_list):
			nbel = 2**cuts
			knotvector = createUniformMaxregularKnotvector(degree, nbel)
			if TYPEOFFIG == 0: quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 
											'type': 'legextra', 'extra':{'nbQPEL':typeQuad}}
			if TYPEOFFIG == 1: quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'wq', 
											'type': 1, 'extra':{'s':typeQuad}}
			if TYPEOFFIG == 2: quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'wq', 
											'type': 1, 'extra':{'r':typeQuad}}
			args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
			model = mechamat1D(args)

			model.activate_mechanical(matArgs)
			model.add_DirichletCondition(table=[1, 0])
			Fext    = np.zeros((model.nbctrlpts, 2*nbSteps + 1))
			Fextref = model.compute_volForce(forceVol(model.qpPhy))
			for i in range(0, nbSteps+1): Fext[:, i] = i/nbSteps*Fextref
			for i in range(nbSteps+1, 2*nbSteps+1): Fext[:, i] = (2*nbSteps - i)/nbSteps*Fextref

			blockPrint()
			# Solve 
			lastStep = nbSteps+1 # -2 or nbSteps+1
			disp_cp, _, stress, _, _ = model.solve(Fext=Fext[:, :lastStep])
			stress_cp = model.L2projectionCtrlpts(stress)

			interp = []
			interp.append(model.interpolateMeshgridField(disp_cp, sampleSize=sampleSize)[0])
			interp.append(model.interpolateMeshgridField(stress_cp, sampleSize=sampleSize)[0])
			enablePrint()

			for l in range(2):
				norm_ref = np.linalg.norm(ref[l][:, :lastStep], axis=0)
				tmp		 = ref[l][:, :lastStep] - interp[l]
				relerror = np.linalg.norm(tmp, axis=0)
				relerror = np.divide(relerror, norm_ref, out=np.zeros_like(relerror), where=np.abs(norm_ref)>1.e-12)*100
				error[k, j, l] = relerror[-1] # Last step 			

		for [l, ax], title in zip(enumerate(axs), ['Displacement', 'Stress']):
			if TYPEOFFIG == 0: ax.semilogy(2**cuts_list, error[:, j, l], color=colorSet[j], marker=markerSet[j], label=r'$p+$' + str(typeQuad) + ' quadPts/el')
			if TYPEOFFIG == 1: ax.semilogy(2**cuts_list, error[:, j, l], color=colorSet[j], marker=markerSet[j], label=str(2+typeQuad) + ' int. points')
			if TYPEOFFIG == 2: ax.semilogy(2**cuts_list, error[:, j, l], color=colorSet[j], marker=markerSet[j], label=r'$p+$' + str(typeQuad) + ' ext. points')
			ax.set_title(title)

	for ax in axs:
		ax.set_ylim(bottom=1e-9, top=1e1)
		ax.set_ylabel('L2 Relative error (\%)')
		ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
		ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
		ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		ax.set_xticks([8, 32, 128, 256, 512])

	axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.tight_layout()
	fig.savefig(folder + 'Fig3ErrorNBel_' + str(quadArgs['quadrule']) + str(TYPEOFFIG) +'.png')

elif FigPlot == 2: 

	degree = 5
	mult_list = [1, 2, 3, 4, 5]
	cuts_list = np.arange(3, 10)
	error = np.zeros((len(cuts_list), 2))

	# First curve
	fig, axs  = plt.subplots(nrows=1, ncols=2, figsize=(11, 4))
	for mult in mult_list:
		for j, cuts in enumerate(cuts_list):
			nbel = 2**cuts
			knotvector = createUniformMaxregularKnotvector(degree, nbel, multiplicity=mult)
			quadArgs   = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
			args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
			model = mechamat1D(args)

			# Add material
			model.activate_mechanical(matArgs)

			# Add boundary condition
			model.add_DirichletCondition(table=[1, 0])

			# Define boundaries conditions
			Fext    = np.zeros((model.nbctrlpts, 2*nbSteps + 1))
			Fextref = model.compute_volForce(forceVol(model.qpPhy))
			for i in range(0, nbSteps+1): Fext[:, i] = i/nbSteps*Fextref
			for i in range(nbSteps+1, 2*nbSteps+1): Fext[:, i] = (2*nbSteps - i)/nbSteps*Fextref

			blockPrint()
			# Solve 
			lastStep = -3 # -2 or nbSteps+1
			disp_cp, _, stress, _, _ = model.solve(Fext=Fext[:, :lastStep])
			stress_cp = model.L2projectionCtrlpts(stress)

			interp = []
			interp.append(model.interpolateMeshgridField(disp_cp, sampleSize=sampleSize)[0])
			interp.append(model.interpolateMeshgridField(stress_cp, sampleSize=sampleSize)[0])
			enablePrint()

			for i in range(2):
				norm_ref = np.linalg.norm(ref[i][:, :lastStep], axis=0)
				tmp		 = ref[i][:, :lastStep] - interp[i]
				relerror = np.linalg.norm(tmp, axis=0)
				relerror = np.divide(relerror, norm_ref, out=np.zeros_like(relerror), where=np.abs(norm_ref)>1.e-12)*100
				error[j, i] = relerror[-1] # Last step 			

		for [i, ax], title in zip(enumerate(axs), ['Displacement', 'Stress']):
			ax.semilogy(2**cuts_list, error[:, i], label='Continuity C'+str(degree-mult))
			ax.set_title(title)

	for ax in axs:
		ax.set_ylabel('L2 Relative error (\%)')
		ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
		ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
		ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		ax.set_xticks([8, 32, 128, 256, 512])

	axs[0].set_ylim(bottom=1e-9, top=1e1)
	axs[1].set_ylim(bottom=1e-6, top=1e1)

	axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.tight_layout()
	fig.savefig(folder + 'Fig3ErrorNBel' + '.png')