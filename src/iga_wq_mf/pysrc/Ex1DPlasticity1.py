from lib.__init__ import *
from lib.lib_base import createKnotVector
from lib.thermomecha1D import mechamat1D, plot_results
from lib.lib_load import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1plasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def run_simulation(degree, knotvector, matArgs, nbSteps, quadrule='iga'):

	# Create geometry
	if   quadrule == 'iga': quadType = 'leg'
	elif quadrule == 'wq' : quadType = 1
	else: raise Warning('Not possible')
	quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': quadrule, 'type': quadType}
	geoArgs   = {'length': 1.e3}
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

	# Solve
	disp_cp, strain_qp, stress_qp, plastic_qp, Cep_qp = model.solve(Fext=Fext)
	print(np.min(stress_qp[:, -1]), np.max(stress_qp[:, -1]))
	strain_cp   = model.interpolate_CntrlPtsField(strain_qp)
	plastic_cp  = model.interpolate_CntrlPtsField(plastic_qp)
	stress_cp 	= model.interpolate_CntrlPtsField(stress_qp)
	plot_results(model.quadRule, geoArgs['length'], disp_cp, plastic_cp, stress_cp, folder=folder, method=quadrule)
	return model, disp_cp, strain_cp, plastic_cp, stress_cp

isReference = True
quadrule = 'wq'

# Set global variables
samplesize = 2500
nbSteps    = 50
matArgs    = {'elastic_modulus':2e5, 'elastic_limit':100,
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}

if isReference:
	degree, nbel = 9, 1024
	knotvector   = createKnotVector(degree, nbel, multiplicity=1)
	model, disp_cp, strain_cp, plastic_cp, stress_cp = run_simulation(degree, knotvector, matArgs, nbSteps, quadrule='iga')
	disp_interp   = model.interpolate_sampleField(disp_cp,   sampleSize=samplesize)[0]
	stress_interp = model.interpolate_sampleField(stress_cp, sampleSize=samplesize)[0]
	np.save(folder+'disp_interp_refcont.npy', disp_interp)
	np.save(folder+'stress_interp_refcont.npy', stress_interp)

else:

	def relativeError(ref, interp):
		norm_ref    = np.linalg.norm(ref, axis=0)
		error       = ref - interp
		norm_error  = np.linalg.norm(error, axis=0)
		relerror    = np.divide(norm_error, norm_ref, out=np.zeros_like(norm_error), where=np.abs(norm_ref)>1.e-12)*100
		return relerror

	ref = []
	ref.append(np.load(folder + 'disp_interp_ref.npy'))
	ref.append(np.load(folder + 'stress_interp_ref.npy'))

	degree = 6
	nbel_list = range(21, 501, 80)
	relerror = np.zeros((len(nbel_list), nbSteps, 3))
	
	for i, nbel in enumerate(nbel_list):
		knotvector = createKnotVector(degree, nbel)
		info = run_simulation(degree, knotvector, matArgs, nbSteps, quadrule=quadrule)
		model, disp_cp, strain_cp, plastic_cp, stress_cp = info
		interp = []
		interp.append(model.interpolate_sampleField(disp_cp, sampleSize=samplesize)[0])
		interp.append(model.interpolate_sampleField(stress_cp, sampleSize=samplesize)[0])

		for j in range(3): relerror[i, :, j] = relativeError(ref[j], interp[j])


	from mpl_toolkits.axes_grid1 import make_axes_locatable
	fig, axs  = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
	nbel_new, steps_new = np.meshgrid(nbel_list, np.arange(0, nbSteps))

	for [i, ax], title in zip(enumerate(axs), ['Displacement', 'Stress']) :
		im = ax.pcolormesh(nbel_new, steps_new, relerror[:, :, i].T,
						norm=mpl.colors.LogNorm(vmin=1.e-11, vmax=1.e-1),
						cmap='viridis', shading='auto')
		ax.set_ylabel('Steps')
		ax.set_xlabel('Number of elements')
		ax.set_title(title + ' error')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes('right', size='5%', pad=0.05)
		cbar = fig.colorbar(im, cax=cax, format='%.1e')
		cbar.ax.set_title('%%')

	fig.tight_layout()
	fig.savefig(folder + 'convergence.png')