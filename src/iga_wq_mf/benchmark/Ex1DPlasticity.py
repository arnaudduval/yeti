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

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformMaxregularKnotvector
from pysrc.lib.thermomecha1D import mechamat1D, plot_results
from pysrc.lib.lib_load import *

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
	plastic_cp  = model.L2projectionCtrlpts(plastic_qp)
	stress_cp 	= model.L2projectionCtrlpts(stress_qp)
	plot_results(model.quadRule, geoArgs['length'], disp_cp, plastic_cp, stress_cp, folder=folder, method=quadrule)
	return model, disp_cp, plastic_cp, stress_cp

isReference = False
quadrule = 'wq'

# Set global variables
samplesize = 2500
nbSteps    = 50
matArgs    = {'elastic_modulus':2e5, 'elastic_limit':100,
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}

if isReference:
	degree, nbel = 4, 32
	knotvector   = createUniformMaxregularKnotvector(degree, nbel, multiplicity=degree)
	model, disp_cp, plastic_cp, stress_cp = run_simulation(degree, knotvector, matArgs, nbSteps, quadrule='iga')
	disp_interp   = model.interpolateMeshgridField(disp_cp,   sampleSize=samplesize)[0]
	stress_interp = model.interpolateMeshgridField(stress_cp, sampleSize=samplesize)[0]

	# basis, knots  = model.quadRule.getSampleBasis(sampleSize=samplesize)
	# strain_interp = basis[1].T @ disp_cp
	# np.save(folder+'disp_interp_ref.npy', disp_interp)
	# np.save(folder+'strain_interp_ref.npy', strain_interp)
	# np.save(folder+'stress_interp_ref.npy', stress_interp)

else:

	ref = []
	ref.append(np.load(folder + 'disp_interp_ref.npy'))
	ref.append(np.load(folder + 'strain_interp_ref.npy'))
	ref.append(np.load(folder + 'stress_interp_ref.npy'))

	degree, nbel = 9, 1024
	knotvector   = createUniformMaxregularKnotvector(degree, nbel, multiplicity=1)
	quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
	geoArgs   = {'length': 1.e3}
	args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
	model = mechamat1D(args)
	qpPhy = np.linspace(0, 1, samplesize)*geoArgs['length']
	plastic = ref[1][:, nbSteps] - ref[2][:, nbSteps]/matArgs['elastic_modulus']

	fig, axs  = plt.subplots(nrows=1, ncols=3, figsize=(15, 4))
	axs[0].plot(qpPhy, ref[1][:, nbSteps]*1e-6)
	axs[1].plot(qpPhy, plastic*1e-6)
	axs[2].plot(qpPhy, ref[2][:, nbSteps])

	for ax in axs:
		ax.set_xlabel('Position (mm)')
		ax.set_ylim(bottom=0)

	for i in range(2):
		axs[i].set_ylim(top=5e-6)
	axs[-1].set_ylim(top=200)
	axs[0].set_ylabel('Strain (m/m)')
	axs[1].set_ylabel('Plastic strain (m/m)')
	axs[2].set_ylabel('Stress (MPa)')
	fig.tight_layout()
	fig.savefig(folder + 'DataMaxLoad' + '.png')