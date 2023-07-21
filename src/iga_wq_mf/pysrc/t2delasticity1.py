"""
.. Test of elastoplasticity 2D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa (200e3)
..      - Length : mm
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.lib_geomdl import Geomdl
from lib.lib_part import part
from lib.lib_material import mechamat
from lib.lib_boundary import boundaryCondition
from lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t2delastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceVolFun1(P:list):
	ref  = np.array([0.0, 0.0])
	prop = np.zeros((2, np.size(P, axis=1)))
	for i in range(2): prop[i, :] = ref[i] 
	return prop

def forceVolFun2(P:list):
	ref  = np.array([0.0, 0.0])
	prop = np.zeros((2, np.size(P, axis=1)))
	for i in range(2): prop[i, :] = ref[i] 
	return prop

def forceSurfFun(P:list):
	ref  = np.array([4e1, 0.0])
	prop = np.zeros((2, np.size(P, axis=1)))
	for i in range(2): prop[i, :] = ref[i] 
	return prop

# Set global variables
name = 'SQ'
ref = np.load(folder + 'disp_interp_refElasticity2.npy')

degree_list = np.arange(3, 7)
cuts_list   = np.arange(3, 9)
error_disp  = np.ones((len(cuts_list), 3))

fig_disp, axs_disp  = plt.subplots(nrows=1, ncols=3, figsize=(14, 4))
for degree in degree_list:
	for j, cuts in enumerate(cuts_list):

		# Create model 
		geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
					'extra':{'Rin':5.e2, 'Rex':1.e3, 
					'XY':np.array([[0.0, 0.0], [1.e3, 0.0], [1.e3, 1.e3], [0.0, 1.e3]])}}
		quadArgs  = {'quadrule': 'wq', 'type': 1}
		# quadArgs  = {'quadrule': 'iga', 'type': 'leg'}

		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		model    = part(modelIGA, quadArgs=quadArgs)

		# Add material 
		matArgs  = {'elastic_modulus':2e5, 'elastic_limit':100, 'poisson_ratio': 0.3}
		material = mechamat(matArgs)

		# Set Dirichlet boundaries
		boundary = boundaryCondition(model.nbctrlpts)
		table = np.zeros((2, 2, 2), dtype=int)
		if name == 'SQ':
			table[0, 0, 0] = 1
			table[0, 0, 1] = 1
		elif name == 'QA':
			# table[1, 0, 0] = 1
			table[1, 1, 0] = 1
			table[1, 0, 1] = 1
		else: raise Warning('Not possible')
		boundary.add_DirichletDisplacement(table=table)

		# Elasticity problem
		problem = mechaproblem(material, model, boundary)
		Fext = problem.eval_surfForce(forceSurfFun, nbFacePosition=1)
		if name == 'SQ': 
			Fext += problem.eval_volForce(forceVolFun1)
		elif name == 'QA': 
			Fext += problem.eval_volForce(forceVolFun2)

		# -------------
		# ELASTICITY
		# -------------
		displacement, resPCG = problem.solveElasticityProblemFT(Fext=Fext)
		interp = problem.part.interpolateField(u_ctrlpts=displacement, nbDOF=2, sampleSize=2500)[-1]
		disp_norm = np.sqrt(interp[0, :]**2+interp[1, :]**2)
		interp = np.vstack([interp, disp_norm])
		
		for i in range(3):
			norm_ref = np.linalg.norm(ref[i])
			tmp		 = ref[i] - interp[i]
			error_disp[j, i] = np.linalg.norm(tmp)/norm_ref*100

	for [i, ax], title in zip(enumerate(axs_disp), ['Disp x', 'Disp y', 'Disp Norm']):
		ax.semilogy(2**cuts_list, error_disp[:, i], marker=markerSet[i], label='degree p='+str(degree))
		ax.set_title(title)

	for ax in axs_disp:
		ax.set_ylabel('L2 Relative error (\%)')
		ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
		ax.set_ylim(bottom=1e-5, top=1e1)
		ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
		ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		ax.set_xticks([8, 64, 128, 256])

	axs_disp[-1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig_disp.tight_layout()
	fig_disp.savefig(folder + 'Fig1Elasticity2DispErrorNBel_' + str(quadArgs['quadrule']) +'.png')
