from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import sigmoid
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
TRACTION, RINT = 1.0, 1.0
YOUNG, POISSON = 1e3, 0.3
GEONAME = 'QA'
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':2, 'poisson_ratio': POISSON, 
			'plasticLaw': {'Isoname':'linear', 'Eiso':YOUNG/10}}
dataExist  = False

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = RINT**2/r_square # Already squared
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = TRACTION/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = TRACTION/2*3*np.sin(3*theta)*(b**2 - b)
	return F

if not dataExist:

	DEGREE, CUTS = 6, 5
	STEP = 30
	ITERMETHODS = ['C', 'JMC', 'TDC', 'WP']

	for PCGmethod in ITERMETHODS:
		filename = folder + 'ResPCGpls_' + GEONAME + '_' + PCGmethod + '.dat'        

		geoArgs = {'name': 'QA', 'degree': DEGREE*np.ones(3, dtype=int), 
						'nb_refinementByDirection': CUTS*np.ones(3, dtype=int), 
						'extra':{'Rin':1.0, 'Rex':4.0}
				}
		quadArgs  = {'quadrule': 'wq', 'type': 1}
		
		blockPrint()
		material = mechamat(MATARGS)
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)
		meshparam = modelPhy.compute_mesh_parameter()

		# Set Dirichlet boundaries
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		table = np.zeros((2, 2, 2), dtype=int)
		table[1, 1, 0] = 1
		table[1, 0, 1] = 1
		boundary.add_DirichletDisplacement(table=table)
		enablePrint()

		# Solve elastic problem
		solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-10, 'PCGmethod': PCGmethod, 'NRThreshold': 1e-9}
		problem = mechaproblem(material, modelPhy, boundary)
		problem.addSolverConstraints(solverArgs=solverArgs)
		Fref = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
		Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
		for k in range(len(TIME_LIST)): Fext_list[:, :, k] = np.sin(TIME_LIST[k])*Fref
		displacement, resPCG = problem.solvePlasticityProblemPy(Fext_list=Fext_list[:, :, :STEP + 1])[0]
		np.savetxt(filename, resPCG)

else:

	pass
	# for geo in geo_list:
	# 	fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

	# 	for i, PCGmethod in enumerate(['C', 'JMC']):
	# 		filename = folder + 'ResPCG_' + geo + '_' + PCGmethod + str(example) + '.dat'
	# 		resPCG   = np.loadtxt(filename)

	# 		if   PCGmethod == "C"  : labelmethod = 'Classic FD\nmethod'
	# 		elif PCGmethod == "JMC": labelmethod = 'This work'

	# 		ind = np.where(resPCG[:, 0]==1)
	# 		newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 		axs[0].semilogy(np.arange(len(newresidue)), newresidue, '-', color='k', label=labelmethod, linewidth=0.5)
			
	# 		newresidue = resPCG[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 		axs[1].semilogy(np.arange(len(newresidue)), newresidue, '-', color='k', label=labelmethod, linewidth=0.5)

	# 		maxStep = int(np.max(resPCG[:, 0]))
	# 		for j in range(2, maxStep):
	# 			opacity = (maxStep - j + 1)*1.0/(maxStep - 1)
	# 			ind = np.where(resPCG[:, 0]==j)

	# 			newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 			axs[0].semilogy(np.arange(len(newresidue)), newresidue, '-', alpha=opacity, 
	# 						color=COLORLIST[i], linewidth=0.5)
				
	# 			newresidue = resPCG[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 			axs[1].semilogy(np.arange(len(newresidue)), newresidue, '-', alpha=opacity, 
	# 						color=COLORLIST[i], linewidth=0.5)

	# 	axs[0].set_title('First NR iterations')
	# 	axs[1].set_title('Last NR iterations')
	# 	for ax in axs:
	# 		ax.set_xlabel('Number of iterations of BiCGSTAB solver')
	# 		ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_2}{||b||_2}$')
	# 		ax.set_ybound(lower=1e-12, upper=10)

	# 	axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
	# 	filename = folder + 'TransientNL_' + geo + str(example) + '.pdf'
	# 	fig.tight_layout()
	# 	fig.savefig(filename)

	# for geo in geo_list:
	# 	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4.7))
	# 	# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))

	# 	for i, PCGmethod in enumerate(IterMethods):
	# 		filename = folder + 'ResPCG_' + geo + '_' + PCGmethod + str(example) + '.dat'
	# 		resPCG   = np.loadtxt(filename)

	# 		if PCGmethod   == "WP" : labelmethod = 'w.o.\npreconditioner'
	# 		elif PCGmethod == "C"  : labelmethod = 'Classic FD\nmethod'
	# 		elif PCGmethod == "JMC": labelmethod = 'This work'
			
	# 		ind = np.where(resPCG[:, 0]==1) # or 
	# 		newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 		ax.semilogy(np.arange(len(newresidue)), newresidue, '-', 
	# 					linewidth=2.5, marker=markerSet[i], label=labelmethod)

	# 	ax.legend(bbox_to_anchor=(-0.25, 1.02, 1.25, 0.2), loc='lower left', mode='expand', ncol=3)
	# 	ax.set_xlabel('Number of iterations of BiCGSTAB solver')
	# 	ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_2}{||b||_2}$')
	# 	ax.set_ybound(lower=1e-12, upper=10)

	# 	filename = folder + 'TransientNL_' + geo + '.pdf'
	# 	fig.tight_layout()
	# 	fig.savefig(filename)