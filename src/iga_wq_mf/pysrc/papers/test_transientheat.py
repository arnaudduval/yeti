from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import sigmoid
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import heatproblem

def conductivityProperty(P:list):
	cst = 10.0
	T   = P[3, :]
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((3, 3, len(T)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j] 
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j]*cst*(1.0 + 2.0/(1.0 + np.exp(-5.0*(T-1.0))))
	return Kprop 

def capacityProperty(P:list):
	cst = 1.0
	T   = P[3, :]
	Cprop = cst*(1 + np.exp(-2.0*abs(T)))
	return Cprop

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/heattransfer/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
dataExist   = False
geo_list    = ['cb', 'vb']
IterMethods = ['C', 'JMC', 'TDC', 'WP']
example     = 2
if   example == 1: nbsteps = 41
elif example == 2: nbsteps = 6

if not dataExist:

	degree, cuts = 6, 5
	time_list    = np.linspace(0, 0.25, nbsteps)  

	for PCGmethod in IterMethods:
		for geo in geo_list:
			filename = folder + 'ResPCG_' + geo + '_' + PCGmethod + str(example) + '.dat'        

			# Create model 
			geoArgs = {'name': geo, 'degree': degree*np.ones(3, dtype=int), 
						'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
			quadArgs  = {'quadrule': 'wq', 'type': 1}

			modelGeo = Geomdl(geoArgs)
			modelIGA = modelGeo.getIGAParametrization()
			modelPhy = part(modelIGA, quadArgs=quadArgs)

			# Add material 
			material = heatmat()
			material.addConductivity(conductivityProperty, isIsotropic=False) 
			material.addCapacity(capacityProperty, isIsotropic=False) 

			# Block boundaries
			boundary = boundaryCondition(modelPhy.nbctrlpts)
			boundary.add_DirichletConstTemperature(table=np.array([[1, 0], [0, 0], [0, 0]]))
			boundary.add_DirichletConstTemperature(table=np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)

			# ---------------------
			# Transient model
			# ---------------------			
			SOLVERARGS = {'preconditioner': PCGmethod}
			problem = heatproblem(material, modelPhy, boundary)
			problem.addSolverConstraints(solverArgs=SOLVERARGS)

			# Create a Dirichlet condition
			Tinout = np.zeros((modelPhy.nbctrlpts_total, len(time_list)))
			for i in range(1, len(time_list)): Tinout[boundary.thdod, i] = boundary.thDirichletBound[boundary.thdod]

			# Add external force 
			Fend = np.zeros((problem.part.nbctrlpts_total, 1))
			Fext = np.kron(Fend, sigmoid(time_list))

			# Solve
			lastStep = 3
			resPCG = problem.solveFourierTransientProblem(Tinout=Tinout[:, :lastStep], Fext_list=Fext[:, :lastStep], 
														time_list=time_list[:lastStep], alpha=1.0)
			np.savetxt(filename, resPCG)

else:

	for geo in geo_list:
		fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

		for i, PCGmethod in enumerate(['C', 'JMC']):
			filename = folder + 'ResPCG_' + geo + '_' + PCGmethod + str(example) + '.dat'
			resPCG   = np.loadtxt(filename)

			if   PCGmethod == "C"  : labelmethod = 'Classic FD\nmethod'
			elif PCGmethod == "JMC": labelmethod = 'This work'

			ind = np.where(resPCG[:, 0]==1)
			newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
			axs[0].semilogy(np.arange(len(newresidue)), newresidue, '-', color='k', label=labelmethod, linewidth=0.5)
			
			newresidue = resPCG[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
			axs[1].semilogy(np.arange(len(newresidue)), newresidue, '-', color='k', label=labelmethod, linewidth=0.5)

			maxStep = int(np.max(resPCG[:, 0]))
			for j in range(2, maxStep):
				opacity = (maxStep - j + 1)*1.0/(maxStep - 1)
				ind = np.where(resPCG[:, 0]==j)

				newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
				axs[0].semilogy(np.arange(len(newresidue)), newresidue, '-', alpha=opacity, 
							color=COLORLIST[i], linewidth=0.5)
				
				newresidue = resPCG[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
				axs[1].semilogy(np.arange(len(newresidue)), newresidue, '-', alpha=opacity, 
							color=COLORLIST[i], linewidth=0.5)

		axs[0].set_title('First NR iterations')
		axs[1].set_title('Last NR iterations')
		for ax in axs:
			ax.set_xlabel('Number of iterations of BiCGSTAB solver')
			ax.set_ylabel('Relative residue')
			ax.set_ybound(lower=1e-12, upper=10)

		axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
		filename = folder + 'TransientNL_' + geo + str(example) + '.pdf'
		fig.tight_layout()
		fig.savefig(filename)

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