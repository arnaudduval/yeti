from lib.__init__ import *
from lib.lib_base import sigmoid
from lib.lib_load import powden
from lib.lib_geomdl import Geomdl
from lib.lib_model import part
from lib.lib_material import thermomat
from lib.lib_step import step
from lib.lib_job import heatproblem

def setKprop(P:list):
	cst = 10.0
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	T = P[3, :]
	# for i in range(3): Kprop[i, i, :] = cst*(1.0 + 2.0/(1.0 + np.exp(-5*(T-1.0))))
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((3, 3, len(x)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j] 
	Kprop[0, 0, :] += 0.75*np.cos(np.pi*y)
	Kprop[1, 1, :] += 2*np.exp(-(z-0.5)**2)
	Kprop[2, 2, :] += 2.5*np.cos(np.pi*x)**2
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j]*cst*(1.0 + 2.0/(1.0 + np.exp(-2.0*(T-1.0))))
	return Kprop 

def setCprop(P:list):
	cst = 1.0
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	T = P[3, :]
	Cprop = cst*(1 + np.exp(-2.0*abs(T)))*np.exp(-(x-0.5)**2)
	return Cprop

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test7/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
dataExist   = True
name_list   = ['CB']
IterMethods = ['WP', 'C', 'JMC']
example     = 1
if   example == 1: nbSteps = 41
elif example == 2: nbSteps = 6
# elif example == 3: nbSteps = 101

if not dataExist:

	degree, cuts = 5, 4
	time_list    = np.linspace(0, 0.25, nbSteps)  

	for name in name_list:
		for PCGmethod in IterMethods:
			filename = folder + 'ResPCG_' + name + '_' + PCGmethod + str(example) + '.dat'        

			# Create model 
			inputs = {'name': name, 'degree':degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
			modelGeo = Geomdl(**inputs)
			modelIGA = modelGeo.getIGAParametrization()
			model    = part(modelIGA)

			# Add material 
			material = thermomat()
			material.addConductivity(setKprop, isIsotropic=False) 
			material.addCapacity(setCprop, isIsotropic=False) 

			# Block boundaries
			boundary = step(model._nbctrlpts)
			boundary.add_DirichletTemperature(table=np.array([[1, 0], [0, 0], [0, 0]]))
			boundary.add_DirichletTemperature(table=np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)

			# ---------------------
			# Transient model
			# ---------------------
			problem = heatproblem(material, model, boundary)
			problem._methodPCG = PCGmethod

			# Create a Dirichlet condition
			Tinout = np.zeros((model._nbctrlpts_total, len(time_list)))
			for i in range(len(time_list)): Tinout[boundary._thdod, i] = boundary._thDirichletBound[boundary._thdod]

			# Add external force 
			Fend = problem.eval_bodyForce(powden)
			Fext = np.kron(np.atleast_2d(Fend).reshape(-1, 1), sigmoid(time_list))

			# Solve
			lastStep = nbSteps
			resPCG = problem.solveNLTransientHeatProblemPy(Tinout=Tinout[:, :lastStep], Fext=Fext[:, :lastStep], 
														time_list=time_list[:lastStep], theta=1.0)
			np.savetxt(filename, resPCG)
			# model.exportResults(u_ctrlpts=Tinout[:, lastStep-1], folder=folder, nbDOF=1)

else:

	# for name in name_list:
	# 	fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

	# 	for i, PCGmethod in enumerate(['C', 'JMC']):
	# 		filename = folder + 'ResPCG_' + name + '_' + PCGmethod + str(example) + '.dat'
	# 		resPCG   = np.loadtxt(filename)

	# 		if   PCGmethod == "C"  : labelmethod = 'Classic FD\nmethod'
	# 		elif PCGmethod == "JMC": labelmethod = 'This work'

	# 		ind = np.where(resPCG[:, 0]==1)
	# 		newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 		ax1.semilogy(np.arange(len(newresidue)), newresidue, '-', color='k', linewidth=0.5)
			
	# 		newresidue = resPCG[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 		ax2.semilogy(np.arange(len(newresidue)), newresidue, '-', color='k', linewidth=0.5)

	# 		maxStep = int(np.max(resPCG[:, 0]))
	# 		for j in range(2, maxStep):
	# 			opacity = (maxStep - j + 1)*1.0/(maxStep - 1)
	# 			ind = np.where(resPCG[:, 0]==j)

	# 			newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 			ax1.semilogy(np.arange(len(newresidue)), newresidue, '-', alpha=opacity, 
	# 						color=colorSet[i], linewidth=0.5)
				
	# 			newresidue = resPCG[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 			ax2.semilogy(np.arange(len(newresidue)), newresidue, '-', alpha=opacity, 
	# 						color=colorSet[i], linewidth=0.5)

	# 	ax1.set_title('First NR iterations')
	# 	ax2.set_title('Last NR iterations')
	# 	for ax in [ax1, ax2]:
	# 		ax.set_xlabel('Number of iterations of BiCGSTAB solver')
	# 		ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
	# 		ax.set_ybound(lower=1e-12, upper=10)

	# 	filename = folder + 'TransientNL_' + name + str(example) + '.png'
	# 	fig.tight_layout()
	# 	fig.savefig(filename)


	for name in name_list:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4.7))

		for i, PCGmethod in enumerate(IterMethods):
			filename = folder + 'ResPCG_' + name + '_' + PCGmethod + str(example) + '.dat'
			resPCG   = np.loadtxt(filename)

			if PCGmethod   == "WP" : labelmethod = 'w.o.\npreconditioner'
			elif PCGmethod == "C"  : labelmethod = 'Classic FD\nmethod'
			elif PCGmethod == "JMC": labelmethod = 'This work'

			# Print the first 
			newresidue = resPCG[0, 2:]; newresidue = newresidue[newresidue>0]
			ax.semilogy(np.arange(len(newresidue)), newresidue, '-', 
						linewidth=2.5, marker=markerSet[i], label=labelmethod)

			# # Print the last
			# newresidue = resPCG[-1, 2:]; newresidue = newresidue[newresidue>0]
			# ax.semilogy(np.arange(len(newresidue)), newresidue, '-', 
			# 			linewidth=2.5, marker=markerSet[i], label=labelmethod)

		ax.legend(bbox_to_anchor=(-0.25, 1.02, 1.25, 0.2), loc='lower left', mode='expand', ncol=3)
		ax.set_xlabel('Number of iterations of BiCGSTAB solver')
		ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
		ax.set_ybound(lower=1e-12, upper=10)

		filename = folder + 'TransientNL_' + name + '.png'
		# fig.tight_layout()
		fig.savefig(filename)