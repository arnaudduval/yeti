from lib.__init__ import *
from lib.lib_base import sigmoid
from lib.lib_load import powden
from lib.lib_geomdl import Geomdl
from lib.lib_model import part
from lib.lib_material import thermomat
from lib.lib_step import step
from lib.lib_job import heatproblem

def setKprop(P:list):
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	T = P[3, :]

	return 

def setCprop(P:list):
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	T = P[3, :]

	return 


# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test7/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
dataExist   = False
name_list   = ['VB']
IterMethods = ['WP', 'C', 'JMC']

if not dataExist:

	degree, cuts = 6, 5
	time_list    = np.linspace(0, 20, 41)  

	for name in name_list:
		for PCGmethod in IterMethods:
			filename = folder + 'ResPCG_' + name + '_' + PCGmethod + '.dat'        

			# Create model 
			inputs = {'name': name, 'degree':degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
			modelGeo = Geomdl(**inputs)
			modelIGA = modelGeo.getIGAParametrization()
			model    = part(modelIGA)

			# Add material 
			mat = thermomat()
			mat.addConductivity(setKprop, isIsotropic=False) 
			mat.addCapacity(setCprop, isIsotropic=False) 

			# Block boundaries
			boundary = step(model._nbctrlpts)
			boundary.add_DirichletTemperature(table=np.array([[1, 1], [0, 0], [0, 0]]))

			# ---------------------
			# Transient model
			# ---------------------
			problem = heatproblem(mat, model, boundary)

			# Create a Dirichlet condition
			Tinout = np.zeros((model._nbctrlpts_total, len(time_list)))
			for i in range(len(time_list)): Tinout[boundary._thdod, i] = 1.0

			# Add external force 
			Fend = problem.eval_heatForce(powden)
			Fext = np.kron(np.atleast_2d(Fend).reshape(-1, 1), sigmoid(time_list))

			# Solve
			lastStep = 40
			problem.solveNLTransientHeatProblemPy(Tinout=Tinout, Fext=Fext[:, :lastStep], 
												time_list=time_list[:lastStep], theta=1.0)
			# # np.savetxt(filename, resPCG)

else:

	for name in name_list:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4.7))

		for i, PCGmethod in enumerate(IterMethods):
			filename = folder + 'ResPCG_' + name + '_' + PCGmethod + '.dat'
			resPCG   = np.loadtxt(filename)

			if PCGmethod   == "WP" : labelmethod = 'w.o.\npreconditioner'
			elif PCGmethod == "C"  : labelmethod = 'Classic FD\nmethod'
			elif PCGmethod == "JMC": labelmethod = 'This work'

			# Print the first 
			newresidue = resPCG[2:, 0]; newresidue = newresidue[newresidue>0]
			ax.semilogy(np.arange(len(newresidue)), newresidue, '-', 
						linewidth=2.5, marker=markerSet[i], label=labelmethod)

			# # Print the last
			# newresidue = resPCG[2:, -1]; newresidue = newresidue[newresidue>0]
			# ax.semilogy(np.arange(len(newresidue)), newresidue, '-', 
			# 			linewidth=2.5, marker=markerSet[i], label=labelmethod)

		ax.legend(bbox_to_anchor=(-0.25, 1.02, 1.25, 0.2), loc='lower left', mode='expand', ncol=3)
		ax.set_xlabel('Number of iterations of BiCGSTAB solver')
		ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
		ax.set_ybound(lower=1e-12, upper=10)

		filename = folder + 'TransientNL_' + name + '.pdf'
		# fig.tight_layout()
		fig.savefig(filename)