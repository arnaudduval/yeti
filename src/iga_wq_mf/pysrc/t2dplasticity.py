"""
.. Test of elastoplasticity 2D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : Pa (210e9)
..      - Length : m
..      - Force  : N
..      - Mass   : kg 
..      - Density: kg/m^3 (7.8e3)
..      - Gravity: m/s^2 (9.8)
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

# Set global variables
degree, cuts = 3, 5
name = 'SQ'
dataExist = True
IterMethods = ['C', 'JMC']

if not dataExist:
	for PCGmethod in IterMethods:
		filename = folder + 'ResPCG_' + name + '_' + PCGmethod + '.dat'  

		# Create model 
		geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
		quadArgs  = {'quadrule': 'iga', 'type': 'leg'}

		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		model    = part(modelIGA, quadArgs=quadArgs)

		# Add material 
		matArgs    = {'density': 7800, 'elastic_modulus':1e9, 'poisson_ratio': 0.3, 'elastic_limit':5e6,
					'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}

		material = mechamat(matArgs)

		# Set Dirichlet boundaries
		boundary = boundaryCondition(model.nbctrlpts)
		table = np.zeros((2, 2, 2), dtype=int)
		table[0, 0, 0] = 1
		table[1, 0, 1] = 1
		boundary.add_DirichletDisplacement(table=table)

		# Elasticity problem
		problem = mechaproblem(material, model, boundary)

		def forceSurfFun(P:list):
			x = P[0, :]
			y = P[1, :]
			ref  = np.array([1.e7, 0.0])
			prop = np.zeros((2, len(x)))
			for i in range(2): prop[i, :] = ref[i] 
			return prop

		Fsurf = problem.eval_surfForce(forceSurfFun, nbFacePosition=1)
		nbStep = 6; dt = 1/nbStep
		Fext   = np.zeros((*np.shape(Fsurf), nbStep+1))
		for i in range(1, nbStep+1): Fext[:, :, i] = i*dt*Fsurf

		displacement, resPCG = problem.solvePlasticityProblemPy(Fext=Fext, thresholdNR=1e-5, methodPCG=PCGmethod)
		np.savetxt(filename, resPCG)

else: 
	fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))

	for i, PCGmethod in enumerate(IterMethods):
		filename = folder + 'ResPCG_' + name + '_' + PCGmethod  + '.dat'
		resPCG   = np.loadtxt(filename)

		if   PCGmethod == "C"  : labelmethod = 'Classic FD\nmethod'
		elif PCGmethod == "JMC": labelmethod = 'This work'

		ind = np.where(resPCG[:, 0]==1)
		newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
		ax1.semilogy(np.arange(len(newresidue)), newresidue, '-', color='k', linewidth=0.5)
		
		newresidue = resPCG[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
		ax2.semilogy(np.arange(len(newresidue)), newresidue, '-', color='k', linewidth=0.5)

		maxStep = int(np.max(resPCG[:, 0]))
		for j in range(2, maxStep):
			opacity = (maxStep - j + 1)*1.0/(maxStep - 1)
			ind = np.where(resPCG[:, 0]==j)

			newresidue = resPCG[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
			ax1.semilogy(np.arange(len(newresidue)), newresidue, '-', alpha=opacity, 
						color=colorSet[i], linewidth=0.5)
			
			newresidue = resPCG[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
			ax2.semilogy(np.arange(len(newresidue)), newresidue, '-', alpha=opacity, 
						color=colorSet[i], linewidth=0.5)

	ax1.set_title('First NR iterations')
	ax2.set_title('Last NR iterations')
	for ax in [ax1, ax2]:
		ax.set_xlabel('Number of iterations of BiCGSTAB solver')
		ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
		ax.set_ybound(lower=1e-12, upper=10)

	filename = folder + 'PlasticityNL_' + name + '.png'
	fig.tight_layout()
	fig.savefig(filename)