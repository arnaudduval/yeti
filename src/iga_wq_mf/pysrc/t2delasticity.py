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

# Set global variables
degree, cuts = 4, 7
name = 'QA'

# Create model 
geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 'extra':{'Rin':5.e2, 'Rex':1.e3}}
quadArgs  = {'quadrule': 'wq', 'type': 1}

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
	ref  = np.array([2.5e1, 0.0])
	prop = np.zeros((2, np.size(P, axis=1)))
	for i in range(2): prop[i, :] = ref[i] 
	return prop

Fext = problem.eval_surfForce(forceSurfFun, nbFacePosition=1)
if name == 'SQ': 
	Fext += problem.eval_volForce(forceVolFun1)
elif name == 'QA': 
	Fext += problem.eval_volForce(forceVolFun2)

# -------------
# ELASTICITY
# -------------
fig, ax = plt.subplots()

# Solve in fortran 
for i, [methodPCG, label] in enumerate(zip(['WP', 'C', 'JMC'], 
							['w.o. preconditioner', 'Fast diag. (FD)', 'This work'])):
    problem.addSolverConstraints(solverArgs={'PCGmethod': methodPCG})
    displacement, resPCG = problem.solveElasticityProblemFT(Fext=Fext)
    resPCG = resPCG[resPCG>0]
    ax.semilogy(np.arange(len(resPCG)), resPCG, '-', label=label, marker=markerSet[i])

model.exportResults(u_ctrlpts=displacement, nbDOF=2, folder=folder)
ax.set_ybound(lower=1e-12, upper=1e1)
ax.legend()
ax.set_xlabel('Number of iterations of BiCGSTAB solver')
ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
fig.tight_layout()
fig.savefig(folder + name + 'ElasRes.png')