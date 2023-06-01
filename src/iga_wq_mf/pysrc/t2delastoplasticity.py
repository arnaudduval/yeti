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
degree, cuts = 2, 3
name = 'QA'

# Create model 
geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
quadArgs  = {'quadrule': 'wq', 'type': 1}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
model    = part(modelIGA, quadArgs=quadArgs)

# Add material 
matArgs  = {'density': 7800, 'elastic_modulus': 1e9, 'poisson_ratio': 0.3, 'elastic_limit': 500e9}
material = mechamat(matArgs)

# Set Dirichlet boundaries
boundary = boundaryCondition(model.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
table[0, 0, 0] = 1
table[0, 1, 1] = 1
boundary.add_DirichletDisplacement(table=table)

# Elasticity problem
problem = mechaproblem(material, model, boundary)

# Set Neumann boundaries
def forceFun(P:list):
	x = P[0, :]
	y = P[1, :]
	ref  = np.array([1e8, 2e8])
	prop = np.zeros((2, len(x)))
	for i in range(2): prop[i, :] = ref[i] 
	return prop

Fsurf = problem.eval_surfForce(forceFun, nbFacePosition=3)

# -------------
# ELASTICITY
# -------------
# Solve in fortran 
displacement, resPCG = problem.solveElasticityProblemFT(Fext=Fsurf, methodPCG='C')


# fig, ax = plt.subplots()

# # Solve in fortran 
# for methodPCG, label in zip(['WP', 'C', 'JMC'], ['w.o. preconditioner', 'Fast diag. (FD)', 'This work']):
#     displacement, resPCG = problem.solveElasticityProblemFT(Fext=Fsurf, methodPCG=methodPCG)
#     resPCG = resPCG[resPCG>0]
#     ax.semilogy(np.arange(len(resPCG)), resPCG, label=label)

# ax.set_ybound(lower=1e-8, upper=10)
# ax.legend()
# ax.set_xlabel('Number of iterations of BiCGSTAB solver')
# ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
# fig.tight_layout()
# fig.savefig(folder + name + 'ElasRes.png')