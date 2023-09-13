"""
.. Test of elastoplasticity 2D
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformMaxregularKnotvector
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem
from pysrc.lib.thermomecha1D import mechamat1D

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2plasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceSurf_rectangle(P:list):
	nnz = np.size(P, axis=1)
	F = np.zeros((2, nnz))
	F[0, :] = 250
	F[1, :] = 0.0
	return F

# Set global variables
nsteps = 20
matArgs = {'elastic_modulus':2e5, 'elastic_limit':100, 'poisson_ratio':0.3, 
			'plasticLaw': {'name':'linear', 'theta':1, 'Hbar':5e3}}
solverArgs = {'nbIterationsPCG':200, 'PCGThreshold':1e-10, 'PCGmethod': 'TDC', 'NRThreshold':1e-9}

degree, cuts = 8, 6
quadArgs = {'quadrule':'iga', 'type':'leg'}
geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'XY': np.array([[0.0, 0.0], [50.0, 0.0], [50.0, 10.0], [0.0, 10.0]])}
}

# ## 1D
# knotvector = createUniformMaxregularKnotvector(degree, int(2**cuts))
# args  = {'quadArgs': {'degree': degree, 'knotvector': knotvector, 'quadrule':'iga', 'type':'leg'}, 'geoArgs':{'length':50.0}}
# model = mechamat1D(args)

# model.activate_mechanical(matArgs)
# model.add_DirichletCondition(table=[1, 0])
# Fend = np.zeros(model.nbctrlpts); Fend[-1] = 250
# Fext_list  = np.zeros((model.nbctrlpts, nsteps + 1))
# for k in range(0, nsteps+1): Fext_list[:, k] = k/nsteps*Fend
# disp_cp, _, stress, _, _ = model.solve(Fext=Fext_list)


## 2D
blockPrint()
material = mechamat(matArgs)
modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)

# Set Dirichlet boundaries
boundary = boundaryCondition(modelPhy.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
table[0, 0, 0] = 1
table[1, 0, 1] = 1
boundary.add_DirichletDisplacement(table=table)
enablePrint()

# Solve elastic problem
problem = mechaproblem(material, modelPhy, boundary)
problem.addSolverConstraints(solverArgs=solverArgs)
Fext_list = problem.compute_surfForce(forceSurf_rectangle, nbFacePosition=1)[0]
Fend = np.zeros((2, modelPhy.nbctrlpts_total, nsteps+1))
for k in range(1, nsteps+1): Fend[:, :, k] = k/nsteps*Fext_list
displacement = problem.solvePlasticityProblemPy(Fext_list=Fend)[0]
