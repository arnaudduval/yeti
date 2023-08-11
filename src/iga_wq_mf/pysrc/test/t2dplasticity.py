"""
.. Test of elastoplasticity 2D
.. We test how elastoplasticity module works
.. Joaquin Cornejo 
"""
import pickle
from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceSurf_infPlate(P:list):
	Tx, a = 5e7, 1.0
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = a**2/r_square # Already squared
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = Tx/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = Tx/2*3*np.sin(3*theta)*(b**2 - b)
	return F

# Set global variables
nsteps = 15
E, nu = 2e11, 0.3
matArgs    = {'elastic_modulus':E, 'elastic_limit':8e15, 'poisson_ratio': nu, 
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-10, 'PCGmethod': 'TDC', 'NRThreshold':1e-5}

degree, cuts = 8, 7
quadArgs = {'quadrule': 'iga', 'type': 'leg'}

geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'Rin':1.0, 'Rex':4.0}
}
blockPrint()
material = mechamat(matArgs)
modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)
with open(folder + 'refpart.pkl', 'wb') as outp:
    pickle.dump(modelPhy, outp, pickle.HIGHEST_PROTOCOL)

# Set Dirichlet boundaries
boundary = boundaryCondition(modelPhy.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
table[1, 1, 0] = 1
table[1, 0, 1] = 1
boundary.add_DirichletDisplacement(table=table)
enablePrint()

# Solve elastic problem
problem = mechaproblem(material, modelPhy, boundary)
problem.addSolverConstraints(solverArgs=solverArgs)
Fend = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, nsteps+1))
for i in range(1, nsteps+1): Fext_list[:, :, i] = i/nsteps*Fend
displacement = problem.solvePlasticityProblemPy(Fext_list=Fext_list)[0]
np.save(folder+'u_ref', displacement)
