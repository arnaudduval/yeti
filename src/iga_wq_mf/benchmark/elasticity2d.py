
"""
.. Test of elasticity 2D
.. Infinite plate with a hole under uniaxial traction. 
.. The analytical solution of this problem is given by Timoshenko
.. The convergence curves are traced for IGA-Legendre and IGA-WQ
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import mechaproblem
import pickle

# Set global variables
TRACTION, RINT, REXT = 1.0, 1.0, 2.0
YOUNG, POISSON = 1e3, 0.3
MATARGS = {'elastic_modulus':YOUNG, 
			'elastic_limit':1e10, 
			'poisson_ratio':POISSON,
			'isoHardLaw': {'name':'none'}}
useElastoAlgo = False

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = RINT**2/r_square
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = TRACTION/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = TRACTION/2*3*np.sin(3*theta)*(b**2 - b)
	return F

def exactDisplacement_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	theta = np.arcsin(y/np.sqrt(r_square))
	b = RINT**2/r_square # Already squared
	c = TRACTION*(1.0 + POISSON)*np.sqrt(r_square)/(2*YOUNG)

	disp = np.zeros((2, nnz))
	disp[0, :] = c*(2*(1-POISSON)*np.cos(theta) + b*(4*(1-POISSON)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	disp[1, :] = c*(-2*POISSON*np.sin(theta) + b*(2*(-1 + 2*POISSON)*np.sin(theta) + np.sin(3*theta)) - b**2*np.sin(3*theta))
	
	return disp

degree, cuts = 8, 8
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':RINT, 'Rex':REXT}
				}
material = mechamat(MATARGS)
modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)
meshparam = modelPhy.compute_global_mesh_parameter()

# Set Dirichlet boundaries
boundary = boundaryCondition(modelPhy.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
table[1, 1, 0] = 1; table[1, 0, 1] = 1
boundary.add_DirichletDisplacement(table=table)

# Solve elastic problem
problem = mechaproblem(material, modelPhy, boundary)
if useElastoAlgo:
	Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, 2))
	Fext_list[:, :, 1] = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
	tmp = np.zeros(np.shape(Fext_list))
	problem.solveElastoPlasticityProblem(tmp, Fext_list)
	displacement = tmp[:, :, -1]
else:
	Fext = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
	displacement = problem._solveLinearizedElasticityProblem(Fext)[0]