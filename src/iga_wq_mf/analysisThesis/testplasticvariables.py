# from pysrc.lib.__init__ import *
# from pysrc.lib.lib_material import mechamat

# # Set global variables
# YOUNG, POISSON = 2400, 0.2
# H, beta = 100, 0.3
# MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':300, 'poisson_ratio': POISSON, 
# 			'isoHardLaw': {'name':'linear', 'Eiso':(1-beta)*H}, 
# 			'kineHardLaw':{'parameters':np.array([[2/3*beta*H, 0]])}
# 			}

# material = mechamat(MATARGS)
# strain = np.zeros((6, 1))
# strain[0, 0] = 0.125 + 0.1; strain[1, 0] = -0.025 - 0.02; strain[2, 0] = -0.025 - 0.02
# pls_n0 = np.zeros((6, 1))
# a_n0   = np.zeros((1, 1))
# b_n0   = np.zeros((1, 6, 1))
# out, _ = material.J2returnMappingAlgorithm3D(strain, pls_n0, a_n0, b_n0)
# lame_mu = material.lame_mu
# dgamma, norm_eta_trial = 0.0948, 180*np.sqrt(6)
# c1 = 4*lame_mu**2/(2*lame_mu+2/3*H)
# c2 = 4*lame_mu**2*dgamma/norm_eta_trial
# print(c1/(2*lame_mu), c2/(2*lame_mu))

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import *
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Set global variables
TRACTION = 400.0
YOUNG, POISSON = 2500, 0.25
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi/2, NBSTEPS)
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':5, 'poisson_ratio': POISSON, 
			'isoHardLaw': {'name':'linear', 'Eiso':500.0}, 
			# 'kineHardLaw':{'parameters':np.array([[500, 0]])}
			}
isReference = True

def forceSurf_infPlate(P:list):
	x = P[0, :]; nnz = np.size(P, axis=1)
	tmp = np.zeros((2, nnz)); tmp[1, :] = (x-0.5)**2-1/4
	F = np.zeros((3, nnz))
	F[2, :] = -TRACTION*(np.min(tmp, axis=0))**2
	return F

def simulate(degree, cuts, quadArgs):
	geoArgs = {'name': 'CB', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
	blockPrint()
	material = mechamat(MATARGS)
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	meshparam = modelPhy.compute_global_mesh_parameter()

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((3, 2, 3), dtype=int); table[2, 0, :] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	problem._thresNL = 1e-6; problem._itersNL = 100
	Fref = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=5)[0]
	Fext_list = np.zeros((3, modelPhy.nbctrlpts_total, NBSTEPS))
	for k in range(len(TIME_LIST)): Fext_list[:, :, k] = np.sin(TIME_LIST[k])*Fref
	displacement = np.zeros(np.shape(Fext_list))
	_, internalVars = problem.solveElastoPlasticityProblem(displacement, Fext_list)
	return problem, displacement, meshparam, internalVars

degree, cuts = 1, 4
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
problem, displacement, _, internalVars = simulate(degree, cuts, quadArgs)
