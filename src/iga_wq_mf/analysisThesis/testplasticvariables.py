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

full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
TRACTION = 400.0
YOUNG, POISSON = 2500, 0.25
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi/2, NBSTEPS)
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e8, 'poisson_ratio': POISSON, 
			'isoHardLaw': {'name':'linear', 'Eiso':500.0}, 
			# 'kineHardLaw':{'parameters':np.array([[500, 0]])}
			}
isReference = True

def forceSurf_infPlate(P:list):
	x = P[0, :]; nnz = np.size(P, axis=1)
	tmp = np.zeros((2, nnz)); tmp[1, :] = x**2-1/4
	F = np.zeros((3, nnz))
	F[2, :] = -TRACTION*(np.min(tmp, axis=0))**2
	return F

def simulate(degree, cuts, quadArgs):
	geoArgs = {'name': 'CB', 'degree': np.array([degree, 3, degree]), 
				'nb_refinementByDirection': np.array([cuts, 0, cuts])}
	blockPrint()
	material = mechamat(MATARGS)
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	meshparam = modelPhy.compute_global_mesh_parameter()

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((3, 2, 3), dtype=int)
	table[0, 0, 0] = 1; table[1, :, 1] = 1; table[2, 0, 2] = 1 
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

degree, cuts = 2, 3
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
problem, displacement, _, internalVars = simulate(degree, cuts, quadArgs)

from pyevtk.vtk import VtkGroup

def run(folder=None):
	assert folder is not None, 'Folder unknown'
	print("Running group...")
	g = VtkGroup(folder)
	for i in range(20):
		g.addFile(filepath = folder + "pls"+str(i)+".vts", sim_time = i)
	g.save()

from pyevtk.vtk import VtkGroup
stress_qp = internalVars.get('stress', None)
alpha_qp = internalVars.get('hardening', None)
plastic_qp = np.where(np.abs(alpha_qp)<1e-6, 0.0, 1.0)
for j, i in enumerate(range(0, NBSTEPS, 4)):
	devstress_qp = computeDeviatoric4All(stress_qp[:, :, i], dim=3)
	vonMises_qp = np.sqrt(3/2)*computeSymTensorNorm4All(devstress_qp, dim=3)
	vonMises_cp = problem.L2projectionCtrlpts(vonMises_qp)
	alpha_cp = problem.L2projectionCtrlpts(alpha_qp[0, :, i])
	plastic_cp = problem.L2projectionCtrlpts(plastic_qp[0, :, i])
	problem.part.exportResultsCP(fields={'stress': vonMises_cp, 'straineq': alpha_cp, 'plastic':plastic_cp}, 
								name='pls3d'+str(j), folder=folder, sampleSize=51)
run(folder=folder)