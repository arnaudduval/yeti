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
from pysrc.lib.lib_material import mechamat, block_dot_product, clean_dirichlet
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem
from pysrc.lib.lib_base import solver

full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path)
if not os.path.isdir(folder): os.mkdir(folder)


# Set global variables
TRACTION, RINT, REXT = 1.0, 1.0, 2.0
YOUNG, POISSON = 1e3, 0.3
GEONAME = 'QA'
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e10, 'poisson_ratio':POISSON,
		'isoHardLaw': {'Isoname':'none'}}
SOLVERARGS  = {'nIterKrylov':250, 'thresholdKrylov':1e-10, 'KrylovPreconditioner': 'JMC'}
isReference = False

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

def simulate(degree, cuts, quadArgs, useElastoAlgo=False):
	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':RINT, 'Rex':REXT}
				}
	blockPrint()
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
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	problem.addSolverConstraints(solverArgs=SOLVERARGS)
	if useElastoAlgo:
		Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, 2))
		Fext_list[:, :, 1] = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
		tmp = np.zeros(np.shape(Fext_list))
		problem.solveElastoPlasticityProblem(tmp, Fext_list)
		displacement = tmp[:, :, -1]
	else:
		Fext = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
		displacement = problem._solveLinearizedElasticityProblem(Fext)[0]
	return problem, displacement, meshparam

# degree_list = np.array([1])
# cuts_list   = np.arange(1, 6)

# fig, ax = plt.subplots(figsize=(8, 7))
# for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
# 	quadArgs = {'quadrule': quadrule, 'type': quadtype}
# 	error_list = np.ones(len(cuts_list))

# 	for i, degree in enumerate(degree_list):
# 		meshparam = np.ones(len(cuts_list))
# 		color = COLORLIST[i]
# 		for j, cuts in enumerate(cuts_list):
# 			problem, displacement, meshparam[j] = simulate(degree, cuts, quadArgs, useElastoAlgo=False)
# 			error_list[j] = problem.normOfError(displacement, isRelative=False, 
# 							normArgs={'type':'L2', 'exactFunction':exactDisplacement_infPlate})
# 		print(error_list)

geoArgs = {'name': GEONAME, 'degree': 1*np.ones(3, dtype=int), 
			'nb_refinementByDirection': 4*np.ones(3, dtype=int), 
			'extra':{'Rin':RINT, 'Rex':REXT}
			}
blockPrint()
material = mechamat(MATARGS)
modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy_iga = part(modelIGA, quadArgs={'quadrule': 'iga'})
modelPhy_wq1 = part(modelIGA, quadArgs={'quadrule': 'wq', 'type':1})
modelPhy_wq2 = part(modelIGA, quadArgs={'quadrule': 'wq', 'type':2})

# Set Dirichlet boundaries
boundary = boundaryCondition(modelPhy_iga.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
table[1, 1, 0] = 1; table[1, 0, 1] = 1
boundary.add_DirichletDisplacement(table=table)
enablePrint()

problem_wq1 = mechaproblem(material, modelPhy_wq1, boundary)
problem_wq1.addSolverConstraints(solverArgs=SOLVERARGS)
Fext_wq1 = problem_wq1.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
disp_wq1 = problem_wq1._solveLinearizedElasticityProblem(Fext_wq1)[0]

problem_wq2 = mechaproblem(material, modelPhy_wq2, boundary)
problem_wq2.addSolverConstraints(solverArgs=SOLVERARGS)
Fext_wq2 = problem_wq2.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
disp_wq2 = problem_wq2._solveLinearizedElasticityProblem(Fext_wq2)[0]

problem_iga = mechaproblem(material, modelPhy_iga, boundary)
problem_iga.addSolverConstraints(solverArgs=SOLVERARGS)
Fext_iga = problem_iga.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]
disp_iga = problem_iga._solveLinearizedElasticityProblem(Fext_iga)[0]

error1 = abs(block_dot_product(disp_wq1-disp_iga, disp_wq1-disp_iga)/block_dot_product(disp_iga, disp_iga))
error2 = abs(block_dot_product(disp_wq2-disp_iga, disp_wq2-disp_iga)/block_dot_product(disp_iga, disp_iga))
error3 = abs(block_dot_product(disp_wq1-disp_wq2, disp_wq1-disp_wq2)/block_dot_product(disp_iga, disp_iga))
print('%.3e, %.3e, %.3e'%(error1, error2, error3))

solv = solver()
for niters in range(1, 170, 20):
	solv._itersLin = niters
	outWQ1 = solv.BiCGSTAB(problem_wq1.compute_mfStiffness, Fext_wq1, 
				dotfun=block_dot_product, cleanfun=clean_dirichlet, dod=boundary.mchdod)

	outWQ2 = solv.BiCGSTAB(problem_wq2.compute_mfStiffness, Fext_wq2, 
				dotfun=block_dot_product, cleanfun=clean_dirichlet, dod=boundary.mchdod)

	outIGA = solv.BiCGSTAB(problem_iga.compute_mfStiffness, Fext_iga, 
				dotfun=block_dot_product, cleanfun=clean_dirichlet, dod=boundary.mchdod)

	solwq1 = outWQ1['sol']; solwq2 = outWQ2['sol']
	error1 = abs(block_dot_product(solwq1-disp_iga, solwq1-disp_iga)/block_dot_product(disp_iga, disp_iga))
	error2 = abs(block_dot_product(solwq2-disp_iga, solwq2-disp_iga)/block_dot_product(disp_iga, disp_iga))
	error3 = abs(block_dot_product(solwq1-solwq2, solwq1-solwq2)/block_dot_product(disp_iga, disp_iga))
	print('%d, %.3e, %.3e, %.3e'%(niters, error1, error2, error3))

# fig, ax = plt.subplots()
# ax.semilogy(resKryIGA, label='IGA')
# ax.semilogy(resKryWQ1, label='WQ1')
# ax.semilogy(resKryWQ2, label='WQ2')
# ax.legend()
# fig.savefig(folder+'resKrylov')
