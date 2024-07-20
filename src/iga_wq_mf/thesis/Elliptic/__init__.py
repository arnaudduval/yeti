from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat, heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import mechaproblem, heatproblem
from pysrc.lib.lib_base import solver

FOLDER2SAVE = os.path.dirname(os.path.realpath(__file__)) + '/results/'
if not os.path.isdir(FOLDER2SAVE): os.mkdir(FOLDER2SAVE)

# Set global variables
GEONAME = 'QA'
TRACTION, RINT, REXT = 1.0, 1.0, 2.0
YOUNG, POISSON = 1e3, 0.3
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e10, 'poisson_ratio':POISSON,
		'isoHardLaw': {'name':'none'}}

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
	b = RINT**2/r_square 
	c = TRACTION*(1.0 + POISSON)*np.sqrt(r_square)/(2*YOUNG)

	disp = np.zeros((2, nnz))
	disp[0, :] = c*(2*(1-POISSON)*np.cos(theta) + b*(4*(1-POISSON)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	disp[1, :] = c*(-2*POISSON*np.sin(theta) + b*(2*(-1 + 2*POISSON)*np.sin(theta) + np.sin(3*theta)) - b**2*np.sin(3*theta))
	
	return disp

def conductivityProperty(args:dict):
	P = args.get('position')
	Kref  = np.array([[1, 0.5],[0.5, 2]])
	Kprop = np.zeros((2, 2, np.size(P, axis=1)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j] 
	return Kprop 

def powerDensity_quartCircle(args:dict):
	""" u = sin(pi*x)*sin(pi*y)*sin(pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	P = args.get('position')
	x = P[0, :]; y = P[1, :]

	f = (3*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- 16*y**2*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 8*x*y*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 4) 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 4) 
	- 8*x**2*np.sin(np.pi*x)*np.sin(np.pi*y)
	)

	return f

def exactTemperature_quartCircle(P: list):
	""" u = sin(pi*x)*sin(pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]

	u = np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 -1)*(x**2 + y**2 - 4)
	return u

def exactDiffTemperature_quartCircle(P: list):
	""" u = sin(pi*x)*sin(pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
		uders = grad(u)
	"""
	x = P[0, :]
	y = P[1, :]

	uders = np.zeros((1, 2, np.size(P, axis=1)))

	uders[0, 0, :] = (2*x*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
				+ 2*x*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
				+ np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
				)

	uders[0, 1, :] = (2*y*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
				+ 2*y*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4)
				+ np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
				)

	return uders


def buildmatrix_el(problem:mechaproblem):

	def krondelta(a,b):
		out = 1 if a==b else 0 
		return out

	def buildsubmatrix(conductivityprop):
		quadrules = problem.part._quadraturerules
		submatrix = sp.csr_matrix((problem.part.nbctrlpts_total, problem.part.nbctrlpts_total))
		for m in range(problem.part.dim):
			beta = np.zeros(problem.part.dim, dtype=int); beta[m] = 1
			tmp1 = sp.kron(quadrules[1]._denseBasis[beta[1]], quadrules[0]._denseBasis[beta[0]]).T
			for l in range(problem.part.dim):
				alpha = np.zeros(problem.part.dim, dtype=int); alpha[l] = 1
				zeta = beta + 2*alpha
				tmp2 = sp.diags(conductivityprop[l, m, :]) @ tmp1
				submatrix += sp.kron(quadrules[1]._denseWeights[zeta[1]], quadrules[0]._denseWeights[zeta[0]]) @ tmp2
		return submatrix

	elastictensor = np.zeros((problem.part.dim,problem.part.dim,problem.part.dim,problem.part.dim))
	for i in range(problem.part.dim):
		for j in range(problem.part.dim):
			for l in range(problem.part.dim):
				for m in range(problem.part.dim):
					elastictensor[i,j,l,m]=(problem.mechamaterial.lame_lambda*krondelta(i,l)*krondelta(j,m)
					+problem.mechamaterial.lame_mu*(krondelta(i,m)*krondelta(j,l)+krondelta(i,j)*krondelta(l,m))
					)

	submatrices = []
	for i in range(problem.part.dim):
		for j in range(problem.part.dim):
			conductivitylike = elastictensor[i,j,:,:]
			prop = np.einsum('ilk,jmk,lm,k->ijk', problem.part.invJ, problem.part.invJ,
						conductivitylike, problem.part.detJ)
			submatrices.append(buildsubmatrix(prop))

	matrix = sp.bmat([[submatrices[0], submatrices[1]],
					[submatrices[2], submatrices[3]]])

	return matrix 

def solvesystem_el(problem:mechaproblem, A, b):
	spilu = sp.linalg.spilu(A)

	def cleanfun(array_in, dod):
		for i in range(problem.part.dim):
			array_in[dod[i]+i*problem.part.nbctrlpts_total] = 0.0
		return
	
	def Afun(array_in):
		array_out = A @ array_in
		return array_out

	def Pfun(array_in):
		array_out = spilu.solve(array_in)
		return array_out

	solv = solver()
	output = solv.GMRES(Afun, b, Pfun=Pfun, cleanfun=cleanfun, dod=problem.boundary.mchdod)
	return output['sol'], output['res']

def simulate_el(degree, cuts, quadArgs=None, preconditioner='JMC'):
	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':RINT, 'Rex':REXT}
				}
	blockPrint()
	if quadArgs is None: quadArgs = {'quadrule':'wq', 'type':1}
	material = mechamat(MATARGS)
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int)
	table[1, 1, 0] = 1; table[1, 0, 1] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	Fext = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=1)[0]

	if preconditioner == 'ilu':
		stiffnessmatrix = buildmatrix_el(problem)
		displacement, residue = solvesystem_el(problem, stiffnessmatrix, np.ravel(Fext))

	else:
		problem._thresLin = 1.e-12; problem._linPreCond = preconditioner
		displacement, residue = problem._solveLinearizedElasticityProblem(Fext=Fext)

	return problem, displacement, residue

def buildmatrix_ht(problem:heatproblem):
	args = {'position':problem.part.qpPhy}
	prop = np.einsum('ilk,jmk,lmk,k->ijk', problem.part.invJ, problem.part.invJ,
				problem.heatmaterial.conductivity(args), problem.part.detJ)
	
	quadrules = problem.part._quadraturerules
	matrix = sp.csr_matrix((problem.part.nbctrlpts_total, problem.part.nbctrlpts_total))
	for j in range(problem.part.dim):
		beta = np.zeros(problem.part.dim, dtype=int); beta[j] = 1
		tmp1 = sp.kron(quadrules[1]._denseBasis[beta[1]], quadrules[0]._denseBasis[beta[0]]).T
		for i in range(problem.part.dim):
			alpha = np.zeros(problem.part.dim, dtype=int); alpha[i] = 1
			zeta = beta + 2*alpha
			tmp2 = sp.diags(prop[i, j, :]) @ tmp1
			matrix += sp.kron(quadrules[1]._denseWeights[zeta[1]], quadrules[0]._denseWeights[zeta[0]]) @ tmp2
	return matrix 

def solvesystem_ht(problem:heatproblem, A, b):
	spilu = sp.linalg.spilu(A)

	def cleanfun(array_in, dod):
		array_in[dod] = 0.0
		return
	
	def Afun(array_in):
		array_out = A @ array_in
		return array_out

	def Pfun(array_in):
		array_out = spilu.solve(array_in)
		return array_out

	solv = solver()
	output = solv.GMRES(Afun, b, Pfun=Pfun, cleanfun=cleanfun, dod=problem.boundary.thdod)
	return output['sol'], output['res']

def simulate_ht(degree, cuts, quadArgs=None, preconditioner='JMC'):
	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':1.0, 'Rex':2.0}
	}
	blockPrint()
	if quadArgs is None: quadArgs = {'quadrule':'wq', 'type':1}
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False)				
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))
	enablePrint()

	# Solve elastic problem
	problem = heatproblem(material, modelPhy, boundary)
	Fext = problem.compute_volForce(powerDensity_quartCircle)

	if preconditioner == 'ilu':
		stiffnessmatrix = buildmatrix_ht(problem)
		temperature, residue = solvesystem_ht(problem, stiffnessmatrix, Fext)
	else:
		problem._thresLin = 1.e-12; problem._linPreCond = preconditioner
		temperature, residue = problem._solveLinearizedSteadyProblem(Fext=Fext, 
									args={'position':problem.part.qpPhy})

	return problem, temperature, residue
