from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import mechaproblem
from pysrc.lib.lib_base import solver

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/preconditioner/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
TRACTION, RINT, REXT = 1.0, 0.5, 1.0
YOUNG, POISSON = 2., 0.0
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
	b = RINT**2/r_square # Already squared
	c = TRACTION*(1.0 + POISSON)*np.sqrt(r_square)/(2*YOUNG)

	disp = np.zeros((2, nnz))
	disp[0, :] = c*(2*(1-POISSON)*np.cos(theta) + b*(4*(1-POISSON)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	disp[1, :] = c*(-2*POISSON*np.sin(theta) + b*(2*(-1 + 2*POISSON)*np.sin(theta) + np.sin(3*theta)) - b**2*np.sin(3*theta))
	
	return disp

def buildmatrix(problem:mechaproblem):

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

def solvesystem(problem:mechaproblem, A, b):
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
	_, output = solv.GMRES(Afun, b, Pfun=Pfun, cleanfun=cleanfun, dod=problem.boundary.mchdod)
	return output['res']

def simulate(degree, cuts, preconditioner='JMC'):
	geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':RINT, 'Rex':REXT}
				}
	material = mechamat(MATARGS)

	blockPrint()
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs={'quadrule':'wq', 'type':1})

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
		stiffnessmatrix = buildmatrix(problem)
		residue = solvesystem(problem, stiffnessmatrix, np.ravel(Fext))

	else:
		problem._thresLin = 1.e-12; problem._linPreCond = preconditioner
		_, residue = problem._solveLinearizedElasticityProblem(Fext=Fext)
		
	return problem, residue

# degree, cuts = 6, 6
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5, 5.5))
# for j, preconditioner in enumerate(['WP', 'ilu', 'C', 'JMC']):
# 	start = time.process_time()
# 	problem, residue = simulate(degree, cuts, preconditioner=preconditioner)
# 	stop = time.process_time()
# 	print('time:%.2e'%(stop-start))

# 	if preconditioner == 'WP': labelfig = 'w.o. preconditioner'
# 	elif preconditioner == 'ilu': labelfig = 'Incomplete LU'
# 	elif preconditioner == 'C' : labelfig = 'Classic FD'
# 	elif preconditioner == 'JMC' : labelfig = 'This work'
# 	ax.semilogy(residue, marker=MARKERLIST[j], label=labelfig)

# ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.2), loc='upper center')
# ax.set_ylabel('Relative residue')
# ax.set_xlabel('Number of iterations (GMRES)')
# ax.set_ylim([1e-12, 1e1])
# ax.set_xlim([0, 100])
# fig.savefig(folder+'preconditioner_el'+'.pdf')


for degree in range(4, 7):
	for cuts in range(6, 9):
		start = time.process_time()
		problem, residue = simulate(degree, cuts, preconditioner='JMC')
		stop = time.process_time()
		print('%d, %d, %.2f, %d' %(degree, cuts, stop-start, len(residue[residue>0.0])))