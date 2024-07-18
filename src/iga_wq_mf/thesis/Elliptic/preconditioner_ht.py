"""
.. Test of steady heat transfer 2D
.. The geometry is a quart of a annulus
.. All the boundary conditions are considered Dirichlet-like 
.. The heat source is computed by f = -grad(k grad T) where T is a given function
.. We compute the relative L2 norm using the exact solution T
.. The convergence curves are traced for IGA-Legendre and IGA-WQ
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import heatproblem
from pysrc.lib.lib_base import solver

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/preconditioner/'
if not os.path.isdir(folder): os.mkdir(folder)

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
	x = P[0, :]
	y = P[1, :]

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

def buildmatrix(problem:heatproblem):
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

def solvesystem(problem:heatproblem, A, b):
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
	_, output = solv.GMRES(Afun, b, Pfun=Pfun, cleanfun=cleanfun, dod=problem.boundary.thdod)
	return output['res']

def simulate(degree, cuts, preconditioner='JMC'):
	geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'Rin':1.0, 'Rex':2.0}
	}
	blockPrint()
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False)				
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs={'quadrule': 'wq', 'type': 1})

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))
	enablePrint()

	# Solve elastic problem
	problem = heatproblem(material, modelPhy, boundary)
	Fext = problem.compute_volForce(powerDensity_quartCircle)

	if preconditioner == 'ilu':
		stiffnessmatrix = buildmatrix(problem)
		residue = solvesystem(problem, stiffnessmatrix, Fext)
	else:
		problem._thresLin = 1.e-12; problem._linPreCond = preconditioner
		_, residue = problem._solveLinearizedSteadyProblem(Fext=Fext, 
									args={'position':problem.part.qpPhy})

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
# fig.savefig(folder+'preconditioner_ht'+'.pdf')


for degree in range(4, 7):
	for cuts in range(6, 9):
		start = time.process_time()
		problem, residue = simulate(degree, cuts, preconditioner='JMC')
		stop = time.process_time()
		print('%d, %d, %.2f, %d' %(degree, cuts, stop-start, len(residue[residue>0.0])))