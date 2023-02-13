"""
This file contains Discontinue Galerkin Finite Elements Methods algorithms
The strong form of the equation to be solved looks like:
du/dt + a df/dx = 0
where f depends on u which depends on x and t
The discretization of the weak form looks like:
duj/dt + a (f_j+1 - f_j-1)/(2*Deltaj) = 0
for each element j of size Deltaj. Here f_k (for k = j-1 and j+1) 
represents the flux at one of the boundaries (left if j-1 and right if j+1)   
"""

from others.__init__ import *

def gaussTable(order):
	" Computes Gauss weights and positions in isoparametric space for a given degree "

	if order == 1:
		pos = [0.0]
		wgt = [2.0]
	elif order == 2:
		pos = [ -0.57735026918962576,
				0.57735026918962576]

		wgt = [ 1.0,
				1.0]
	elif order == 3:
		pos = [ -0.77459666924148337,
				0.0,
				0.77459666924148337]

		wgt = [ 5.0 / 9.0,
				8.0 / 9.0,
				5.0 / 9.0]
	elif order == 4:
		pos = [ -0.8611363115940526,
				-0.3399810435848563,
				0.3399810435848563, 
				0.8611363115940526]

		wgt = [ 0.3478548451374539,
				0.6521451548625461,
				0.6521451548625461,
				0.3478548451374539]
	else: raise Warning('Not degree found')

	# Change type of arrays
	pos, wgt = np.array(pos), np.array(wgt)

	return pos, wgt

def find_span(array, x, threshold=1e-8):
	"Find the span of an array where the value x is located"
	isnode    = False
	isspecial = False

	# Special cases
	if (abs(array[0]-x)<threshold): 
		isnode = True; isspecial = True; span = 0
	if (abs(array[-1]-x)<threshold):
		isnode = True; isspecial = True; span = len(array)-2

	# Node cases
	if isspecial is False:
		for i, node in enumerate(array[1:-1]):
			if abs(node-x) < threshold:
				isnode = True; span = [i, i+1]; break

	# General case
	if isnode is False:
		span = 1
		size = len(array)
		while (span<size) & (array[span]-x<threshold):
			span = span + 1
		span = span -1

	return np.atleast_1d(span)

def dgEvalPolyBasis(xis, degree):
	""" Returns the values of the basis functions evaluated at xi.
		These functions are defined in the parametric element [-1, 1].
	"""
	
	b0, b1 = [], []
	for xi in xis:
		temp = [1.0, 0.5*xi]
		temp.extend([(xi/2.0)**n for n in range(2, degree+1)])
		b0.append(temp)
		temp = [0.0, 0.5]
		temp.extend([0.5*n*(0.5*xi)**(n-1) for n in range(2, degree+1)])
		b1.append(temp)
	
	return [np.array(b0), np.array(b1)]

def dgGetMass(degree, nodes):

	def dgGetElementMass(degree, nodes, j):
		"Considers a uniform discretization. Maybe later it will be adapted to use knotvector"

		xij, xfj = nodes[j], nodes[j+1]
		deltaj   = xfj - xij
		xj       = 0.5*(xfj + xij)

		# Get quadrature points and weights in parametric space
		gaussPos, gaussWeight = gaussTable(degree+1)
		b0, b1 = dgEvalPolyBasis(gaussPos, degree)

		# Get quadrature points in physical space
		x = []
		for xi in gaussPos:
			temp = 0.5*xi*deltaj + xj
			x.append(temp)
		x = np.array(x)
		
		# Compute elementary mass
		mass  = (b0.T @ np.diag(gaussWeight) @ b0)*deltaj*0.5

		return mass

	nbel  = len(nodes)-1
	mass  = np.zeros((nbel*(degree+1), nbel*(degree+1)))
	for j in range(nbel):
		massel = dgGetElementMass(degree, nodes, j)
		rang = range(j*(degree+1), (j+1)*(degree+1))
		mass[np.ix_(rang, rang)] += massel
	return mass

def dgGenFunFUx(xs, kwargs:dict):
	"Defines the function F which depends on U. At the same time U depends on x"
	
	def fun_Ux(degree, xi, xf, x, Udof):
		# Interpolate x in parametric space
		xj = 0.5*(xi + xf)
		delta = xf - xi
		xi = 2.0*(x - xj)/delta
		b0 = dgEvalPolyBasis([xi], degree)[0]
		Ufun = Udof @ b0
		return Ufun
	
	degree   = kwargs.get('degree')
	nodes    = kwargs.get('nodes')
	Uctrlpts = kwargs.get('Uctrlpts')
	funFU    = kwargs.get('funFU')

	FUx = []
	for x in xs:
		span = find_span(nodes, x)
		sumo = 0.0
		for sp in span:
			xi = nodes[sp]
			xf = nodes[sp+1]
			Uctrlpt = Uctrlpts[sp*(degree+1):(sp+1)*(degree+1)]
			Uinterp = fun_Ux(degree, xi, xf, x, Uctrlpt)
			sumo   += funFU(Uinterp)
		sumo = sumo/len(span)
		FUx.append(sumo)

	return np.atleast_1d(FUx)

def dgGetForce(degree, nodes, Uctrlpts, funFU):

	def dgGetElementForce(degree, nodes, j, kwargs):
		"Considers a uniform discretization. Maybe later it will be adapted to use knotvector"

		xij, xfj = nodes[j], nodes[j+1]
		deltaj   = xfj - xij
		xj       = 0.5*(xfj + xij)

		# Get quadrature points and weights in parametric space
		gaussPos, gaussWeight = gaussTable(degree+1)
		b0, b1 = dgEvalPolyBasis(gaussPos, degree)

		# Get quadrature points in physical space
		x = []
		for xi in gaussPos:
			temp = 0.5*xi*deltaj + xj
			x.append(temp)
		x = np.atleast_1d(x)
		
		# Compute elementary mass and force
		force = b1.T @ (gaussWeight*dgGenFunFUx(x, kwargs))*2.0/deltaj
		force = np.atleast_2d(force)

		# Add elements of the forces
		basisLeft = dgEvalPolyBasis([-1], degree)[0]
		basisRight = dgEvalPolyBasis([1], degree)[0]

		force += dgGenFunFUx([xij], kwargs)*basisLeft - dgGenFunFUx([xfj], kwargs)*basisRight

		return force

	nbel  = len(nodes)-1
	force = np.zeros(nbel*(degree+1)) 

	kwargs = {'degree': degree, 'nodes': nodes, 'Uctrlpts': Uctrlpts, 'funFU':funFU}
	for j in range(nbel):
		forceel = dgGetElementForce(degree, nodes, j, kwargs)
		rang = range(j*(degree+1), (j+1)*(degree+1))
		force[rang] += np.ndarray.flatten(forceel)
	return force

def dgGetU0(degree, nodes, funU0):

	def dgGetElementU0(degree, nodes, j, funU0):
		"Considers a uniform discretization. Maybe later it will be adapted to use knotvector"

		xij, xfj = nodes[j], nodes[j+1]
		deltaj   = xfj - xij
		xj       = 0.5*(xfj + xij)

		# Get quadrature points and weights in parametric space
		gaussPos, gaussWeight = gaussTable(degree+1)
		b0, b1 = dgEvalPolyBasis(gaussPos, degree)

		# Get quadrature points in physical space
		x = []
		for xi in gaussPos:
			temp = 0.5*xi*deltaj + xj
			x.append(temp)
		x = np.array(x)
		
		# Compute elementary U0
		U0  = (b0.T @ np.diag(gaussWeight)*funU0(x))*deltaj*0.5
		return U0
	
	nbel  = len(nodes)-1
	U0  = np.zeros(nbel*(degree+1))
	for j in range(nbel):
		U0el = dgGetElementU0(degree, nodes, j, funU0)
		rang = range(j*(degree+1), (j+1)*(degree+1))
		U0[rang] += U0el
	return U0

def rungekutta3(Utilde0, Tspan, kwargs:dict):
	""" This algorithm solves a partial diferential equation given by:
		M dU/dt = L(U), M U(0) = Utilde0 in th interval [0, Tspan].
		It uses Runge-Kutta order 3 method
	"""
	nbSteps = kwargs.get('nsteps', 2)  # At least 2 steps
	degree  = kwargs.get('degree')
	nodes   = kwargs.get('nodes')
	funFU   = kwargs.get('funFU')

	# Get mass matrix
	M = dgGetMass(degree, nodes)

	# Get real U0
	U  = np.zeros((nbel*(degree+1), nbSteps))
	U[:, 0] = np.linalg.solve(M, Utilde0)

	if nbSteps < 2: raise Warning('At least 2 steps')
	dt = Tspan/(nbSteps-1)
	for step in range(1, nbSteps):
		U_n0 = U[:, step-1]
		FU_n_tilde = dgGetForce(degree, nodes, U_n0, funFU)
		FU_n = np.linalg.solve(M, FU_n_tilde)
		w1   = U_n0 + dt*FU_n
		Fw1_tilde  = dgGetForce(degree, nodes, w1, funFU)
		Fw1  = np.linalg.solve(M, Fw1_tilde)
		w2   = 0.75*U_n0 + 0.25*(w1 + dt*Fw1)
		Fw2_tilde  = dgGetForce(degree, nodes, w2, funFU)
		Fw2  = np.linalg.solve(M, Fw2_tilde)
		U_n1 = 1.0/3.0*U_n0 + 2.0/3.0*(w2 + dt*Fw2)
		U[:, step] = U_n1

	return U

def interpolate(degree, nodes, Uctrlpts):

	return

# Problem data 
length = 1.0
alpha  = 2.0

# Galerkin discretisation
nbSteps = 10 
Tspan   = 0.2
degree, nbel = 2, 4
nodes  = np.linspace(0, length, nbel+1)

# Define function of F(U) and U(x)
funFU  = lambda u: alpha*u
funU0  = lambda x: x*(x-1)
Utilde0  = dgGetU0(degree, nodes, funU0)
kwargs   = {'nsteps': nbSteps, 'degree': degree, 'nodes': nodes, 'funFU': funFU}
Uctrlpts = rungekutta3(Utilde0, Tspan, kwargs) 
