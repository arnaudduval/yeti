"""
This file contains Discontinue Galerkin Finite Elements Methods algorithms
The strong form of the equation to be solved looks like:
du/dt + a df/dx = 0
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

def genfun_FUx(xs, kwargs:dict):
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
			Uctrlpt = Uctrlpts[sp, :]
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
		force = b1.T @ (gaussWeight*genfun_FUx(x, kwargs))*2.0/deltaj
		force = np.atleast_2d(force)

		# Add elements of the forces
		basisLeft = dgEvalPolyBasis([-1], degree)[0]
		basisRight = dgEvalPolyBasis([1], degree)[0]

		force += genfun_FUx([xij], kwargs)*basisLeft - genfun_FUx([xfj], kwargs)*basisRight

		return force

	nbel  = len(nodes)-1
	force = np.zeros(nbel*(degree+1)) 

	kwargs = {'degree': degree, 'nodes': nodes, 'Uctrlpts': Uctrlpts, 'funFU':funFU}
	for j in range(nbel):
		forceel = dgGetElementForce(degree, nodes, j, kwargs)
		rang = range(j*(degree+1), (j+1)*(degree+1))
		force[rang] += np.ndarray.flatten(forceel)
	return force

degree, nbel = 2, 4
length = 1.0
nodes  = np.linspace(0, length, nbel+1)
span = find_span(nodes, 1.0)
