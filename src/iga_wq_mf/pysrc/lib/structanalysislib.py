from lib.__init__ import *
from lib.thermomecha1D import model1D

class Timoshenko(model1D):
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.set_geometry()
		self.activate_mechanical()
		self.compute_mechaproperties()
		self._dod, self._gdod, self._dof = [], [], []
		return
	
	def set_geometry(self):
		super().set_geometry()
		kwargs  = self._kwargs
		name    = kwargs.get('section', 'square').lower()
		if  name == 'square': 
			L = kwargs.get('width', 1.0)
			self._create_square(L)
		elif name == 'rectangle': 
			b = kwargs.get('width', 1.0)
			h = kwargs.get('height', 0.1)
			self._create_rectangle(b, h)
		elif name == 'circle':
			R = kwargs.get('radius', 1.0)
			self._create_circle(R)
		else: raise Warning('Geometry not found')

		return

	def activate_mechanical(self):
		super().activate_mechanical()
		self._shearmodulus = self._elasticmodulus/(2.0*(1.0+self._poissonratio))
		return

	def _create_square(self, L):
		self._area    = L**2.0
		self._inertia = L**4.0/12.0
		self._Ks      = 5.0/6.0
		return

	def _create_circle(self, R):
		self._area    = np.pi*R**2.0/4.0
		self._inertia = np.pi*R**4.0/4.0
		self._Ks      = 5.0/6.0
		return

	def _create_rectangle(self, b, h):
		self._area    = b*h
		self._inertia = b*h**3.0/12.0
		self._Ks      = 5.0/6.0
		return

	def compute_mechaproperties(self):
		self._EA   = self._elasticmodulus*self._area*np.ones(self._nbqp)
		self._EI   = self._elasticmodulus*self._inertia*np.ones(self._nbqp)
		self._GAKs = self._shearmodulus*self._area*self._Ks*np.ones(self._nbqp)
		return 
	
	def _find_freeCntrlPts(self): 
		set_dof = set([i for i in range(3*self._nbctrlpts)])
		set_dod = set(self._dod)
		diff    = set_dof.difference(set_dod)
		set_dof = list(diff)
		return set_dof

	def set_DirichletCondition(self, bound, values=[0., 0., 0.], table=[True, True, True]):
		if bound in self._dod: raise Warning('Enter non-repeated boundary')
		if   bound == 0:
			for i in range(3): 
				if table[i]: self._dod.append(i*self._nbctrlpts)
		elif bound == -1:
			for i in range(3): 
				if table[i]: self._dod.append(self._nbctrlpts + i*self._nbctrlpts - 1)
		else: raise Warning('Only possible block first or last control point')
		for i in range(3): 
			if table[i]: self._gdod.append(values[i])
		self._dof = self._find_freeCntrlPts()
		return

	def compute_stiffness(self, dw):
		" Compute Timoshenko stiffness matrix "
		
		EAdw  = self._EA * dw
		EAdw2 = EAdw * dw

		# Compute sub-matrices
		K11 = compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, self._EA)
		K12 = 0.5*compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, EAdw)
		# % K13 = 0
		K21 = compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, EAdw)
		K22 = (compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, self._GAKs) 
				+ 0.5*compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, EAdw2))
		K23 = compute_sub_advention_matrix(self._basis, self._weights, self._GAKs)
		# % K31 = 0
		K32 = K23.T
		K33 = (compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, self._EI) 
				+ compute_sub_mass_matrix(self._detJ, self._basis, self._weights, self._GAKs))

		# Assemble matrix 
		nr   = self._nbctrlpts
		ind1 = range(nr)
		ind2 = range(nr, 2*nr)
		ind3 = range(2*nr, 3*nr)
		K = np.zeros((3*nr, 3*nr))
		K[np.ix_(ind1, ind1)] = K11
		K[np.ix_(ind1, ind2)] = K12
		K[np.ix_(ind2, ind1)] = K21
		K[np.ix_(ind2, ind2)] = K22
		K[np.ix_(ind2, ind3)] = K23
		K[np.ix_(ind3, ind2)] = K32
		K[np.ix_(ind3, ind3)] = K33

		return K

	def compute_tangentMatrix(self, du, dw):
		" Compute Timoshenko tangent stiffness matrix "
		
		EAdw    = self._EA * dw
		EAdw2   = EAdw * dw
		EAdudw2 = EAdw2 + self._EA*du

		# Compute sub-matrices
		S11 = compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, self._EA)
		S12 = compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, EAdw)
		# % K13 = 0
		S21 = compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, EAdw)
		S22 = (compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, self._GAKs) 
				+ 0.5*compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, EAdw2)
				+ compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, EAdudw2)
		)
		S23 = compute_sub_advention_matrix(self._basis, self._weights, self._GAKs)
		# % K31 = 0
		S32 = S23.T
		S33 = (compute_sub_stiffness_matrix(self._detJ, self._basis, self._weights, self._EI) 
				+ compute_sub_mass_matrix(self._detJ, self._basis, self._weights, self._GAKs))

		# Assemble matrix 
		nr   = self._nbctrlpts
		ind1 = range(nr)
		ind2 = range(nr, 2*nr)
		ind3 = range(2*nr, 3*nr)
		S = np.zeros((3*nr, 3*nr))
		S[np.ix_(ind1, ind1)] = S11
		S[np.ix_(ind1, ind2)] = S12
		S[np.ix_(ind2, ind1)] = S21
		S[np.ix_(ind2, ind2)] = S22
		S[np.ix_(ind2, ind3)] = S23
		S[np.ix_(ind3, ind2)] = S32
		S[np.ix_(ind3, ind3)] = S33

		return S

	def compute_volforce(self, p, q):
		" Compute Timoshenko force vector "
		if np.isscalar(p): p = p*np.ones(self._nbqp)
		if np.isscalar(q): q = q*np.ones(self._nbqp)
		if len(p) != self._nbqp or len(q) != self._nbqp: raise Warning('Not possible')
		F1 = compute_sub_force_vector(self._detJ, self._weights, p)
		F2 = compute_sub_force_vector(self._detJ, self._weights, q)

		nr = self._nbctrlpts
		F  = np.zeros(3*nr)
		F[range(nr)]       = F1 
		F[range(nr, 2*nr)] = F2   
		return F

	def compute_surfForce(self, bound, values=[0, 0, 0]):
		dod = []
		if   bound == 0:
			for i in range(3): dod.append(i*self._nb_ctrlpts)
		elif bound == -1:
			for i in range(3): dod.append(self._nb_ctrlpts + i*self._nb_ctrlpts - 1)
		else: raise Warning('Only possible first or last control point')
		F = np.zeros(3*self._nb_ctrlpts)
		F[dod] = values
		self._dof = self._find_freeCntrlPts()
		return F

	def solve(self, Fext, tol=1e-9, nbIterNL=100):
		" Solves Timoshenko method "

		if not len(self._dod): raise Warning('Add boundary conditions')
		
		nr  = self._nbctrlpts
		sol = np.zeros(3*nr)
		sol[self._dod] = self._gdod
		dof = self._dof
		
		for i in range(nbIterNL):
			# Compute derivative of w
			dwdx = compute_derivative(self._detJ, self._basis, sol[range(nr, 2*nr)])

			# Compute stiffness matrix
			K = self.compute_stiffness(dwdx)

			# Compute dF
			Fint    = K @ sol
			dF      = Fext[dof] - Fint[dof]
			errorNL = np.sqrt(np.dot(dF, dF))
			print('Iteration %s, error: %.5e' %(i, errorNL))
			if errorNL <= tol: break

			# Compute derivative of u
			dudx = compute_derivative(self._detJ, self._basis, sol[range(nr)])

			# Compute stiffness
			S = self.compute_tangentMatrix(dudx, dwdx)[np.ix_(dof, dof)]
			delta_sol = np.linalg.solve(S, dF)
			sol[dof] += delta_sol

		# Select important data
		u     = sol[range(nr)]
		w     = sol[range(nr, 2*nr)]
		theta = sol[range(2*nr, 3*nr)]

		return u, w, theta

def compute_sub_stiffness_matrix(detJ, basis, weights, prop):
	" Computes sub-stiffness matrix in Timoshenko theory "
	coefs = prop / detJ
	S     = weights[-1] @ np.diag(coefs) @ basis[1].T
	return S

def compute_sub_mass_matrix(detJ, basis, weights, prop): 
	" Compute sub-mass matrix in Timoshenko theory "
	coefs = prop * detJ
	M     = weights[0] @ np.diag(coefs) @ basis[0].T
	return M

def compute_sub_advention_matrix(basis, weights, prop):
	" Compute sub-advention matrix in Timoshenko theory "
	A     = weights[-1] @ np.diag(prop) @ basis[0].T
	return A

def compute_sub_force_vector(detJ, weights, prop):
	" Compute sub-force vector in Timoshenko theory "
	coefs = prop * detJ  
	F     = weights[0] @ coefs
	return F

def compute_derivative(detJ, basis, u_ctrlpts):
	" Compute the derivative of u "
	der = basis[1].T @ u_ctrlpts / detJ
	return der
