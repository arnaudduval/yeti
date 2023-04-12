# from lib.__init__ import *

# class Timoshenko():
# 	def __init__(self, name=None, properties=dict):

# 		self.read_material(prop=properties)
# 		self.read_bspline(prop=properties)

# 		# Select geometry
# 		self._Jpq = properties.get('length', 1.0)
# 		if  name == 'square': 
# 			L = properties.get('width', 1.0)
# 			self.create_square(L)
# 		elif name == 'rectangle': 
# 			b = properties.get('width', 1.0)
# 			h = properties.get('height', 0.1)
# 			self.create_rectangle(b, h)
# 		elif name == 'circle':
# 			R = properties.get('radius', 1.0)
# 			self.create_circle(R)
# 		else: raise Warning('Geometry not found')

# 		self.set_bspline_basis_weights()
# 		self._EA, self._EI, self._GAKs = None, None, None
# 		self._dod, self._gdod          = [], []

# 		return

# 	def read_material(self, prop:dict):
# 		self._Emod = prop.get('E', 1.e8)
# 		self._nu   = prop.get('nu', 0.3)
# 		self._Gmod = self._Emod/(2.0*(1.0+self._nu))
# 		return

# 	def read_bspline(self, prop:dict):
# 		self._degree     = prop.get('degree', 3)
# 		self._knotvector = prop.get('knotvector', np.array([0, 0, 0, 0, 1, 1, 1, 1]))
# 		self._size_kv    = len(self._knotvector)
# 		self._nb_ctrlpts = self._size_kv - (self._degree + 1)
# 		return

# 	def create_square(self, L):
# 		self._area    = L**2.0
# 		self._inertia = L**4.0/12.0
# 		self._Ks      = 5.0/6.0
# 		return

# 	def create_circle(self, R):
# 		self._area    = np.pi*R**2.0/4.0
# 		self._inertia = np.pi*R**4.0/4.0
# 		self._Ks      = 5.0/6.0
# 		return

# 	def create_rectangle(self, b, h):
# 		self._area    = b*h
# 		self._inertia = b*h**3.0/12.0
# 		self._Ks      = 5.0/6.0
# 		return

# 	def set_bspline_basis_weights(self):
# 		self._DB, self._DW = [], []
# 		qp_position, B0, B1, qp_weight = iga_find_basis_weights_opt(self._degree, self._knotvector)
# 		self._nb_qp  = len(qp_position)
# 		self._qp_cgg = qp_position
# 		self._DB.append(B0); self._DB.append(B1); self._DW = qp_weight
# 		return

# 	def compute_mechaprop(self):
# 		self._EA   = self._Emod*self._area*np.ones(self._nb_qp)
# 		self._EI   = self._Emod*self._inertia*np.ones(self._nb_qp)
# 		self._GAKs = self._Gmod*self._area*self._Ks*np.ones(self._nb_qp)
# 		return 

# 	def verify_mechaprop(self):
# 		item_list = [self._EA, self._EI, self._GAKs]
# 		if any([item is None for item in item_list]): self.compute_mechaprop()
# 		return

# 	def block_boundary(self, bound, values=[0, 0, 0], table=[True, True, True]):
# 		if bound in self._dod: raise Warning('Enter non-repeated boundary')
# 		if   bound == 0:
# 			for i in range(3): 
# 				if table[i]: self._dod.append(i*self._nb_ctrlpts)
# 		elif bound == -1:
# 			for i in range(3): 
# 				if table[i]: self._dod.append(self._nb_ctrlpts + i*self._nb_ctrlpts - 1)
# 		else: raise Warning('Only possible block first or last control point')
# 		for i in range(3): 
# 			if table[i]: self._gdod.append(values[i])
# 		return

# 	def find_free_ctrlpts(self): 
# 		set_dof = set([i for i in range(3*self._nb_ctrlpts)])
# 		set_dod = set(self._dod)
# 		diff    = set_dof.difference(set_dod)
# 		set_dof = list(diff)
# 		return set_dof

# 	def compute_timoshenko_stiffness(self, dw):
# 		" Compute Timoshenko stiffness matrix "
		
# 		EAdw  = self._EA * dw
# 		EAdw2 = EAdw * dw

# 		# Compute sub-matrices
# 		K11 = compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, self._EA)
# 		K12 = 0.5*compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, EAdw)
# 		# % K13 = 0
# 		K21 = compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, EAdw)
# 		K22 = (compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, self._GAKs) 
# 				+ 0.5*compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, EAdw2))
# 		K23 = compute_sub_advention_matrix(self._DB, self._DW, self._GAKs)
# 		# % K31 = 0
# 		K32 = K23.T
# 		K33 = (compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, self._EI) 
# 				+ compute_sub_mass_matrix(self._Jpq, self._DB, self._DW, self._GAKs))

# 		# Assemble matrix 
# 		nr   = self._nb_ctrlpts
# 		ind1 = range(nr)
# 		ind2 = range(nr, 2*nr)
# 		ind3 = range(2*nr, 3*nr)
# 		K = np.zeros((3*nr, 3*nr))
# 		K[np.ix_(ind1, ind1)] = K11
# 		K[np.ix_(ind1, ind2)] = K12
# 		K[np.ix_(ind2, ind1)] = K21
# 		K[np.ix_(ind2, ind2)] = K22
# 		K[np.ix_(ind2, ind3)] = K23
# 		K[np.ix_(ind3, ind2)] = K32
# 		K[np.ix_(ind3, ind3)] = K33

# 		return K

# 	def compute_timoshenko_tangent_stiffness(self, du, dw):
# 		" Compute Timoshenko tangent stiffness matrix "
		
# 		EAdw    = self._EA * dw
# 		EAdw2   = EAdw * dw
# 		EAdudw2 = EAdw2 + self._EA*du

# 		# Compute sub-matrices
# 		S11 = compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, self._EA)
# 		S12 = compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, EAdw)
# 		# % K13 = 0
# 		S21 = compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, EAdw)
# 		S22 = (compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, self._GAKs) 
# 				+ 0.5*compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, EAdw2)
# 				+ compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, EAdudw2)
# 		)
# 		S23 = compute_sub_advention_matrix(self._DB, self._DW, self._GAKs)
# 		# % K31 = 0
# 		S32 = S23.T
# 		S33 = (compute_sub_stiffness_matrix(self._Jpq, self._DB, self._DW, self._EI) 
# 				+ compute_sub_mass_matrix(self._Jpq, self._DB, self._DW, self._GAKs))

# 		# Assemble matrix 
# 		nr   = self._nb_ctrlpts
# 		ind1 = range(nr)
# 		ind2 = range(nr, 2*nr)
# 		ind3 = range(2*nr, 3*nr)
# 		S = np.zeros((3*nr, 3*nr))
# 		S[np.ix_(ind1, ind1)] = S11
# 		S[np.ix_(ind1, ind2)] = S12
# 		S[np.ix_(ind2, ind1)] = S21
# 		S[np.ix_(ind2, ind2)] = S22
# 		S[np.ix_(ind2, ind3)] = S23
# 		S[np.ix_(ind3, ind2)] = S32
# 		S[np.ix_(ind3, ind3)] = S33

# 		return S

# 	def compute_timoshenko_force(self, p, q):
# 		" Compute Timoshenko force vector "
# 		F1 = compute_sub_force_vector(self._Jpq, self._DB, self._DW, p)
# 		F2 = compute_sub_force_vector(self._Jpq, self._DB, self._DW, q)

# 		nr = self._nb_ctrlpts
# 		F  = np.zeros(3*nr)
# 		F[range(nr)]       = F1 
# 		F[range(nr, 2*nr)] = F2   
# 		return F

# 	def add_boundary_force(self, bound, values=[0, 0, 0]):
# 		dod = []
# 		if   bound == 0:
# 			for i in range(3): dod.append(i*self._nb_ctrlpts)
# 		elif bound == -1:
# 			for i in range(3): dod.append(self._nb_ctrlpts + i*self._nb_ctrlpts - 1)
# 		else: raise Warning('Only possible first or last control point')
# 		F = np.zeros(3*self._nb_ctrlpts)
# 		F[dod] = values
# 		return F

# 	def solve_timoshenko(self, Fext, tol=1e-9, nbIterNL=100):
# 		" Solves Timoshenko method "

# 		self.verify_mechaprop()
# 		if not len(self._dod): raise Warning('Add boundary conditions')
		
# 		nr  = self._nb_ctrlpts
# 		sol = np.zeros(3*nr)
# 		sol[self._dod] = self._gdod
# 		dof = self.find_free_ctrlpts()
		
# 		for i in range(nbIterNL):
# 			# Compute derivative of w
# 			dwdx = compute_derivative(self._Jpq, self._DB, sol[range(nr, 2*nr)])

# 			# Compute stiffness matrix
# 			K = self.compute_timoshenko_stiffness(dwdx)

# 			# Compute dF
# 			Fint    = K @ sol
# 			dF      = Fext[dof] - Fint[dof]
# 			errorNL = np.sqrt(np.dot(dF, dF))
# 			print('Iteration %s, error: %.5e' %(i, errorNL))
# 			if errorNL <= tol: break

# 			# Compute derivative of u
# 			dudx = compute_derivative(self._Jpq, self._DB, sol[range(nr)])

# 			# Compute stiffness
# 			S = self.compute_timoshenko_tangent_stiffness(dudx, dwdx)[np.ix_(dof, dof)]
# 			delta_sol = np.linalg.solve(S, dF)
# 			sol[dof] += delta_sol

# 		# Select important data
# 		u     = sol[range(nr)]
# 		w     = sol[range(nr, 2*nr)]
# 		theta = sol[range(2*nr, 3*nr)]

# 		return u, w, theta

# def compute_sub_stiffness_matrix(J, DB, W, prop):
# 	" Computes sub-stiffness matrix in Timoshenko theory "
# 	coefs = prop / J * W
# 	S     = DB[1] @ np.diag(coefs) @ DB[1].T
# 	return S

# def compute_sub_mass_matrix(J, DB, W, prop): 
# 	" Compute sub-mass matrix in Timoshenko theory "
# 	coefs = prop * J * W
# 	M     = DB[0] @ np.diag(coefs) @ DB[0].T
# 	return M

# def compute_sub_advention_matrix(DB, W, prop):
# 	" Compute sub-advention matrix in Timoshenko theory "
# 	coefs = prop * W
# 	A     = DB[1] @ np.diag(coefs) @ DB[0].T
# 	return A

# def compute_sub_force_vector(J, DB, W, prop):
# 	" Compute sub-force vector in Timoshenko theory "
# 	coefs = prop * J * W 
# 	F     = DB[0] @ coefs
# 	return F

# def compute_derivative(J, DB, u_ctrlpts):
# 	" Compute the derivative of u "
# 	der = DB[1].T @ u_ctrlpts / J
# 	return der
