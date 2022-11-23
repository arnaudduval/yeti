"""
.. This module helps to read information from a given geometry (YETI format)
.. and to create thermo-mechanical model (material, boundary conditions, etc.)
.. Joaquin Cornejo
"""

from .__init__ import *

# My libraries
from .D3viscoplasticity import *
from .base_functions import eval_basis_fortran
from .create_geomdl import geomdlModel

class thermoMechaModel(): 

	def __init__(self, modelIGA, material={}, Dirichlet={}, Neumann={}):
		
		print('\nInitializing thermo-mechanical model')

		print('Setting B-spline properties')
		self._sample_size = 101
		self._r_ = 2
		self._name = self.read_name(modelIGA)
		self._dim = self.read_dimensions(modelIGA)
		self._degree = self.read_degree(modelIGA)
		self._knotvector, self._size_kv = self.read_knotvector(modelIGA)
		self._ctrlpts = self.read_ControlPoints(modelIGA)
		self._set_parametric_properties()
		
		self._set_material(material)
		self._DirichletBound = np.zeros(self._nb_ctrlpts_total)
		self._set_dirichlet_boundaries(Dirichlet)
		self._set_neumann_condition(Neumann)

		return

	def _set_parametric_properties(self):
		" Sets B-spline properties "

		# Number of elements in each dimension
		self._nb_el = np.ones(3, dtype=int)  
		for dim in range(self._dim):
			self._nb_el[dim] = len(np.unique(self._knotvector[dim])) - 1
		
		# Whether or not Bspline respect IGA-WQ conditions 
		if any(self._nb_el[:self._dim]<2): 
			raise Warning('Model need at least 2 elements in each dimension')

		# Total number of elements
		self._nb_el_total = np.product(self._nb_el)

		# Set number of control points in each dimension
		self._nb_ctrlpts = np.ones(3, dtype=int)  
		for dim in range(self._dim):
			self._nb_ctrlpts[dim] = self._size_kv[dim] - self._degree[dim] - 1

		# Total number of control points  
		self._nb_ctrlpts_total = np.product(self._nb_ctrlpts)

		# Set number of quadrature points
		# In WQ approach
		self._nb_qp_wq = np.ones(3, dtype=int) 
		self._nb_qp_wq[:self._dim] = 2*(self._degree[:self._dim] 
										+ self._nb_el[:self._dim] + self._r_) - 5 #!! Eventually it can change
		# In IGA approach
		self._nb_qp_cgg = np.ones(3, dtype= int)
		self._nb_qp_cgg[:self._dim] = (self._degree[:self._dim] + 1)*self._nb_el[:self._dim]

		# Set total number of quadrature points
		self._nb_qp_cgg_total = np.prod(self._nb_qp_cgg)
		self._nb_qp_wq_total = np.prod(self._nb_qp_wq)

		print('B-spline properties :\n\
				-Dimensions: %d\n\
				-Degree: %s\n\
				-Number of elements: %s\n\
				-Number of control points : %s\n\
				-Number of WQ quadrature points: %s\n\
				-Number of Gauss quadrature points: %s\n\
				'
				%(
				self._dim, 
				self._degree,
				self._nb_el,
				self._nb_ctrlpts, 
				self._nb_qp_wq, 
				self._nb_qp_cgg
				)
			)

		return

	def _set_material(self, material:dict): 
		" Define material properties (thermal and mechanical) "

		# Set thermal properties
		self._conductivity = material.get('conductivity', None)
		self._capacity = material.get('capacity', None)

		# Set mechanical properties
		self._density = material.get('density', None)
		self._poissonCoef = material.get('poisson', None)
		self._youngModule = material.get('young', None)
		self._hardening = material.get('hardening', None)
		self._betaHard = material.get('betahard', None)
		self._sigmaY = material.get('sigmaY', None)

		# Initialize others properties
		self._Ctensor = None

		return
	
	def _set_extra_mechanical_properties(self):
		" Compute mechanical properties from E and nu"

		if self._youngModule is None or self._poissonCoef is None: pass
		else:
			# Create tensor
			identity = create_fourth_order_identity(self._dim)
			onekronone = create_one_kron_one(self._dim)

			# Get material properties
			E = self._youngModule
			nu = self._poissonCoef
			lamb = nu*E/((1+nu)*(1-2*nu))
			mu = E/(2*(1+nu))
			bulk = lamb + 2.0/3.0*mu

			# Create material tensor
			Idev = identity - 1.0/3.0*onekronone
			C = lamb*onekronone + 2*mu*identity
			S = 1.0/(9.0*bulk)*onekronone + 1.0/(2.0*mu)*Idev

			# Update
			self._Ctensor = C
			self._Stensor = S
			self._Idev = Idev
			self._lame_lambda = lamb
			self._lame_mu = mu
			self._lame_bulk = bulk

		return

	def _set_dirichlet_boundaries(self, Dirichlet:dict):
		" Gets free and blocked control points "
		
		# Thermal 
		try: 
			TTable = Dirichlet['thermal']
			Tdof, Tdod = self.get_boundaries_nodes(table=TTable)
			Tdof, Tdod = Tdof[0], Tdod[0]
		except: 
			TTable, Tdof, Tdod = None, None, None
		self._thermalDirichlet = TTable
		self._thermal_dof = Tdof
		self._thermal_dod = Tdod

		# Mechanical 
		try: 
			MTable = Dirichlet['mechanical']
			Mdof, Mdod = self.get_boundaries_nodes(table=MTable)
		except: 
			MTable, Mdof, Mdod = None, None, None
		self._mechanicalDirichlet = MTable
		self._mechanical_dof = Mdof
		self._mechanical_dod = Mdod

		return 

	def _add_thermal_IBC(self, table=None, temperature=0.0):
		"This function is first tentative of adding constant boundary conditions "
		if self._thermalDirichlet is None: raise Warning('Please define first total Dirichlet boundaries')
		if np.any(self._thermalDirichlet-table<0): raise Warning('Table is not well defined')
		if not np.any(table>0): raise Warning('At least one blocked face is needed')
		dod = self.get_boundaries_nodes(table=table)[1][0]
		self._DirichletBound[dod] += temperature*np.ones(len(dod))
		return 

	def _get_thermal_IBC(self): 
		if self._thermalDirichlet is None: raise Warning('Please define first total Dirichlet boundaries')
		return self._DirichletBound[self._thermal_dod]

	def _set_neumann_condition(self, Neumann:dict):
		" Gets Neumann control points and quadrature points"
		
		# Thermal 
		try: ForcesTable = Neumann['thermal']
		except: ForcesTable = None
		self._thermalNeumann = ForcesTable

		# Mechanical 
		try: ForcesTable = Neumann['mechanical']
		except: ForcesTable = None
		self._mechanicalNeumann = ForcesTable

		return 

	def _clear_material(self): 
		" Clears material "
		self._conductivity = None
		self._capacity = None
		self._poissonCoef = None
		self._youngModule = None
		self._hardening = None
		self._betaHard = None
		self._sigmaY = None

		return

	def _clear_Dirichlet(self):
		" Clears blocked boundaries "
		self._thermalDirichlet = None
		self._thermal_dof = None
		self._thermal_dod = None
		self._mechanicalDirichlet = None
		self._mechanical_dof = None
		self._mechanical_dod = None
		return

	def _clear_Neumann(self):
		" Clears Neumann boundaries "
		self._thermalNeumann = None
		self._mechanicalNeumann = None
		return

	def _verify_mechanics(self): 
		" Verifies if mechanical properties exits "
		proplist = [self._youngModule, self._hardening, self._betaHard, self._poissonCoef, self._sigmaY]
		if any([prop is None for prop in proplist]): raise Warning('Mechanics not well defined')
		if self._Ctensor is None: self._set_extra_mechanical_properties()
		return

	def _verify_thermal(self): 
		" Verifies if thermal properties exits "
		if self._conductivity is None: raise Warning('Conductivity not defined')
		if self._capacity is None: raise Warning('Capacity not defined')
		return
	
	def _compute_velocity(self, X, time_list, theta=1):
		"Computes velocity defined as dX/dt from X given and time"

		if np.shape(X)[1] != len(time_list): raise Warning('Different size')

		# Compute initial velocity
		if len(time_list) == 2:
			dt = time_list[1] - time_list[0]
			dX = 1.0/dt*(X[:, 1]-X[:, 0])
		elif len(time_list) > 2:
			dt = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt
			dX = (1.0/(dt*(factor - factor**2))*
				(X[:, 2] - (factor**2)*X[:, 1] - (1 - factor**2)*X[:, 0]))
		else:
			raise Warning("It is not possible")

		# Initialize V
		V = np.zeros(np.shape(X)); V[:, 0] = dX

		for i in range(1, len(time_list)):
			dt = time_list[i] - time_list[i-1]
			dX1 = 1.0/theta*(1.0/dt*(X[:, i] - X[:, i-1]) - (1 - theta)*dX)
			V[:, 0] = dX1; dX = dX1

		return V

	# =======================
	# READ FILE
	# =======================
	
	def read_name(self, modelIGA): 
		" Reads name from model "
		if isinstance(modelIGA, IGAparametrization): 
			try: name = modelIGA._name
			except: name = 'IGAparametrization'
		elif isinstance(modelIGA, geomdlModel): name = modelIGA._name
		return name

	def read_degree(self, modelIGA): 
		" Reads degree from model "
		if isinstance(modelIGA, IGAparametrization): degree = modelIGA._Jpqr.flatten()
		elif isinstance(modelIGA, geomdlModel): degree = modelIGA._degree
		if any(p == 1 for p in degree[:self._dim]): 
			raise Warning('Model must have at least degree p = 2')
		return degree

	def read_dimensions(self, modelIGA):
		" Reads dimensions from model "
		if isinstance(modelIGA, IGAparametrization): dim = modelIGA._dim[0]
		elif isinstance(modelIGA, geomdlModel): dim = modelIGA._dim
		# if dim != 3: raise Warning('Model must be 3D')
		if dim != 3: print("WARNING: Some functions have not been well-implemented for 2D geometries")
		return dim

	def read_knotvector(self, modelIGA):
		" Reads knot-vector from model "
		if isinstance(modelIGA, IGAparametrization): 
			knotvector = modelIGA._Ukv[0]
			size_kv = modelIGA._Nkv.flatten()
		elif isinstance(modelIGA, geomdlModel): 
			knotvector = modelIGA._knotvector
			size_kv = modelIGA._size_kv
		return knotvector, size_kv

	def read_ControlPoints(self, modelIGA): 
		" Reads control points from model "
		if isinstance(modelIGA, IGAparametrization): 
			ctrlpts = modelIGA._COORDS[:self._dim, :]
		elif isinstance(modelIGA, geomdlModel): 
			ctrlpts = modelIGA._ctrlpts[:self._dim, :]
		return ctrlpts

	# ========================
	# LOCAL FUNCTIONS
	# ========================

	def get_INC_table(self, nnz_dim): 
		""" Sets topology table, also known as INC: NURBS coordinates
		"""

		# Create INC: NURBS coordinates
		nnz_total = np.prod(nnz_dim)
		INC = np.zeros((nnz_total, 3), dtype= int)
		for i3 in range(nnz_dim[2]): 
			for i2 in range(nnz_dim[1]): 
				for i1 in range(nnz_dim[0]):
					genPos = i1 + i2*nnz_dim[0] + i3*nnz_dim[0]*nnz_dim[1]
					INC[genPos, :] = [i1, i2, i3]

		return INC

	def get_boundaries_nodes(self, table): 
		" Gets the indices of the blocked and free control points from table"

		# The table of dirichlet boundaries must be at least 3D
		table = np.atleast_3d(table)

		# Get number of degree of freedom (DOF) per node
		nbDOF = np.shape(table)[2]
		if np.shape(table)[0] < self._dim or np.shape(table)[1] != 2:
			raise Warning('Table is not well defined')

		# Find nodes
		nb_ctrlpts = self._nb_ctrlpts
		nb_ctrlpts_total = self._nb_ctrlpts_total
		INC = self.get_INC_table(nb_ctrlpts)
		dod_total = []; dof_total = []
		for i in range(nbDOF):
			dod = []
			for dim in range(self._dim):
				block_bound_dim = table[dim, :, i]
				
				if block_bound_dim[0]: 
					dod.extend(np.where(INC[:, dim] == 0)[0])

				if block_bound_dim[1]:
					dod.extend(np.where(INC[:, dim] == nb_ctrlpts[dim]-1)[0])

			# Rearrange
			dod = np.unique(dod)
			dof = set(np.arange(nb_ctrlpts_total, dtype= int)) - set(dod)
			dod_total.append(list(dod))
			dof_total.append(list(dof))

		return dof_total, dod_total

	def array2coo_matrix(self, data, indi, indj):
		" Computes coo sparse matrix  "

		nb_rows = max(indi) + 1
		nb_cols = max(indj) + 1
		sparse_matrix = sp.coo_matrix((data, (indi, indj)), shape=(nb_rows, nb_cols))
										
		return sparse_matrix

	def array2csr_matrix(self, data, indi, indj):
		" Computes csr sparse matrix "

		nb_rows = len(indi) - 1
		nb_cols = max(indj) + 1
		sparse_matrix = sp.csr_matrix((data, indj, indi), shape=(nb_rows, nb_cols))
										
		return sparse_matrix

	# ===========================
	# SPECIAL  FUNCTIONS 
	# ===========================

	def eval_jacobien_physicalPosition(self, dim, nnz, ctrlpts, DB): 
		""" Computes jacobien matrix and find the position in physical space 
			of points defined in parametric space. 
			This functions is deprecated. Please use Fortran for faster computations 
		"""
		print('Evaluating jacobien and physical position')
		start = time.process_time()
		
		J = np.zeros((dim, dim, nnz))
		PPS = np.zeros((dim, nnz))
		invJ = np.zeros((dim, dim, nnz))
		detJ = np.zeros(nnz)

		# Evaluate jacobien
		for j in range(dim):
			alpha = np.zeros(dim, dtype = int); alpha[j] = 1

			B = 1
			for d in range(dim):
				at = alpha[d] 
				B = sp.kron(DB[d][at], B)

			for i in range(dim):
				J[i, j, :] = sp.coo_matrix.dot(B.T, ctrlpts[i, :])
			
		# Evaluate position in physical space
		for i in range(dim):
			B = 1
			for d in range(dim):
				B = sp.kron(DB[d][0], B)

			PPS[i, :] = sp.coo_matrix.dot(B.T, ctrlpts[i, :])

		# Evaluate inverse and determinant
		for i in range(nnz):
			invJ[:, :, i] = np.linalg.inv(J[:, :, i])
			detJ[i] = np.linalg.det(J[:, :, i])
		
		stop = time.process_time()
		print('\t Time jacobien: %.5f s' %(stop-start))
		
		return J, detJ, invJ, PPS

	def eval_conductivity_coefficient(self, invJ, detJ, Kprop, isFortran=True): 
		" Computes conductivity coefficients "

		print('Getting conductivity coefficients')
		start = time.process_time()
		
		Kprop = np.atleast_3d(Kprop)
		nnz = len(detJ)
		coefs = np.zeros(np.shape(invJ))

		if isFortran:
			coefs, info = assembly.eval_conductivity_coefficient(invJ, detJ, Kprop)
			if info == 0: raise Warning('It is not possible to compute coefficients')
		else:
			if np.shape(Kprop)[2] == 1:
				for i in range(nnz): 
					inv_J = invJ[:,:,i]
					coefs[:, :, i] = inv_J @ Kprop[:, :, 0] @ inv_J.T * detJ[i]

			elif np.shape(Kprop)[2] == nnz:
				for i in range(nnz): 
					inv_J = invJ[:,:,i]
					coefs[:, :, i] = inv_J @ Kprop[:, :, i] @ inv_J.T * detJ[i]

			else: raise Warning('It is not possible to compute coefficients')

		stop = time.process_time()
		print('\tConductivity coefficients in : %.5f s' %(stop-start))

		return coefs

	def eval_capacity_coefficient(self, detJ, Cprop, isFortran=True): 
		" Computes capacity coefficients "

		print('Getting capacity coefficients')
		start = time.process_time()

		Cprop = np.atleast_1d(Cprop)
		nnz = np.size(detJ)
		coefs = np.zeros(nnz)

		if isFortran:
			coefs, info = assembly.eval_capacity_coefficient(detJ, Cprop)
			if info == 0: raise Warning('It is not possible to compute coefficients')
		else:
			if len(Cprop) == 1:
				for i in range(nnz): 
					coefs[i] = Cprop[0] * detJ[i]

			elif len(Cprop) == nnz:
				for i in range(nnz): 
					coefs[i] = Cprop[i] * detJ[i]

			else: raise Warning('It is not possible to compute coefficients')

		stop = time.process_time()
		print('\tCapacity coefficients in : %.5f s' %(stop-start))

		return coefs

	def eval_source_coefficient(self, fun, qp, detJ): 
		" Computes source coefficients "

		print('Getting source coefficients')
		start = time.process_time()

		qp = np.atleast_2d(qp)
		coefs = fun(qp)*detJ

		stop = time.process_time()
		print('\tSource coefficients in : %.5f s' %(stop-start))

		return coefs

	def eval_elastic_coefficient(self, invJ, detJ):
		""" Computes elasto-plastic coefficients.
			This function only consider linear isotropic case 
		"""

		# Set shape
		d = self._dim
		nnz = len(detJ)

		# Create tensors
		EE = create_incidence_matrix(d)
		coefs = np.zeros((d*d, d*d, nnz))
		
		# Create material tensor
		if self._Ctensor is None: raise Warning('C tensor not define')

		for k in range(nnz):
			invJt = invJ[:, :, k]
			detJt = detJ[k]
			for i in range(d):
				for j in range(d):
					ETCE = EE[:, :, i].T @ self._Ctensor @ EE[:, :, j]
					Dij = invJt @ ETCE @ invJt.T
					coefs[i*d:(i+1)*d, j*d:(j+1)*d, k] = Dij*detJt

		return coefs

	# ===========================
	# POST-PROCESSING 
	# ===========================

	def interpolate_field(self, samplesize=None, u_ctrlpts=None, nbDOF=3):
		" Interpolates the input field. It also returns the jacobien "

		# Get basis using fortran
		if samplesize == None: samplesize = self._sample_size
		knots = np.linspace(0, 1, samplesize)
		data, indices = [], []
		for dim in range(self._dim):  
			B, indi, indj = eval_basis_fortran(self._degree[dim], self._knotvector[dim], knots)
			data.append(B); indices.append(indi); indices.append(indj)

		# Get position and determinant 
		inputs = [*self._dim*[samplesize], *indices, *data, self._ctrlpts]
		if self._dim == 2:
			JJ_interp, detJJ_interp, _ = assembly.eval_jacobien_2d(*inputs)
			position_interp = assembly.interpolate_fieldphy_2d(*inputs)
		elif self._dim == 3: 
			JJ_interp, detJJ_interp, _ = assembly.eval_jacobien_3d(*inputs)
			position_interp = assembly.interpolate_fieldphy_3d(*inputs)

		# Get interpolation
		if u_ctrlpts is not None:
			u_temp = np.atleast_2d(u_ctrlpts)
			inputs = [*self._dim*[samplesize], *indices, *data, u_temp]

			if self._dim == 2: u_interp = assembly.interpolate_fieldphy_2d(*inputs)    
			elif self._dim == 3: u_interp = assembly.interpolate_fieldphy_3d(*inputs)
			if nbDOF == 1: u_interp = np.ravel(u_interp)
	
		else: u_interp = None

		return JJ_interp, position_interp, detJJ_interp, u_interp
	
	def export_results(self, u_ctrlpts=None, folder=None, nbDOF=3): 
		""" Export solution in VTK format. 
			It is possible to use Paraview to visualize data
		"""

		if folder == None: 
			full_path = os.path.realpath(__file__)
			dirname = os.path.dirname
			folder = dirname(dirname(full_path)) + '/results/'
		if not os.path.isdir(folder): 
			os.mkdir(folder)
		print("File saved in %s" %folder)

		if u_ctrlpts is None: pass
		elif isinstance(u_ctrlpts, np.ndarray): 
			if np.size(u_ctrlpts)%self._nb_ctrlpts_total != 0: 
				raise Warning('Not enough control points')
		else: raise Warning('Solution must be ndarray type')

		# ------------------
		# Get interpolation
		# ------------------
		_, qp_PS, detJ, u_interp = self.interpolate_field(u_ctrlpts=u_ctrlpts, nbDOF=nbDOF)
		mean_detJ = statistics.mean(detJ)
		detJ /= mean_detJ

		# ------------------
		# Export results
		# ------------------
		shape_pts = [1, 1, 1]
		for dim in range(self._dim): shape_pts[dim] = self._sample_size
		shape_pts = tuple(shape_pts)
		X1, X2, X3 = np.zeros(shape_pts), np.zeros(shape_pts), np.zeros(shape_pts)
		U = np.zeros((nbDOF, *shape_pts))
		DET = np.zeros(shape_pts)

		for k in range(shape_pts[2]):
			for j in range(shape_pts[1]):
				for i in range(shape_pts[0]):
					pos = i + j * self._sample_size + k * self._sample_size**2
					X1[i,j,k] = qp_PS[0, pos]
					X2[i,j,k] = qp_PS[1, pos]
					DET[i,j,k] = detJ[pos]
					if self._dim == 3: X3[i,j,k] = qp_PS[2, pos]
					if u_interp is not None: 
						u_interp = np.atleast_2d(u_interp)
						for l in range(nbDOF):
							U[l,i,j,k] = u_interp[l, pos]
		
		# Create point data 
		pointData = {}
		if u_interp is not None: 
			for l in range(nbDOF):
				varname = 'U' + str(l+1)
				pointData[varname] = U[l,:,:,:]
		pointData['detJ'] = DET

		# Export geometry
		name = folder + self._name
		gridToVTK(name, X1, X2, X3, pointData=pointData)
		
		return
