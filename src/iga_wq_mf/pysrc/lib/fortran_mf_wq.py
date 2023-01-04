"""
.. This module is an adaptation of Matrix Free (MF) 
.. Weighted Quadrature (WQ) for IsoGeometric Analysis (MF-WQ-IGA)
.. Joaquin Cornejo 
"""

from .__init__ import *

# My libraries
from .base_functions import (compute_eig_diag, 
							eigen_decomposition, 
							erase_rows_csr, 
							fast_diagonalization, 
							wq_find_basis_weights_fortran,
							compute_mean_3d
)
from .create_model import thermoMechaModel
from .D3viscoplasticity import *

class fortran_mf_wq(thermoMechaModel):

	def __init__(self, modelIGA, material={}, Dirichlet={}, Neumann={}):

		super().__init__(modelIGA, material=material, Dirichlet=Dirichlet, Neumann=Neumann)

		self._nb_qp, self._nb_qp_total = np.ones(self._dim, dtype=int), None
		self.eval_basis_weigths()
		self.eval_jacobien_physicalPosition()

		return

	def eval_basis_weigths(self): 
		" Computes basis and weights "

		print('Evaluating basis and weights')
		start = time.process_time()

		self._nnz_I, self._qp_dim, self._DB, self._DW, self._indices = [], [], [], [], []
		for dim in range(self._dim):  
			nnz_I, qp_position, basis, \
			weights, indi, indj = wq_find_basis_weights_fortran(self._degree[dim], self._knotvector[dim])
			self._nb_qp[dim] = len(qp_position)
			self._nnz_I.append(nnz_I); self._qp_dim.append(qp_position)
			self._DB.append(basis); self._DW.append(weights)
			self._indices.append(indi); self._indices.append(indj)

		# Update number of quadrature points
		self._nb_qp_total = np.prod(self._nb_qp)

		stop = time.process_time()
		print('\tBasis and weights in : %.5f s' %(stop-start))

		return

	def eval_jacobien_physicalPosition(self): 
		" Computes jacobien and physical position "

		print('Evaluating jacobien and physical position')
		start = time.process_time()
		
		inputs = [*self._nb_qp, *self._indices, *self._DB, self._ctrlpts]
		if self._dim == 2:
			self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_2d(*inputs)
			self._qp_PS = assembly.interpolate_fieldphy_2d(*inputs)
		if self._dim == 3:
			self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_3d(*inputs)
			self._qp_PS = assembly.interpolate_fieldphy_3d(*inputs)
		stop = time.process_time()
		print('\t Time jacobien: %.5f s' %(stop-start))

		return

	# ----------------------
	# WEIGHTED QUADRATURE
	# ----------------------

	def eval_source_coefficient(self, fun): 
		" Computes source coefficients "
		coefs = super().eval_source_coefficient(fun, self._qp_PS, self._detJ)
		return coefs

	def eval_bodyforce_coefficient(self, fun):
		" Computes body force coefficients "
		factor = self._detJ*self._density
		coefs = super().eval_source_coefficient(fun, self._qp_PS, factor)
		return coefs

	def eval_capacity_matrix(self, indi=None, indj=None): 
		" Computes capacity matrix "

		if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
		if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
		
		super()._verify_thermal()
		coefs = super().eval_capacity_coefficient(self._detJ, self._capacity)
		inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW, *self._nnz_I]

		start = time.process_time()
		if self._dim == 2: val_, indi_, indj_ = assembly.wq_get_capacity_2d(*inputs)
		if self._dim == 3: val_, indi_, indj_ = assembly.wq_get_capacity_3d(*inputs)
		matrix = super().array2csr_matrix(val_, indi_, indj_).tocsc()[indi, :][:, indj]
		stop = time.process_time()

		print('Capacity matrix assembled in : %.5f s' %(stop-start))
		
		return  matrix

	def eval_conductivity_matrix(self, indi= None, indj= None): 
		" Computes conductivity matrix "

		if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
		if indj is None: indj = np.arange(self._nb_ctrlpts_total, dtype=int)
		
		super()._verify_thermal()
		coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
		inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW, *self._nnz_I]

		start = time.process_time()
		if self._dim == 2: val_, indi_, indj_ = assembly.wq_get_conductivity_2d(*inputs)
		if self._dim == 3: val_, indi_, indj_ = assembly.wq_get_conductivity_3d(*inputs)
		matrix = super().array2csr_matrix(val_, indi_, indj_).tocsc()[indi, :][:, indj]
		stop = time.process_time()

		print('Conductivity matrix assembled in : %5f s' %(stop-start))
		
		return matrix

	def eval_stiffness_matrix(self, coefs=None): 
		" Computes conductivity matrix "
		
		if self._dim != 3: raise Warning('Not yet')
		super()._verify_mechanics()
		
		if coefs is None: coefs = super().eval_elastic_coefficient(self._invJ, self._detJ)
		inputs = [coefs, *self._nb_qp, *self._indices, *self._DB, *self._DW, *self._nnz_I]
		
		start = time.process_time()
		val_, indi_, indj_ = assembly.wq_get_stiffness_3d(*inputs)
		matrix = super().array2coo_matrix(val_, indi_, indj_).tocsc()
		stop = time.process_time()

		print('Stiffness matrix assembled in : %5f s' %(stop-start))
		
		return matrix

	def eval_Ku(self, u, table=None): 
		" Computes K u where K is conductivity matrix "

		super()._verify_thermal()
		coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
		inputs = self.get_input4MatrixFree(table=table)
	
		start = time.process_time()
		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: result = solver.mf_wq_get_ku_3d_py(coefs, *inputs, u)
		stop = time.process_time()
		timeCPU = stop - start

		return result, timeCPU

	def eval_Cu(self, u, table=None): 
		" Computes C u where C is capacity matrix "

		super()._verify_thermal()
		coefs = super().eval_capacity_coefficient(self._detJ, self._capacity)
		inputs = self.get_input4MatrixFree(table=table)
		
		start = time.process_time()
		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: result = solver.mf_wq_get_cu_3d_py(coefs, *inputs, u)
		stop = time.process_time()
		timeCPU = stop - start

		return result, timeCPU

	def eval_KCu(self, u, table=None, alpha=1.0, beta=1.0): 
		" Computes C u where C is capacity matrix "

		super()._verify_thermal()
		Kcoefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
		Ccoefs = super().eval_capacity_coefficient(self._detJ, self._capacity)
		inputs = self.get_input4MatrixFree(table=table)
		
		start = time.process_time()
		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: result = solver.mf_wq_get_kcu_3d_py(Ccoefs, Kcoefs, *inputs, u, alpha, beta)
		stop = time.process_time()
		timeCPU = stop - start

		return result, timeCPU

	def eval_Su(self, u):
		" Computes S u where S is stiffness matrix "
		
		if self._dim != 3: raise Warning('Until now not done')
		super()._verify_mechanics()
		prop = [self._youngModule, self._poissonCoef]
		inputs = [*self._nb_qp, *self._indices, *self._DB, *self._DW, self._invJ, self._detJ, prop]
		result = elastoplasticity.mf_wq_get_su_3d_py(*inputs, u)

		return result

	def eval_source_vector(self, fun, indi=None): 
		" Computes source vector "

		if indi is None: indi = np.arange(self._nb_ctrlpts_total, dtype=int)
		coefs = self.eval_source_coefficient(fun)
		inputs = [coefs, *self._nb_qp, *self._indices, *self._DW]

		start = time.process_time()
		if self._dim == 2: vector = assembly.wq_get_source_2d(*inputs)[indi]
		if self._dim == 3: vector = assembly.wq_get_source_3d(*inputs)[indi]
		stop = time.process_time()

		print('Source vector assembled in : %.5f s' %(stop-start))

		return vector

	def eval_force_surf(self):
		" Returns force vector at the surface. In 3D: surface integrals. "
		## For now, we only consider constants forces at the boundaries

		if self._dim != 3: raise Warning('Method only for 3D geometries')
		if self._mechanicalNeumann is None: raise Warning('Define Neumann condition')

		def get_info(nb):
			direction = int(np.floor(nb/2))
			if nb%2 == 1: side = 1
			else: side = 0
			return direction, side

		Ftemp = np.zeros((self._dim+1, self._nb_ctrlpts_total))

		# Get INC of control points and INC of quadrature points
		INC_CP = super().get_INC_table(self._nb_ctrlpts)
		INC_QP = super().get_INC_table(self._nb_qp)

		for dim in range(self._dim*2):
			direction, side = get_info(dim)
			force = self._mechanicalNeumann[dim]

			if np.array_equal(force, np.zeros(self._dim)): continue
			else:
				# Get control points and quadrature points list
				if side == 0: 
					CPList = np.where(INC_CP[:, direction] == 0)[0]
					QPList = np.where(INC_QP[:, direction] == 0)[0]

				elif side == 1: 
					CPList = np.where(INC_CP[:, direction] == self._nb_ctrlpts[direction]-1)[0]
					QPList = np.where(INC_QP[:, direction] == self._nb_qp[direction]-1)[0]
				
				CPList = list(np.sort(CPList))
				QPList = list(np.sort(QPList))

				# Modify Jacobien matrix
				valrange = [i for i in range(self._dim)]
				valrange.pop(direction)
				JJ = self._Jqp[:, :, QPList]
				JJ = JJ[:, valrange, :]

				# Compute surface force
				nnz, indices, data_W = [], [], []
				for _ in valrange:
					nnz.append(self._nb_qp[_]); data_W.append(self._DW[_])
					indices.append(self._indices[2*_]); indices.append(self._indices[2*_+1]) 
				
				FSurf = elastoplasticity.wq_get_forcesurf_3d(force, JJ, *nnz, *indices, *data_W)
				Ftemp[:-1, CPList] += FSurf
				Ftemp[-1, CPList] += 1

		# Update vector (average)
		F = Ftemp[:-1, :]
		# FList = np.where(Ftemp[-1, :] >= 2)[0]
		# F[:, FList] /= Ftemp[-1, FList]
		# To modify 

		return F

	def eval_force_body(self, fun):
		" Computes body (volumetric) force "

		if self._dim != 3: raise Warning('Method only for 3D geometries')
		coefs = self.eval_bodyforce_coefficient(fun)
		inputs = [coefs, *self._nb_qp, *self._indices, *self._DW]

		start = time.process_time()
		vector = elastoplasticity.wq_get_forcevol_3d(*inputs)
		stop = time.process_time()

		print('Body force vector assembled in : %.5f s' %(stop-start))

		return vector

	# ----------------------------------
	# HEAT TRANSFER (IN FORTRAN)
	# ----------------------------------   

	def get_input4MatrixFree(self, table=None):
		" Returns necessary inputs to compute the product between a matrix and a vector "
		
		if table is None: table = np.asarray([[0, 0], [0, 0], [0, 0]])
		indices, data_B, data_W = [], [], []
		for dim in range(self._dim):
			# Select data
			if np.array_equal(table[dim, :], [0, 0]): rows2erase = []
			if np.array_equal(table[dim, :], [0, 1]): rows2erase = [-1]
			if np.array_equal(table[dim, :], [1, 0]): rows2erase = [0]
			if np.array_equal(table[dim, :], [1, 1]): rows2erase = [0, -1]
			indi_t, indj_t, data_t = erase_rows_csr(rows2erase, 
									self._indices[2*dim], self._indices[2*dim+1],  
									[self._DB[dim], self._DW[dim]])
			
			# Extract data and append to list
			[dB, dW] = data_t
			indices.append(indi_t); indices.append(indj_t) 
			data_B.append(dB); data_W.append(dW)

		inputs = [*self._nb_qp, *indices, *data_B, *data_W]

		return inputs

	def MFsteadyHeat(self, b, nbIterPCG=100, threshold=1e-12, methodPCG='FDC'): 
		" Solves steady heat problems using directly substitution method "

		if self._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')

		super()._verify_thermal()
		coefs = super().eval_conductivity_coefficient(self._invJ, self._detJ, self._conductivity)
		inputs_tmp = self.get_input4MatrixFree(table=self._thermalDirichlet)
		inputs = [coefs, *inputs_tmp, b, nbIterPCG, threshold, methodPCG]

		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: sol, residue = solver.mf_wq_steady_heat_3d(*inputs)

		return sol, residue

	def MFtransientHeatNL(self, Fext, G=None, time_list=None, theta=1, 
						table_Kprop=None, table_Cprop=None, methodPCG='FDC', pos=-1):
		" Solves transient heat problem "

		if self._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
		if G is None or time_list is None: raise Warning('Important information missing')

		if table_Kprop is None: 
			table_Kprop = np.array([[0, 1]])
			print('WARNING: Default conductivity = 1.0')
		if table_Cprop is None: 
			table_Cprop = np.array([[0, 1]])
			print('WARNING: Default capacity = 1.0')

		dod = np.copy(self._thermal_dod); dod += 1
		inputs_tmp = self.get_input4MatrixFree(table=self._thermalDirichlet)
		inputs = [*inputs_tmp, *self._indices, *self._DB, *self._DW, time_list, Fext, dod, G, 
					table_Kprop, table_Cprop, self._invJ, self._detJ, theta, methodPCG]

		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: sol, resPCG, Kprop, Cprop = solver.mf_wq_transient_nonlinear_3d(*inputs)
		resPCG = resPCG[:, resPCG[0, :]>0]
		Cprop  = Cprop[:, :np.shape(resPCG)[1]]
		Kprop  = Kprop[:, :np.shape(resPCG)[1]]

		properties_first = []
		properties_first.append(self.interpolate_ControlPoints(datafield=Cprop[:, 0]))
		properties_first.append(self.interpolate_ControlPoints(datafield=Kprop[:, 0]))
		properties_first = np.array(properties_first)

		properties_last = []
		properties_last.append(self.interpolate_ControlPoints(datafield=Cprop[:, pos]))
		properties_last.append(self.interpolate_ControlPoints(datafield=Kprop[:, pos]))
		properties_last = np.array(properties_last)

		return sol, resPCG, properties_first, properties_last

	def interpolate_ControlPoints(self, funfield=None, datafield=None, nbIterPCG=100, threshold=1e-14):
		" Interpolation from parametric space to physical space "

		coefs = None
		if datafield is not None: coefs = datafield * self._detJ
		if funfield is not None: coefs = funfield(self._qp_PS) * self._detJ
		if coefs is None: raise Warning('Missing data')

		# Calculate vector
		inputs = [coefs, *self._nb_qp, *self._indices, *self._DW]
		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: vector = assembly.wq_get_source_3d(*inputs)

		# Solve linear system with fortran
		inputs = [self._detJ, *self._nb_qp, *self._indices, *self._DB, *self._DW, vector, nbIterPCG, threshold]
		start = time.process_time()
		u_interp, relres = solver.mf_wq_interpolate_cp_3d(*inputs)
		stop = time.process_time()
		res_end = relres[np.nonzero(relres)][-1]
		print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, res_end))

		return u_interp

	# ----------------------------------
	# ELASTO-PLASTICITY (IN FORTRAN)
	# ----------------------------------
	
	def MFelasticity_fortran(self, Fext=None, indi=None, nbIterPCG=100, threshold=1e-8, methodPCG='JMC'):
		" Solves a elasticity problem "
		
		if self._mechanicalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
		if indi is None or Fext is None: raise Warning('Impossible')
		super()._verify_mechanics()

		dod = deepcopy(indi)
		for i in range(len(dod)):
			dod_t = np.array(dod[i])
			dod_t += 1
			dod[i] = list(dod_t)

		prop = [self._youngModule, self._poissonCoef]
		inputs = [*self._nb_qp, *self._indices, *self._DB, *self._DW, Fext, *dod, self._mechanicalDirichlet, 
					self._invJ, self._detJ, prop, nbIterPCG, threshold, methodPCG]
		displacement, residue = elastoplasticity.mf_wq_elasticity_3d_py(*inputs)

		return displacement, residue
	
	def MFplasticity_fortran(self, Fext=None, indi=None):
		" Solves a plasticity problem "

		if self._dim != 3: raise Warning('Not yet')
		super()._verify_mechanics()
		if self._mechanicalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
		if indi is None or Fext is None: raise Warning('Impossible')

		dod = deepcopy(indi)
		for i in range(len(dod)):
			dod_t = np.array(dod[i])
			dod_t += 1
			dod[i] = list(dod_t)

		prop = np.array([self._youngModule, self._hardening, self._betaHard, self._poissonCoef, self._sigmaY])       
		inputs = [*self._nb_qp, *self._indices, *self._DB, *self._DW, Fext, *dod, 
					self._mechanicalDirichlet, self._invJ, self._detJ, prop]
		displacement, stress_vm = elastoplasticity.mf_wq_plasticity_3d(*inputs)

		return displacement, stress_vm

	# ----------------------------------
	# ELASTO-PLASTICITY (IN PYTHON)
	# ---------------------------------- 

	def compute_eigen_all(self, table=None, nddl=3, coefs=np.ones((3, 3))):
		""" Computes the eigen values and vectors considering Robin condition
			If Dirichlet condition is needed, then is better to try other function
		"""

		Deig     = np.zeros((3, self._nb_ctrlpts_total))
		eigvec_U = np.zeros((self._nb_ctrlpts[0], self._nb_ctrlpts[0], 3))
		eigvec_V = np.zeros((self._nb_ctrlpts[1], self._nb_ctrlpts[1], 3))
		eigvec_W = np.zeros((self._nb_ctrlpts[2], self._nb_ctrlpts[2], 3))

		for iddl in range(nddl):
			list_eig, list_vectors = [], []
			for dim in range(self._dim):
				indit = self._indices[2*dim]
				indjt = self._indices[2*dim+1]
				B = self._DB[dim]
				W = self._DW[dim]
				datat = [B[:, 0], B[:, 1], W[:, 0], W[:, -1]]
				robin = table[dim, :, iddl]
				eigt, Ut = eigen_decomposition(indit, indjt, datat, robin_condition=robin)
				list_eig.append(eigt)
				list_vectors.append(Ut)

			Deig[iddl, :] = compute_eig_diag(*list_eig, coefs=coefs[iddl, :])
			eigvec_U[:, :, iddl] = list_vectors[0]
			eigvec_V[:, :, iddl] = list_vectors[1]
			eigvec_W[:, :, iddl] = list_vectors[2]
		
		DU = [eigvec_U, eigvec_V, eigvec_W, Deig]

		return DU

	def compute_strain(self, u=None):
		" Compute strain field from displacement field "

		if self._dim != 3: raise Warning('Not yet')
		if u is None: raise Warning('Insert displacement')

		inputs = [*self._nb_qp, *self._indices, *self._DB, self._invJ, u, True]
		eps = elastoplasticity.interpolate_strain_3d(*inputs)

		return eps

	def compute_internal_force(self, coefs=None):
		"Compute internal force using sigma coefficients "

		if self._dim != 3: raise Warning('Not yet')
		if coefs is None: raise Warning('Insert coefficients')

		inputs = [coefs, *self._nb_qp, *self._indices, *self._DW]
		Fint = elastoplasticity.wq_get_forceint_3d_py(*inputs)

		return Fint

	def MFelasticity_py(self, coefs=None, DU=None, Fext=None, indi=None, 
							nbIterPCG=100, threshold=1e-8, isPrecond=True):
		" Solve linear system using Bi-CG Stab algorithm for elasticity problems "

		if self._dim != 3: raise Warning('Not yet')
		super()._verify_mechanics()
		if coefs is None: coefs = super().eval_elastic_coefficient(self._invJ, self._detJ)

		x = np.zeros(np.shape(Fext)); r = Fext
		clean_dirichlet_3d(r, indi)
		rhat, p = r, r
		rsold   = block_dot_product(self._dim, r, rhat)
		norm2b  = np.amax(np.absolute(r))

		if not isPrecond: # Without preconditioner 
			
			for i in range(nbIterPCG):
				Ap = self.eval_Su(p)
				clean_dirichlet_3d(Ap, indi)

				alpha = rsold/block_dot_product(self._dim, Ap, rhat)
				s = r -alpha*Ap

				As = self.eval_Su(s)
				clean_dirichlet_3d(As, indi)

				omega = block_dot_product(self._dim, As, s)/block_dot_product(self._dim, As, As)
				x += alpha*p + omega*s
				r = s - omega*As

				resPCG = np.amax(np.absolute(r))/norm2b
				print(resPCG)
				if resPCG <= threshold: break

				rsnew = block_dot_product(self._dim, r, rhat)
				beta = (alpha/omega)*(rsnew/rsold)
				p = r + beta*(p - omega*Ap)
				rsold = rsnew

		else: # With preconditioner

			if DU is None: 
				eigencoefs = np.ones((self._dim, self._dim))
				for i in range(3):
					su = compute_mean_3d(*self._nb_qp, coefs[i*self._dim, i*self._dim, :])
					sv = compute_mean_3d(*self._nb_qp, coefs[i*self._dim+1, i*self._dim+1, :])
					sw = compute_mean_3d(*self._nb_qp, coefs[i*self._dim+2, i*self._dim+2, :])
					eigencoefs[i, :] = [su, sv, sw] 
				DU = self.compute_eigen_all(table=self._mechanicalDirichlet, coefs=eigencoefs)
			U, V, W, D = DU[0], DU[1], DU[2], DU[3]

			for i in range(nbIterPCG):
				ptilde = fast_diagonalization(U, V, W, D, p, fdtype='elastic')
				clean_dirichlet_3d(ptilde, indi)

				Aptilde = self.eval_Su(ptilde)
				clean_dirichlet_3d(Aptilde, indi)

				alpha = rsold/block_dot_product(self._dim, Aptilde, rhat)
				s = r - alpha*Aptilde

				stilde = fast_diagonalization(U, V, W, D, s, fdtype='elastic')
				clean_dirichlet_3d(stilde, indi)

				Astilde = self.eval_Su(stilde)
				clean_dirichlet_3d(Astilde, indi)

				omega = block_dot_product(self._dim, Astilde, s)/block_dot_product(self._dim, Astilde, Astilde)
				x += alpha*ptilde + omega*stilde
				r = s - omega*Astilde

				resPCG = np.amax(np.absolute(r))/norm2b
				if resPCG <= threshold: break

				rsnew = block_dot_product(self._dim, r, rhat)
				beta = (alpha/omega)*(rsnew/rsold)
				p = r + beta*(p - omega*Aptilde)
				rsold = rsnew

		print('After %d iteration, the relative residue is %e' %(i, resPCG))

		return x

	def MFplasticity_py(self, Fext=None, indi=None, threshold=1e-8, nbIterNL=10, d=3):
		" Solves plasticity problem "

		if self._dim != 3: raise Warning('Only for 3D')
		super()._verify_mechanics()

		ddl = int(d*(d+1)/2)
		ep_n0 = np.zeros((ddl, self._nb_qp_total))
		a_n0  = np.zeros(self._nb_qp_total)
		b_n0  = np.zeros((ddl, self._nb_qp_total))
		ep_n1 = np.zeros((ddl, self._nb_qp_total))
		a_n1  = np.zeros(self._nb_qp_total)
		b_n1  = np.zeros((ddl, self._nb_qp_total))
		sigma = np.zeros((ddl, self._nb_qp_total))
		Cep   = np.zeros((ddl, ddl, self._nb_qp_total))
		disp  = np.zeros(np.shape(Fext))

		one      = create_second_order_identity(d)
		identity = create_fourth_order_identity(d)
		inputs   = [self._Ctensor, self._sigmaY, self._lame_bulk, self._lame_mu, 
					self._betaHard, self._hardening, identity, self._Idev, one]

		DU = self.compute_eigen_all(table=self._mechanicalDirichlet)

		for i in range(1, np.shape(Fext)[2]):

			ddisp = np.zeros(np.shape(disp[:, :, i-1]))
			Fstep = Fext[:, :, i]

			for j in range(nbIterNL): # Solver Newton-Raphson
				print('Step %d, iteration %d' %(i+1, j+1))

				# Compute strain as function of displacement
				d_n1 = disp[:, :, i-1] + ddisp
				eps = self.compute_strain(d_n1)
	
				# Closest point projection in perfect plasticity
				for k in range(self._nb_qp_total):
					sigmat, ep_n1t, a_n1t, b_n1t, Cept = cpp_combined_hardening(inputs, eps[:, k], ep_n0[:, k], a_n0[k], b_n0[:, k])
					sigma[:, k], ep_n1[:, k], a_n1[k], b_n1[:, k], Cep[:, :, k] = sigmat, ep_n1t, a_n1t, b_n1t, Cept

				# Compute coefficients to compute Fint and Stiffness
				coef_Fint, coef_Stiff = compute_plasticity_coef(sigma, Cep, self._invJ, self._detJ, d=d)

				# Compute Fint 
				Fint = self.compute_internal_force(coef_Fint)
				dF = Fstep - Fint
				clean_dirichlet_3d(dF, indi) 
				prod1 = block_dot_product(d, dF, dF)
				resNL = np.sqrt(prod1)
				print('Relative error: %.5e' %resNL)
				if resNL <= threshold: break
				
				delta_disp = self.MFelasticity_py(coefs=coef_Stiff, DU=DU, indi=indi, Fext=dF)
				ddisp += delta_disp 
		
			disp[:, :, i] = d_n1			
			ep_n0 = np.copy(ep_n1)
			a_n0 = np.copy(a_n1)
			b_n0 = np.copy(b_n1)

		return disp
