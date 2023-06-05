from lib.__init__ import *
from lib.lib_base import eraseRowsCSR, array2csr_matrix
from lib.lib_material import thermomat, mechamat, clean_dirichlet, block_dot_product
from lib.lib_part import part
from lib.lib_boundary import boundaryCondition, get_INCTable

class problem():
	def __init__(self, part:part, boundary:boundaryCondition, solverArgs:dict):
		self.material = None
		self.part     = part
		self.boundary = boundary
		self._nbIterPCG    = solverArgs.get('nbIterationsPCG', 100)
		self._nbIterNR     = solverArgs.get('nbIterationsNR', 20)
		self._thresholdPCG = solverArgs.get('PCGThreshold', 1e-12)
		self._thresholdNR  = solverArgs.get('NRThreshold', 1e-8)
		self._methodPCG    = solverArgs.get('PCGmethod', 'JMC')
		return

class heatproblem(problem):
	def __init__(self, material:thermomat, part:part, boundary:boundaryCondition, solverArgs=None):
		if solverArgs is None: solverArgs = dict()
		super().__init__(part, boundary, solverArgs)
		self.material = material
		return
	
	# Matrix free functions
	def get_input4MatrixFree(self, table=None):
		" Returns necessary inputs to compute the product between a matrix and a vector "
		
		if table is None: table = self.boundary.thDirichletTable
		indices, basis, weights = [], [], []
		for i in range(self.part.dim):
			# Select data
			if np.array_equal(table[i, :], [0, 0]): rows2erase = []
			if np.array_equal(table[i, :], [0, 1]): rows2erase = [-1]
			if np.array_equal(table[i, :], [1, 0]): rows2erase = [0]
			if np.array_equal(table[i, :], [1, 1]): rows2erase = [0, -1]
			indi_t, indj_t, data_t = eraseRowsCSR(rows2erase, 
									self.part.indices[2*i], self.part.indices[2*i+1],  
									[self.part.basis[i], self.part.weights[i]])
			
			# Extract data and append to list
			[basist, weightst] = data_t
			indices.append(indi_t); indices.append(indj_t) 
			basis.append(basist); weights.append(weightst)

		inputs = [*self.part.nbqp[:self.part.dim], *indices, *basis, *weights]

		return inputs
	
	def assemble_capacity(self, coefs=None, args=None):
		if coefs is None: coefs = self.material.eval_capacityCoefficients(self.part.detJ, args)
		nnz_I_list, nnz = np.array([-1, -1, -1], dtype=np.int32), 1
		inputs = [coefs, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, *self.part.weights]
		if self.part.dim == 2: 
			assembly.wq_get_capacity_2d(*inputs, nnz_I_list, nnz)
			nnz = np.prod(nnz_I_list)
			val, indi, indj = assembly.wq_get_capacity_2d(*inputs, nnz_I_list, nnz)
		if self.part.dim == 3: 
			assembly.wq_get_capacity_3d(*inputs, nnz_I_list, nnz)
			nnz = np.prod(nnz_I_list)
			val, indi, indj = assembly.wq_get_capacity_3d(*inputs, nnz_I_list, nnz)
		matrix = array2csr_matrix(val, indi, indj)
		return matrix
	
	def assemble_conductivity(self, coefs=None, args=None):
		if coefs is None: coefs = self.material.eval_conductivityCoefficients(self.part.detJ, args)
		nnz_I_list, nnz = np.array([-1, -1, -1], dtype=np.int32), 1
		inputs = [coefs, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, *self.part.weights]
		if self.part.dim == 2: 
			assembly.wq_get_conductivity_2d(*inputs, nnz_I_list, nnz)
			nnz = np.prod(nnz_I_list)
			val, indi, indj = assembly.wq_get_conductivity_2d(*inputs, nnz_I_list, nnz)
		if self.part.dim == 3: 
			assembly.wq_get_conductivity_3d(*inputs, nnz_I_list, nnz)
			nnz = np.prod(nnz_I_list)
			val, indi, indj = assembly.wq_get_conductivity_3d(*inputs, nnz_I_list, nnz)
		matrix = array2csr_matrix(val, indi, indj)
		return matrix

	def eval_mfConductivity(self, u, coefs=None, table=None, args=None):
		inputs = self.get_input4MatrixFree(table=table)
		if coefs is None: coefs  = self.material.eval_conductivityCoefficients(self.part.invJ, self.part.detJ, args)
		if self.part.dim == 2: raise Warning('Until now not done')
		if self.part.dim == 3: result = heatsolver.mf_wq_get_ku_3d(coefs, *inputs, u)
		return result
	
	def eval_mfCapacity(self, u, coefs=None, table=None, args=None): 
		inputs = self.get_input4MatrixFree(table=table)
		if coefs is None: coefs  = self.material.eval_capacityCoefficients(self.part.detJ, args)
		if self.part.dim == 2: raise Warning('Until now not done')
		if self.part.dim == 3: result = heatsolver.mf_wq_get_cu_3d(coefs, *inputs, u)
		return result

	def eval_volForce(self, fun, indi=None): 
		if indi is None: indi = np.arange(self.part.nbctrlpts_total, dtype=int)
		coefs = self.material.eval_heatForceCoefficients(fun, self.part.detJ, self.part.qpPhy)
		inputs = [coefs, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights]
		if self.part.dim == 2: vector = heatsolver.wq_get_heatvol_2d(*inputs)[indi]
		if self.part.dim == 3: vector = heatsolver.wq_get_heatvol_3d(*inputs)[indi]
		return vector
	
	def eval_surfForce(self, fun, nbFacePosition):
		if self.part.dim != 3: raise Warning('Method only for 3D geometries')

		def get_faceInfo(nb):
			direction = int(np.floor(nb/2))
			if nb%2 == 1: side = 1
			else: side = 0
			return direction, side

		vector = np.zeros(self.part.nbctrlpts_total)
		INC_ctrlpts = get_INCTable(self.part.nbctrlpts)
		INC_quadpts = get_INCTable(self.part.nbqp)
		direction, side = get_faceInfo(nbFacePosition)

		# Get control points and quadrature points list
		if side == 0: 
			CPList = np.where(INC_ctrlpts[:, direction] == 0)[0]
			QPList = np.where(INC_quadpts[:, direction] == 0)[0]

		elif side == 1: 
			CPList = np.where(INC_ctrlpts[:, direction] == self.part.nbctrlpts[direction]-1)[0]
			QPList = np.where(INC_quadpts[:, direction] == self.part.nbqp[direction]-1)[0]
		
		CPList = list(np.sort(CPList))
		QPList = list(np.sort(QPList))

		# Modify Jacobien matrix
		valrange = [i for i in range(self.part.dim)]
		valrange.pop(direction)
		JJ = self.part.Jqp[:, :, QPList]
		JJ = JJ[:, valrange, :]

		# Get force values at quadrature points
		qpPhy = self.part.qpPhy[:, QPList]
		coefs = fun(qpPhy)

		# Compute surface force
		nnz, indices, weights = [], [], []
		for _ in valrange:
			nnz.append(self.part.nbqp[_]); weights.append(self.part.weights[_])
			indices.append(self.part.indices[2*_]); indices.append(self.part.indices[2*_+1]) 
		
		tmp = plasticitysolver.wq_get_heatsurf_3d(coefs, JJ, *nnz, *indices, *weights)
		vector[:, CPList] = tmp

		return vector

	def interpolate_temperature(self, uctrlpts):
		basis   = self.part.basis
		indices = self.part.indices
		nbqp    = self.part.nbqp[:self.part.dim]

		# Get position and determinant 
		inputs = [*nbqp, *indices, *basis, np.atleast_2d(uctrlpts)]
		if self.part.dim == 2:
			uinterp = geophy.interpolate_fieldphy_2d(*inputs)
		elif self.part.dim == 3: 
			uinterp = geophy.interpolate_fieldphy_3d(*inputs)
		uinterp = np.ravel(uinterp)
		return uinterp

	# Solve using fortran
	def solveInterpolationProblemFT(self, funfield=None, datafield=None, nbIterPCG=None):
		coefs = None
		if datafield is not None: coefs = datafield*self.part.detJ
		if funfield is not None:  coefs = funfield(self.part.qpPhy)*self.part.detJ
		if coefs is None: raise Warning('Missing data')

		# Calculate vector
		inputs = [coefs, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights]
		if self.part.dim == 2: raise Warning('Until now not done')
		if self.part.dim == 3: vector = heatsolver.wq_get_heatvol_3d(*inputs)

		# Solve linear system with fortran
		if nbIterPCG is None: nbIterPCG = self._nbIterPCG
		inputs = [self.part.detJ, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, 
	    		 *self.part.weights, vector, nbIterPCG, self._thresholdPCG]
		start = time.process_time()
		u_interp, relres = heatsolver.mf_wq_interpolate_cp_3d(*inputs)
		stop = time.process_time()
		res_end = relres[np.nonzero(relres)][-1]
		print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, res_end))
		return u_interp

	def solveSteadyHeatProblemFT(self, b, coefs=None, nbIterPCG=None, methodPCG=None):
		if coefs is None: coefs  = self.material.eval_conductivityCoefficients(self.part.invJ, 
															self.part.detJ, self.part.qpPhy)
		tmp = self.get_input4MatrixFree(table=self.boundary.thDirichletTable)
		if nbIterPCG is None: nbIterPCG = self._nbIterPCG
		if methodPCG is None: methodPCG = self._methodPCG
		inputs = [coefs, *tmp, b, nbIterPCG, self._thresholdPCG, methodPCG]

		if self.part.dim == 2: raise Warning('Until now not done')
		if self.part.dim == 3: sol, residue = heatsolver.mf_wq_steady_heat_3d(*inputs)

		return sol, residue
	
	def solveLinearTransientHeatProblemFT(self, dt, b, theta=1.0, Ccoefs=None, Kcoefs=None, 
										nbIterPCG=None, methodPCG=None, args={}):
		if Ccoefs is None: 
			temperature = args.get('temperature')
			inpt = self.part.qpPhy
			inpt = np.row_stack((inpt, temperature*np.ones(np.size(inpt, axis=1))))
			Ccoefs  = self.material.eval_capacityCoefficients(self.part.detJ, inpt)
			
		if Kcoefs is None:
			temperature = args.get('temperature')
			inpt = self.part.qpPhy
			inpt = np.row_stack((inpt, temperature*np.ones(np.size(inpt, axis=1))))
			Kcoefs  = self.material.eval_conductivityCoefficients(self.part.invJ, self.part.detJ, inpt)
		
		tmp = self.get_input4MatrixFree(table=self.boundary.thDirichletTable)
		if nbIterPCG is None: nbIterPCG = self._nbIterPCG
		if methodPCG is None: methodPCG = self._methodPCG
		inputs = [Ccoefs, Kcoefs, *tmp, b, theta*dt, nbIterPCG, self._thresholdPCG, methodPCG]

		if self.part.dim == 2: raise Warning('Until now not done')
		if self.part.dim == 3: sol, residue = heatsolver.mf_wq_lineartransient_heat_3d(*inputs)

		return sol, residue

	# Solve using python
	def solveNLTransientHeatProblemPy(self, Tinout, time_list, Fext, theta=1.0, thresholdNR=None):
		if thresholdNR is None: thresholdNR = self._thresholdNR
		m, n = np.shape(Tinout)
		nbSteps         = len(time_list)
		nbctrlpts_total = self.part.nbctrlpts_total
		if n != nbSteps: raise Warning('Not possible')
		if m != nbctrlpts_total: raise Warning('Not possible')
		dod, _, dof = self.boundary.getThermalBoundaryConditionInfo()

		VVn0 = np.zeros(nbctrlpts_total)
		resPCG_list = []
		if nbSteps == 2:
			dt = time_list[1] - time_list[0]
			VVn0[dod] = 1.0/dt*(Tinout[dod, 1] - Tinout[dod, 0])
		elif nbSteps > 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			VVn0[dod] = 1.0/(dt1*(factor - factor**2))*(Tinout[dod, 2] 
					- (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
		else:
			raise Warning('At least 2 steps')
		
		for i in range(1, nbSteps):
			# Get delta time
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			TTn0 = np.copy(Tinout[:, i-1])

			# Get approximative values of new step
			TTn1 = TTn0 + dt*(1-theta)*VVn0; TTn1[dod] = np.copy(Tinout[dod, i])
			TTn10 = np.copy(TTn1); VVn1 = np.zeros(np.shape(VVn0))
			VVn1[dod] = 1.0/theta*(1.0/dt*(Tinout[dod, i]-Tinout[dod, i-1]) - (1-theta)*VVn0[dod])
			Fstep = Fext[:, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR):

				# Compute temperature and properties at each quadrature point
				TTinterp = self.interpolate_temperature(TTn1)
			
				# Compute internal force
				inpt = self.part.qpPhy
				inpt = np.row_stack((inpt, TTinterp*np.ones(np.size(inpt, axis=1))))
				Ccoefs  = self.material.eval_capacityCoefficients(self.part.detJ, inpt)
				Kcoefs  = self.material.eval_conductivityCoefficients(self.part.invJ, self.part.detJ, inpt)

				CdTemp = self.eval_mfCapacity(VVn1, coefs=Ccoefs, table=np.zeros((3, 2), dtype=bool))
				KTemp  = self.eval_mfConductivity(TTn1, coefs=Kcoefs, table=np.zeros((3, 2), dtype=bool))
				Fint   = KTemp + CdTemp

				# Compute residue
				ddFF    = Fstep - Fint
				ddFFdof = ddFF[dof]
				resNL   = np.sqrt(np.dot(ddFFdof, ddFFdof))
				print('NR error: %.5e' %resNL)
				if resNL <= thresholdNR: break

				# Iterative solver
				resPCG = np.array([i, j+1])
				ddVV, resPCGt = self.solveLinearTransientHeatProblemFT(dt, ddFFdof, Ccoefs=Ccoefs, Kcoefs=Kcoefs, theta=theta)
				resPCG = np.append(resPCG, resPCGt)
				resPCG_list.append(resPCG)

				# Update values
				VVn1[dof] += ddVV
				TTn1[dof] = TTn10[dof] + theta*dt*VVn1[dof]

			Tinout[:, i] = np.copy(TTn1)
			VVn0 = np.copy(VVn1)

		return resPCG_list

class mechaproblem(problem):
	def __init__(self, material:mechamat, part:part, boundary:boundaryCondition, solverArgs=None):
		if solverArgs is None: solverArgs = dict()
		super().__init__(part, boundary, solverArgs)
		self.material = material
		return
	
	# Matrix free functions
	def eval_mfStiffness(self, displacement, tensorArgs=None):
		if self.part.dim != 3: raise Warning('Until now not done')
		if tensorArgs is None: 
			dimen  = self.part.dim
			nvoigt = int(dimen*(dimen+1)/2)
			tensorArgs = np.zeros((nvoigt+3, self.part.nbqp_total))
			tensorArgs[0, :] = self.material.lame_lambda
			tensorArgs[1, :] = self.material.lame_mu
		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, 
	    			*self.part.basis, *self.part.weights, 
					self.part.invJ, self.part.detJ, tensorArgs]
		if   self.part.dim == 2:  result = plasticitysolver.mf_wq_get_su_2d(*inputs, displacement)
		elif self.part.dim == 3:  result = plasticitysolver.mf_wq_get_su_3d(*inputs, displacement)
		return result
	
	def eval_volForce(self, fun):
		coefs = self.material.eval_volForceCoefficients(fun, self.part.detJ, self.part.qpPhy)
		inputs = [coefs, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights]
		if   self.part.dim == 2: vector = plasticitysolver.wq_get_forcevol_2d(*inputs)
		elif self.part.dim == 3: vector = plasticitysolver.wq_get_forcevol_3d(*inputs)
		return vector
	
	def eval_surfForce(self, fun, nbFacePosition):
		" Returns force vector at the surface. In 3D: surface integrals. "

		def get_faceInfo(nb):
			direction = int(np.floor(nb/2))
			if nb%2 == 1: side = 1
			else: side = 0
			return direction, side

		vector = np.zeros((self.part.dim, self.part.nbctrlpts_total))
		INC_ctrlpts = get_INCTable(self.part.nbctrlpts)
		INC_quadpts = get_INCTable(self.part.nbqp)
		direction, side = get_faceInfo(nbFacePosition)
		if direction>=2*self.part.dim: raise Warning('Not possible')

		# Get control points and quadrature points list
		if side == 0: 
			CPList = np.where(INC_ctrlpts[:, direction] == 0)[0]
			QPList = np.where(INC_quadpts[:, direction] == 0)[0]

		elif side == 1: 
			CPList = np.where(INC_ctrlpts[:, direction] == self.part.nbctrlpts[direction]-1)[0]
			QPList = np.where(INC_quadpts[:, direction] == self.part.nbqp[direction]-1)[0]
		
		CPList = list(np.sort(CPList))
		QPList = list(np.sort(QPList))

		# Modify Jacobien matrix
		valrange = [i for i in range(self.part.dim)]
		valrange.pop(direction)
		JJ = self.part.Jqp[:, :, QPList]
		JJ = JJ[:, valrange, :]

		# Get force values at quadrature points
		qpPhy = self.part.qpPhy[:, QPList]
		coefs = fun(qpPhy)

		# Compute surface force
		nnz, indices, weights = [], [], []
		for _ in valrange:
			nnz.append(self.part.nbqp[_]); weights.append(self.part.weights[_])
			indices.append(self.part.indices[2*_]); indices.append(self.part.indices[2*_+1]) 
		
		if   self.part.dim == 2: tmp = plasticitysolver.wq_get_forcesurf_2d(coefs, JJ, *nnz, *indices, *weights)
		elif self.part.dim == 3: tmp = plasticitysolver.wq_get_forcesurf_3d(coefs, JJ, *nnz, *indices, *weights)
		vector[:, CPList] = tmp

		return vector

	# Solve using fortran
	def solveElasticityProblemFT(self, Fext, tensorArgs=None, nbIterPCG=None, methodPCG=None):
		dod_total = deepcopy(self.boundary.mchdod)
		for i, dod in enumerate(dod_total):
			tmp = dod + 1; dod_total[i] = tmp

		dimen  = self.part.dim
		nvoigt = int(dimen*(dimen+1)/2)
		prop = [self.material.elasticmodulus, self.material.poissonratio, self.material.elasticlimit]
		if tensorArgs is None:
			tensorArgs = np.zeros((nvoigt+3, self.part.nbqp_total))
			tensorArgs[0, :] = self.material.lame_lambda
			tensorArgs[1, :] = self.material.lame_mu
		if nbIterPCG is None: nbIterPCG = self._nbIterPCG
		if methodPCG is None: methodPCG = self._methodPCG
		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, 
				*self.part.weights, Fext, *dod_total, self.boundary.mchDirichletTable, 
				self.part.invJ, self.part.detJ, prop, tensorArgs, nbIterPCG, self._thresholdPCG, methodPCG]
		if   self.part.dim == 2: displacement, residue = plasticitysolver.mf_wq_elasticity_2d(*inputs)
		elif self.part.dim == 3: displacement, residue = plasticitysolver.mf_wq_elasticity_3d(*inputs)

		return displacement, residue

	# Solve using python
	def compute_strain(self, displacement, isVoigt=False):
		" Compute strain field from displacement field "

		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, self.part.invJ, displacement, isVoigt]
		if   self.part.dim == 2: eps = plasticitysolver.interpolate_strain_2d(*inputs)
		elif self.part.dim == 3: eps = plasticitysolver.interpolate_strain_3d(*inputs)

		return eps
	
	def compute_intForce(self, stress):
		"Compute internal force using sigma coefficients "

		inputs = [stress, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.invJ, self.part.detJ]
		if   self.part.dim == 2: Fint = plasticitysolver.wq_get_intforce_2d(*inputs)
		elif self.part.dim == 3: Fint = plasticitysolver.wq_get_intforce_3d(*inputs)

		return Fint
	
	def solvePlasticityProblemPy(self, Fext, nbIterPCG=None, methodPCG=None, thresholdNR=None): 
		if thresholdNR is None: thresholdNR = self._thresholdNR
		if not self.material._isPlasticityPossible: raise Warning('Plasticity not defined')
		if nbIterPCG is None: nbIterPCG = self._nbIterPCG
		if methodPCG is None: methodPCG = self._methodPCG

		d     = self.part.dim
		ddl   = int(d*(d+1)/2)

		nbqp_total = self.part.nbqp_total
		pls_n0 = np.zeros((ddl, nbqp_total))
		a_n0  = np.zeros(nbqp_total)
		b_n0  = np.zeros((ddl, nbqp_total))
		pls_n1 = np.zeros((ddl, nbqp_total))
		a_n1  = np.zeros(nbqp_total)
		b_n1  = np.zeros((ddl, nbqp_total))
		stress = np.zeros((ddl, nbqp_total))
		Cep   = np.zeros((ddl+3, nbqp_total))
		disp  = np.zeros(np.shape(Fext))
		resPCG_list = []

		for i in range(1, np.shape(Fext)[2]):

			ddisp = np.zeros(np.shape(disp[:, :, i-1]))
			Fstep = Fext[:, :, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR): # Solver Newton-Raphson

				# Compute strain as function of displacement
				d_n1 = disp[:, :, i-1] + ddisp
				strain = self.compute_strain(d_n1)
	
				# Closest point projection in perfect plasticity
				for k in range(nbqp_total):
					sigmat, pls_n1t, a_n1t, b_n1t, Cept = self.material.returnMappingAlgorithm(strain[:, k], pls_n0[:, k], a_n0[k], b_n0[:, k])
					stress[:, k], pls_n1[:, k], a_n1[k], b_n1[:, k], Cep[:, k] = sigmat, pls_n1t, a_n1t, b_n1t, Cept

				# Compute Fint 
				Fint = self.compute_intForce(stress)
				dF   = Fstep - Fint
				clean_dirichlet(dF, self.boundary.mchdod) 
				prod1 = block_dot_product(d, dF, dF)
				resNL = np.sqrt(prod1)
				print('NR error: %.5e' %resNL)
				if resNL <= thresholdNR: break
				
				resPCG = np.array([i, j+1])
				vtmp, resPCGt = self.solveElasticityProblemFT(Fext=dF, tensorArgs=Cep, nbIterPCG=nbIterPCG, methodPCG=methodPCG)
				resPCG = np.append(resPCG, resPCGt)
				resPCG_list.append(resPCG)

				ddisp  += vtmp
		
			disp[:, :, i] = d_n1			
			pls_n0 = np.copy(pls_n1)
			a_n0 = np.copy(a_n1)
			b_n0 = np.copy(b_n1)
		return disp, resPCG_list

