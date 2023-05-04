from lib.__init__ import *
from lib.lib_base import eraseRowsCSR
from lib.lib_material import *
from lib.lib_model import *
from lib.lib_step import *

class heatproblem():

	def __init__(self, material: thermomat, model: part, boundary: step, **kwargs):
		self._material = material
		self._model    = model
		self._boundary = boundary
		self._kwargs   = kwargs
		self.__extractInfo()
		return
	
	def __extractInfo(self):
		kwargs = self._kwargs
		self._nbIterPCG    = kwargs.get('nbIterationsPCG', 100)
		self._nbIterNR     = kwargs.get('nbIterationsNR', 20)
		self._thresholdPCG = kwargs.get('PCGThreshold', 1e-12)
		self._thresholdNR  = kwargs.get('NRThreshold', 1e-6)
		self._methodPCG    = kwargs.get('PCGmethod', 'FDC')
		return
	
	# Matrix free functions
	def get_input4MatrixFree(self, table=None):
		" Returns necessary inputs to compute the product between a matrix and a vector "
		
		if table is None: table = self._boundary._thDirichletTable
		indices, basis, weights = [], [], []
		for i in range(self._model._dim):
			# Select data
			if np.array_equal(table[i, :], [0, 0]): rows2erase = []
			if np.array_equal(table[i, :], [0, 1]): rows2erase = [-1]
			if np.array_equal(table[i, :], [1, 0]): rows2erase = [0]
			if np.array_equal(table[i, :], [1, 1]): rows2erase = [0, -1]
			indi_t, indj_t, data_t = eraseRowsCSR(rows2erase, 
									self._model._indices[2*i], self._model._indices[2*i+1],  
									[self._model._basis[i], self._model._weights[i]])
			
			# Extract data and append to list
			[basist, weightst] = data_t
			indices.append(indi_t); indices.append(indj_t) 
			basis.append(basist); weights.append(weightst)

		inputs = [*self._model._nbqp, *indices, *basis, *weights]

		return inputs

	def eval_mfConductivity(self, u, coefs=None, table=None, **kwargs):

		inputs = self.get_input4MatrixFree(table=table)
		if coefs is None: 
			inpt   = kwargs.get(['input'])
			coefs  = self._material.eval_conductivityCoefficients(self._model._invJ, self._model._detJ, inpt)
		if self._model._dim == 2: raise Warning('Until now not done')
		if self._model._dim == 3: result = heatsolver.mf_wq_get_ku_3d(coefs, *inputs, u)

		return result
	
	def eval_mfCapacity(self, u, coefs=None, table=None, **kwargs): 

		inputs = self.get_input4MatrixFree(table=table)
		if coefs is None: 
			inpt   = kwargs.get(['input'])
			coefs  = self._material.eval_capacityCoefficients(self._model._detJ, inpt)
		if self._model._dim == 2: raise Warning('Until now not done')
		if self._model._dim == 3: result = heatsolver.mf_wq_get_cu_3d(coefs, *inputs, u)

		return result

	def eval_bodyForce(self, fun, indi=None): 

		if indi is None: indi = np.arange(self._model._nbctrlpts_total, dtype=int)
		coefs = self._material.eval_heatForceCoefficients(fun, self._model._detJ, self._model._qpPhy)
		inputs = [coefs, *self._model._nbqp, *self._model._indices, *self._model._weights]
		if self._model._dim == 2: raise Warning('Not done yet')
		if self._model._dim == 3: vector = heatsolver.wq_get_bodyheat_3d(*inputs)[indi]

		return vector
	
	def eval_surfForce(self, fun, nbFacePosition):
		if self._model._dim != 3: raise Warning('Method only for 3D geometries')

		def get_faceInfo(nb):
			direction = int(np.floor(nb/2))
			if nb%2 == 1: side = 1
			else: side = 0
			return direction, side

		vector = np.zeros(self._model._nbctrlpts_total)
		INC_ctrlpts = self._boundary.__get_INCTable(self._model._nbctrlpts)
		INC_quadpts = self._boundary.__get_INCTable(self._model._nbqp)
		direction, side = get_faceInfo(nbFacePosition)

		# Get control points and quadrature points list
		if side == 0: 
			CPList = np.where(INC_ctrlpts[:, direction] == 0)[0]
			QPList = np.where(INC_quadpts[:, direction] == 0)[0]

		elif side == 1: 
			CPList = np.where(INC_ctrlpts[:, direction] == self._model._nbctrlpts[direction]-1)[0]
			QPList = np.where(INC_quadpts[:, direction] == self._model._nbqp[direction]-1)[0]
		
		CPList = list(np.sort(CPList))
		QPList = list(np.sort(QPList))

		# Modify Jacobien matrix
		valrange = [i for i in range(self._model._dim)]
		valrange.pop(direction)
		JJ = self._model._Jqp[:, :, QPList]
		JJ = JJ[:, valrange, :]

		# Get force values at quadrature points
		qpPhy = self._model._qpPhy[:, QPList]
		coefs = fun(qpPhy)

		# Compute surface force
		nnz, indices, weights = [], [], []
		for _ in valrange:
			nnz.append(self._model._nbqp[_]); weights.append(self._model._weights[_])
			indices.append(self._model._indices[2*_]); indices.append(self._model._indices[2*_+1]) 
		
		tmp = plasticitysolver.wq_get_heatflux_3d(coefs, JJ, *nnz, *indices, *weights)
		vector[:, CPList] = tmp

		return vector

	def interpolateTemperature(self, uctrlpts):
		basis   = self._model._basis
		indices = self._model._indices
		nbqp    = self._model._nbqp

		# Get position and determinant 
		inputs = [*nbqp, *indices, *basis, np.atleast_2d(uctrlpts)]
		if self._model._dim == 2:
			uinterp = interpolation.interpolate_fieldphy_2d(*inputs)
		elif self._model._dim == 3: 
			uinterp = interpolation.interpolate_fieldphy_3d(*inputs)
		uinterp = np.ravel(uinterp)
		return uinterp

	# Solve using fortran
	def solveInterpolationProblemFT(self, funfield=None, datafield=None):
		coefs = None
		if datafield is not None: coefs = datafield * self._model._detJ
		if funfield is not None:  coefs = funfield(self._model._qpPhy) * self._model._detJ
		if coefs is None: raise Warning('Missing data')

		# Calculate vector
		inputs = [coefs, *self._model._nbqp, *self._model._indices, *self._model._weights]
		if self._model._dim == 2: raise Warning('Until now not done')
		if self._model._dim == 3: vector = heatsolver.wq_get_bodyheat_3d(*inputs)

		# Solve linear system with fortran
		inputs = [self._model._detJ, *self._model._nbqp, *self._model._indices, *self._model._basis, 
	    		 *self._model._weights, vector, self._nbIterPCG, self._thresholdPCG]
		start = time.process_time()
		u_interp, relres = heatsolver.mf_wq_interpolate_cp_3d(*inputs)
		stop = time.process_time()
		res_end = relres[np.nonzero(relres)][-1]
		print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, res_end))
		return u_interp

	def solveSteadyHeatProblemFT(self, b, coefs=None):
		if not self._material._isConductivityIsotropic: raise Warning('Not possible by now')
		if coefs is None: 
			inpt = self._model._qpPhy
			coefs  = self._material.eval_conductivityCoefficients(self._model._invJ, 
															self._model._detJ, inpt)
		tmp    = self.get_input4MatrixFree(table=self._boundary._thDirichletTable)
		inputs = [coefs, *tmp, b, self._nbIterPCG, self._thresholdPCG, self._methodPCG]

		if self._model._dim == 2: raise Warning('Until now not done')
		if self._model._dim == 3: sol, residue = heatsolver.mf_wq_steady_heat_3d(*inputs)

		return sol, residue
	
	def solveLinearTransientHeatProblemFT(self, dt, b, theta=1.0, Ccoefs=None, Kcoefs=None, **kwargs):
		
		if Ccoefs is None: 
			temperature = kwargs.get('temperature')
			inpt = self._model._qpPhy
			inpt = np.row_stack((inpt, temperature*np.ones(np.size(inpt, axis=1))))
			Ccoefs  = self._material.eval_capacityCoefficients(self._model._detJ, inpt)
			
		if Kcoefs is None:
			temperature = kwargs.get('temperature')
			inpt = self._model._qpPhy
			inpt = np.row_stack((inpt, temperature*np.ones(np.size(inpt, axis=1))))
			Kcoefs  = self._material.eval_conductivityCoefficients(self._model._invJ, 
						self._model._detJ, inpt)
		
		tmp    = self.get_input4MatrixFree(table=self._boundary._thDirichletTable)
		inputs = [Ccoefs, Kcoefs, *tmp, b, theta*dt, self._nbIterPCG, self._thresholdPCG, self._methodPCG]

		if self._model._dim == 2: raise Warning('Until now not done')
		if self._model._dim == 3: sol, residue = heatsolver.mf_wq_lineartransient_heat_3d(*inputs)

		return sol, residue

	# Solve using python
	def solveNLTransientHeatProblemPy(self, Tinout, time_list, Fext, theta=1.0):
		m, n = np.shape(Tinout)
		nbSteps         = len(time_list)
		nbctrlpts_total = self._model._nbctrlpts_total
		if n != nbSteps: raise Warning('Not possible')
		if m != nbctrlpts_total: raise Warning('Not possible')
		dod, _, dof = self._boundary.getThermalBoundaryConditionInfo()

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
				TTinterp = self.interpolateTemperature(TTn1)
			
				# Compute internal force
				inpt = self._model._qpPhy
				inpt = np.row_stack((inpt, TTinterp*np.ones(np.size(inpt, axis=1))))
				Ccoefs  = self._material.eval_capacityCoefficients(self._model._detJ, inpt)
				Kcoefs  = self._material.eval_conductivityCoefficients(self._model._invJ, self._model._detJ, inpt)

				CdTemp = self.eval_mfCapacity(VVn1, coefs=Ccoefs, table=np.zeros((3, 2), dtype=bool))
				KTemp  = self.eval_mfConductivity(TTn1, coefs=Kcoefs, table=np.zeros((3, 2), dtype=bool))
				Fint   = KTemp + CdTemp

				# Compute residue
				ddFF    = Fstep - Fint
				ddFFdof = ddFF[dof]
				resNL   = np.sqrt(np.dot(ddFFdof, ddFFdof))
				print('NR error: %.5f' %resNL)
				if resNL <= self._thresholdNR: break

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

class mechaproblem():
	def __init__(self, material: mechamat, model: part, boundary: step, **kwargs):
		self._material = material
		self._model    = model
		self._boundary = boundary
		self._kwargs   = kwargs
		self.__extractInfo()
		return
	
	def __extractInfo(self):
		kwargs = self._kwargs
		self._nbiterPCG    = kwargs.get('nbIterations', 100)
		self._thresholdPCG = kwargs.get('PCGThreshold', 1e-12)
		self._methodPCG    = 'JMC'
		return
	
	# Matrix free functions
	def eval_mfStiffness(self, u):
		if self._model._dim != 3: raise Warning('Until now not done')
		self._material.verifyMechanicalProperties()
		prop = [self._material._elasticmodulus, self._material._poissonratio]
		inputs = [*self._model._nbqp, *self._model._indices, 
	    			*self._model._basis, *self._model._weights, 
					self._model._invJ, self._model._detJ, prop]
		result = plasticitysolver.mf_wq_get_su_3d(*inputs, u)
		return result
	
	def eval_bodyForce(self, fun):
		if self._model._dim != 3: raise Warning('Method only for 3D geometries')
		coefs = self._material.eval_volForceCoefficients(fun, self._model._detJ, self._model._qpPhy)
		inputs = [coefs, *self._model._nbqp, *self._model._indices, *self._model._weights]
		vector = plasticitysolver.wq_get_forcevol_3d(*inputs)
		return vector
	
	def eval_surfForce(self, fun, nbFacePosition):
		" Returns force vector at the surface. In 3D: surface integrals. "

		if self._model._dim != 3: raise Warning('Method only for 3D geometries')

		def get_faceInfo(nb):
			direction = int(np.floor(nb/2))
			if nb%2 == 1: side = 1
			else: side = 0
			return direction, side

		vector = np.zeros((self._model._dim, self._model._nbctrlpts_total))
		INC_ctrlpts = self._boundary.__get_INCTable(self._model._nbctrlpts)
		INC_quadpts = self._boundary.__get_INCTable(self._model._nbqp)
		direction, side = get_faceInfo(nbFacePosition)

		# Get control points and quadrature points list
		if side == 0: 
			CPList = np.where(INC_ctrlpts[:, direction] == 0)[0]
			QPList = np.where(INC_quadpts[:, direction] == 0)[0]

		elif side == 1: 
			CPList = np.where(INC_ctrlpts[:, direction] == self._model._nbctrlpts[direction]-1)[0]
			QPList = np.where(INC_quadpts[:, direction] == self._model._nbqp[direction]-1)[0]
		
		CPList = list(np.sort(CPList))
		QPList = list(np.sort(QPList))

		# Modify Jacobien matrix
		valrange = [i for i in range(self._model._dim)]
		valrange.pop(direction)
		JJ = self._model._Jqp[:, :, QPList]
		JJ = JJ[:, valrange, :]

		# Get force values at quadrature points
		qpPhy = self._model._qpPhy[:, QPList]
		coefs = fun(qpPhy)

		# Compute surface force
		nnz, indices, weights = [], [], []
		for _ in valrange:
			nnz.append(self._model._nbqp[_]); weights.append(self._model._weights[_])
			indices.append(self._model._indices[2*_]); indices.append(self._model._indices[2*_+1]) 
		
		tmp = plasticitysolver.wq_get_forcesurf_3d(coefs, JJ, *nnz, *indices, *weights)
		vector[:, CPList] = tmp

		return vector
	
	# Solve using fortran
	def solveElasticityProblemFT(self, Fext):
		self._material.verifyMechanicalProperties()
		dod_total = deepcopy(self._boundary._mchdod)
		for i, dod in enumerate(dod_total):
			tmp = dod + 1
			dod_total[i] = tmp

		prop = [self._material._elasticmodulus, self._material._poissonratio]
		inputs = [*self._model._nbqp, *self._model._indices, *self._model._basis, 
	    			*self._model._weights, Fext, *dod_total, self._boundary._mchDirichletTable, 
					self._model._invJ, self._model._detJ, prop, self._nbiterPCG, 
					self._thresholdPCG, self._methodPCG]
		displacement, residue = plasticitysolver.mf_wq_elasticity_3d(*inputs)

		return displacement, residue

	# Solve using python
	def solvePlasticityProblem(self, Fext): 
		return

