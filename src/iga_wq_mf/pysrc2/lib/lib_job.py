from lib.__init__ import *
from lib.lib_base import eraseRowsCSR
from lib.lib_material import *
from lib.lib_model import *
from lib.lib_step import *

class heatproblem():

	def __init__(self, mat: thermomat, geo: part, boundcond: step, **kwargs):
		self._material = mat
		self._geometry = geo
		self._boundary = boundcond
		self._kwargs   = kwargs
		self.__extractInfo()
		return
	
	def __extractInfo(self):
		kwargs = self._kwargs
		self._nbIterPCG    = kwargs.get('nbIterations', 100)
		self._thresholdPCG = kwargs.get('PCGThreshold', 1e-12)
		self._methodPCG    = 'FDC'
		return
	
	# Matrix free functions
	def get_input4MatrixFree(self, table=None):
		" Returns necessary inputs to compute the product between a matrix and a vector "
		
		if table is None: table = self._boundary._thDirichletTable
		indices, basis, weights = [], [], []
		for i in range(self._geometry._dim):
			# Select data
			if np.array_equal(table[i, :], [0, 0]): rows2erase = []
			if np.array_equal(table[i, :], [0, 1]): rows2erase = [-1]
			if np.array_equal(table[i, :], [1, 0]): rows2erase = [0]
			if np.array_equal(table[i, :], [1, 1]): rows2erase = [0, -1]
			indi_t, indj_t, data_t = eraseRowsCSR(rows2erase, 
									self._geometry._indices[2*i], self._geometry._indices[2*i+1],  
									[self._geometry._basis[i], self._geometry._weights[i]])
			
			# Extract data and append to list
			[basist, weightst] = data_t
			indices.append(indi_t); indices.append(indj_t) 
			basis.append(basist); weights.append(weightst)

		inputs = [*self._geometry._nbqp, *indices, *basis, *weights]

		return inputs

	def eval_mfConductivity(self, input, u, table=None):

		inputs = self.get_input4MatrixFree(table=table)
		coefs  = self._material.eval_conductivityCoefficients(self._geometry._invJ, self._geometry._detJ, input)
		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: result = solver.mf_wq_get_ku_3d_py(coefs, *inputs, u)

		return result
	
	def eval_mfCapacity(self, input, u, table=None): 

		coefs = self._material.eval_capacityCoefficients(self._geometry._detJ, input)
		inputs = self.get_input4MatrixFree(table=table)
		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: result = solver.mf_wq_get_cu_3d_py(coefs, *inputs, u)

		return result
	
	def eval_mfCondCap(self, Kinput, Cinput, u, table=None, alpha=1.0, beta=1.0):

		Kcoefs = self._material.eval_conductivityCoefficients(self._geometry._invJ, self._geometry._detJ, Kinput)
		Ccoefs = self._material.eval_capacityCoefficients(self._geometry._detJ, Cinput)
		inputs = self.get_input4MatrixFree(table=table)
		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: result = solver.mf_wq_get_kcu_3d_py(Ccoefs, Kcoefs, *inputs, u, alpha, beta)

		return result

	def eval_heatForce(self, fun, indi=None): 

		if indi is None: indi = np.arange(self._geometry._nbctrlpts, dtype=int)
		coefs = self._material.eval_heatForceCoefficients(fun, self._geometry._detJ, self._geometry._qpPhy)
		inputs = [coefs, *self._geometry._nbqp, *self._geometry._indices, *self._geometry._weights]
		if self._dim == 2: vector = assembly.wq_get_source_2d(*inputs)[indi]
		if self._dim == 3: vector = assembly.wq_get_source_3d(*inputs)[indi]

		return vector
	
	# Solve using fortran
	def solveInterpolationProblem(self, funfield=None, datafield=None):
		coefs = None
		if datafield is not None: coefs = datafield * self._geometry._detJ
		if funfield is not None:  coefs = funfield(self._geometry._qpPhy) * self._geometry._detJ
		if coefs is None: raise Warning('Missing data')

		# Calculate vector
		inputs = [coefs, *self._geometry._nbqp, *self._geometry._indices, *self._geometry._weights]
		if self._dim == 2: raise Warning('Until now not done')
		if self._dim == 3: vector = assembly.wq_get_source_3d(*inputs)

		# Solve linear system with fortran
		inputs = [self._geometry._detJ, *self._geometry._nbqp, *self._geometry._indices, *self._geometry._basis, 
	    		 *self._geometry._weights, vector, self._nbIterPCG, self._thresholdPCG]
		start = time.process_time()
		u_interp, relres = solver.mf_wq_interpolate_cp_3d(*inputs)
		stop = time.process_time()
		res_end = relres[np.nonzero(relres)][-1]
		print('Interpolation in: %.3e s with relative residue %.3e' %(stop-start, res_end))
		return u_interp

	def solveSteadyHeatProblem(self):
		
		return
	
	def solveTransientHeatProblem(self):

		return

class mechaproblem():
	def __init__(self, mat: mechamat, geo: part, boundcond: step, **kwargs):
		self._material = mat
		self._geometry = geo
		self._boundary = boundcond
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
		if self._dim != 3: raise Warning('Until now not done')
		self._material.verifyMechanicalProperties()
		prop = [self._material._elasticmodulus, self._material._poissonratio]
		inputs = [*self._geometry._nbqp, *self._geometry._indices, 
	    			*self._geometry._basis, *self._geometry._weights, 
					self._geometry._invJ, self._geometry._detJ, prop]
		result = elastoplasticity.mf_wq_get_su_3d_py(*inputs, u)
		return result
	
	def eval_bodyForce(self, fun):
		if self._dim != 3: raise Warning('Method only for 3D geometries')
		coefs = self._material.eval_volForceCoefficients(fun, self._geometry._detJ, self._geometry._qpPhy)
		inputs = [coefs, *self._geometry._nbqp, *self._geometry._indices, *self._geometry._weights]
		vector = elastoplasticity.wq_get_forcevol_3d(*inputs)
		return vector
	
	# Solve using fortran
	def solveElasticityProblem(self, Fext):
		self._material.verifyMechanicalProperties()
		dod = deepcopy(self._boundary._mchdod)
		for i in range(len(dod)):
			tmp = np.array(dod[i]); tmp += 1
			dod[i] = list(tmp)

		prop = [self._material._elasticmodulus, self._material._poissonratio]
		inputs = [*self._geometry._nbqp, *self._geometry._indices, *self._geometry._basis, 
	    			*self._geometry._weights, Fext, *dod, self._boundary._mchDirichletTable, 
					self._geometry._invJ, self._geometry._detJ, prop, self._nbiterPCG, 
					self._thresholdPCG, self._methodPCG]
		displacement, residue = elastoplasticity.mf_wq_elasticity_3d_py(*inputs)
		return displacement, residue

	def solvePlasticityProblem(self, Fext): 
		return

