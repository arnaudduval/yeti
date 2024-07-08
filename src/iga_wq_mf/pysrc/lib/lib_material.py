from .__init__ import *

def clean_dirichlet(A, dod):
	""" Set to 0 (Dirichlet condition) the values of an array using the indices in each dimension
		A is actually a vector arranged following each dimension [Au, Av, Aw]
	"""
	dim = np.size(A, axis=0)
	for i in range(dim): A[i, dod[i]] = 0.0
	return

def block_dot_product(A, B):
	""" Computes dot product of A and B. 
		Both are actually vectors arranged following each dimension
		A = [Au, Av, Aw] and B = [Bu, Bv, Bw]. Then A.B = Au.Bu + Av.Bv + Aw.Bw
	"""
	d = min(np.size(A, axis=0), np.size(B, axis=0))
	result = 0.0
	for i in range(d): result += A[i, :] @ B[i, :]
	return result

class material():
	""" WARNING: Depending on the algorithm (for heat transfer and elasticity), 
		the functions for isotropic materials may not work.
	"""

	def __init__(self):
		self.density = None
		return		
	
	def setScalarProperty(self, inpt, isIsotropic=False):
		if isIsotropic:
			# Isotropic material 
			assert np.isscalar(inpt), 'Not possible'
			prop = lambda x: inpt*np.ones(len(x))

		elif callable(inpt):
			# Anisotropic material using a continuous function
			prop = lambda x: inpt(x)
		else:
			raise Warning('Not implemented')
		return prop
	
	def setTensorProperty(self, inpt, shape=3, isIsotropic=False):

		def create3ArrayFrom2Array(inpt, x):
			y = np.zeros((*np.shape(inpt), len(x)))
			for i in range(np.shape(inpt)[0]):
				for j in range(np.shape(inpt)[1]):
					y[i, j, :] = inpt[i, j]
			return y

		if isIsotropic:
			# Isotropic material 
			if np.isscalar(inpt):
				prop = lambda x: create3ArrayFrom2Array(inpt*np.eye(shape), x)
			else:
				prop = lambda x: create3ArrayFrom2Array(inpt, x)
		elif callable(inpt):
			# Anisotropic material using a continuous function
			prop = lambda x: inpt(x)
		else:
			raise Warning('Not implemented')
		return prop
	
	def addDensity(self, inpt, isIsotropic):
		self.density = self.setScalarProperty(inpt, isIsotropic=isIsotropic)
		return

class heatmat(material):
	""" In our work we consider nonlinar materials, 
	meaning that its thermal properties could change depending on the position, temperature, etc.
	"""
	def __init__(self, matArgs=dict()):
		material.__init__(self)
		self.capacity     = None
		self.conductivity = None
		self.relaxation   = None
		self._isCapacityIsotropic     = False
		self._isConductivityIsotropic = False
		self._isRelaxationIsotropic = False
		self.refTemp = matArgs.get('ref_temp', 300)
		return
	
	def addCapacity(self, inpt, isIsotropic):
		if isIsotropic: self._isCapacityIsotropic = True
		self.capacity = self.setScalarProperty(inpt, isIsotropic=isIsotropic)
		return
	
	def addConductivity(self, inpt, isIsotropic, shape=3):
		if isIsotropic: self._isConductivityIsotropic = True
		self.conductivity = self.setTensorProperty(inpt, shape=shape, isIsotropic=isIsotropic)
		return
	
	def addRelaxation(self, inpt, isIsotropic):
		if isIsotropic: self._isRelaxationIsotropic = True
		self.relaxation = self.setScalarProperty(inpt, isIsotropic=isIsotropic)
		return
	
	def addCapacityDers(self, inpt, isIsotropic):
		self.capacityDers = self.setScalarProperty(inpt, isIsotropic=isIsotropic)
		return
	
	def addConductivityDers(self, inpt, isIsotropic, shape=3):
		self.conductivityDers = self.setTensorProperty(inpt, shape=shape, isIsotropic=isIsotropic)
		return

class mechamat(material):
	""" In our work we only consider isotropic material (with Lamé parameters). 
		For plasticity, we consider the J2 model with nonlinear isotropic hardening and 
		a Chaboche kinematic hardening. In that way, linear and Armstrong-Frederick hardening are subcases. 
	"""

	def __init__(self, matArgs:dict):
		material.__init__(self)
		self.thermalExpansion = matArgs.get('thermal_expansion', 1.0)
		self.elasticModulus = matArgs.get('elastic_modulus', None)
		self.poissonRatio   = matArgs.get('poisson_ratio', None)
		self.elasticLimit   = matArgs.get('elastic_limit', 1e15)
		if any(prop is None for prop in [self.elasticModulus, self.elasticLimit, self.poissonRatio]): raise Warning('Mechanics not well defined')
		self.__setExtraMechanicalProperties()

		isohardArgs = matArgs.get('isoHardLaw', {})
		self._isoHardening = self.isotropicHardening(self.elasticLimit, isohardArgs=isohardArgs)
		kineLaw = matArgs.get('kineHardLaw', {})
		chabocheTable = kineLaw.get('parameters', np.array([[0, 0]]))
		# By default if only one parameter we consider a linear kinematic hardening
		# Armstrong Frederick hardening is only considering a single row [[c, d]]
		self._chabocheTable = np.atleast_2d(chabocheTable)
		self._chabocheNBparameters = np.size(self._chabocheTable, axis=0)

		self.parametersPreCalc1D = self.__parametersPreCalc1D
		self.parametersPreCalc3D = self.__parametersPreCalc3D
		self.consistentTangentAlgorithm3D = self.__consistentTangentAlgorithm3D
		self.consistentTangentAlgorithm1D = self.__consistentTangentAlgorithm1D		

		return
		
	def __setExtraMechanicalProperties(self):
		E  = self.elasticModulus
		nu = self.poissonRatio
		self.lame_lambda, self.lame_mu, self.lame_bulk = None, None, None
		if E is not None and nu is not None:
			lamb = nu*E/((1+nu)*(1-2*nu))
			mu = E/(2*(1+nu))
			bulk = lamb + 2.0/3.0*mu
			self.lame_lambda = lamb
			self.lame_mu = mu
			self.lame_bulk = bulk
		return
	
	def evalElasticStress(self, strain):
		nvoigt = np.size(strain, axis=0)
		assert nvoigt == 6, 'Try another method'
		traceStrain = computeTrace4All(strain)
		stress = 2*self.lame_mu*strain
		for i in range(3): stress[i, :] += self.lame_lambda*traceStrain
		return stress
	
	class isotropicHardening():
		def __init__(self, elasticlimit, isohardArgs:dict):
			self.__elasticlimit = elasticlimit
			isoname = isohardArgs.get('name', 'none').lower()
			if   isoname == 'linear': self.__setIsoLinearModel(isohardArgs)
			elif isoname == 'swift' : self.__setIsoSwiftModel(isohardArgs)
			elif isoname == 'voce'  : self.__setIsoVoceModel(isohardArgs)
			elif isoname == 'none'  : self.__setIsoNoneModel()
			else: raise Warning('Unknown method')
			return	
		
		def __setIsoNoneModel(self):
			self._isohardfun = lambda a: 1e6*self.__elasticlimit*np.ones(np.size(a))
			self._isohardfunders = lambda a: np.zeros(np.size(a))
			return
		
		def __setIsoLinearModel(self, args:dict):
			Eiso = args.get('Eiso', None)
			self._isohardfun = lambda a: self.__elasticlimit + Eiso*a 
			self._isohardfunders = lambda a: Eiso*np.ones(np.size(a))
			return
		
		def __setIsoSwiftModel(self, args:dict):
			e0 = args.get('e0', None)
			n = args.get('n', None)
			self._isohardfun = lambda a: self.__elasticlimit*(1 + a/e0)**n
			self._isohardfunders = lambda a: self.__elasticlimit*n/e0*(1 + a/e0)**(n - 1)
			return
		
		def __setIsoVoceModel(self, args:dict):
			ssat = args.get('ssat', None)
			beta = args.get('beta', None)
			self._isohardfun = lambda a: self.__elasticlimit +  ssat*(1.0 - np.exp(-beta*a))
			self._isohardfunders = lambda a: ssat*beta*np.exp(-beta*a) 
			return
	
	# 3D (or 2D even if within their name one could find '3D')
	
	def __sumOverChabocheTable(self, dgamma, back):
		meanback = np.zeros(np.shape(back[0, :, :])); hatback = np.zeros(np.shape(back[0, :, :]))
		const1, const2 = np.zeros(shape=np.shape(dgamma)), np.zeros(shape=np.shape(dgamma))
		for i in range(self._chabocheNBparameters):
			[ci, di] = self._chabocheTable[i, :]
			meanback += back[i, :, :]/(1 + di*dgamma)
			hatback += di*back[i, :, :]/(1 + di*dgamma)**2
			const1 += ci/(1 + di*dgamma)
			const2 += ci/(1 + di*dgamma)**2
		return meanback, hatback, const1, const2

	def __parametersPreCalc3D(self, stress_trial, back_n0, plseq_n0, nbIter=50, threshold=1e-8):

		dgamma = np.zeros(shape=np.shape(plseq_n0)); straineq_n1 = np.copy(plseq_n0)
		theta = np.zeros(shape=np.shape(plseq_n0)); thetatilde = np.zeros(shape=np.shape(plseq_n0))
		for k in range(nbIter):
			meanback, hatback, const1, const2 = self.__sumOverChabocheTable(dgamma, back_n0)
			hatshifted = computeDeviatoric4All(stress_trial - meanback)
			normhatshifted = computeSymTensorNorm4All(hatshifted)
			funyield = (np.sqrt(3.0/2.0)*normhatshifted - 3.0/2.0*(2*self.lame_mu + const1)*dgamma 
						- self._isoHardening._isohardfun(straineq_n1)
			)
			resNLj = np.sqrt(np.dot(np.ravel(funyield), np.ravel(funyield)))
			if k == 0: resNL0 = resNLj
			if resNLj <= max([threshold*resNL0, 1e-12]): break
			normal = hatshifted/normhatshifted
			dersfunyield = (np.sqrt(3.0/2.0)*computeSymDoubleContraction4All(normal, hatback) 
				- 3.0/2.0*(2*self.lame_mu + const2) - self._isoHardening._isohardfunders(straineq_n1)
			)	
			dgamma -= funyield/dersfunyield; straineq_n1 = plseq_n0 + dgamma
			thetatilde = -3*self.lame_mu/dersfunyield
			theta = 2*self.lame_mu*dgamma/(normhatshifted)*np.sqrt(3.0/2.0)

		return dgamma, hatback, normal, theta, thetatilde
	
	def __consistentTangentAlgorithm3D(self, nnz, isElastic, plsVars={}):
		if isElastic:
			mechArgs = np.zeros((2, nnz))
			mechArgs[0, :] = self.lame_lambda; mechArgs[1, :] = self.lame_mu
		else:
			plsInd = plsVars['plsInd']; normal = plsVars['normal']; hatback = plsVars['hatback']
			theta  = np.ravel(plsVars['theta']); thetatilde = np.ravel(plsVars['thetatilde'])
			nvoigt = np.size(hatback, axis=0)
			mechArgs = np.zeros((4+2*nvoigt, nnz))
			mechArgs[0, :] = self.lame_lambda; mechArgs[0, plsInd] += 2.0/3.0*self.lame_mu*theta 
			mechArgs[1, :] = self.lame_mu; mechArgs[1, plsInd] = self.lame_mu*(1 - theta)
			mechArgs[2, plsInd] = -2*self.lame_mu*(thetatilde - theta)
			mechArgs[3, plsInd] = -np.sqrt(2.0/3.0)*theta*thetatilde
			mechArgs[4:4+nvoigt, plsInd] = normal
			mechArgs[4+nvoigt:, plsInd]  = hatback
		return mechArgs

	def J2returnMappingAlgorithm3D(self, strain_n1, plasticstrain_n0, plseq_n0, back_n0, isElasticMatrix=False, threshold=1e-8, nvoigtreal=6):
		""" Return mapping algorithm for multidimensional rate-independent plasticity. 
		"""		
		nnz = np.shape(strain_n1)[1]
		isElasticLoad = True; output = {}

		# Compute trial stress
		strain_trial = strain_n1 - plasticstrain_n0
		stress_trial = self.evalElasticStress(strain_trial)

		# Compute shifted stress
		shifted_trial = computeDeviatoric4All(stress_trial - np.sum(back_n0, axis=0))

		# Check yield status
		norm_shifted_trial = computeSymTensorNorm4All(shifted_trial)
		yield_trial = np.sqrt(3.0/2.0)*norm_shifted_trial - self._isoHardening._isohardfun(plseq_n0)
		stress_n1 = np.copy(stress_trial); plasticstrain_n1= np.copy(plasticstrain_n0); 
		plseq_n1 = np.copy(plseq_n0); back_n1 = np.copy(back_n0); plsVars = {}

		plsInd = np.nonzero(np.ravel(yield_trial)>threshold*self.elasticLimit)[0]
		if np.size(plsInd) > 0:

			isElasticLoad = False

			# Compute plastic-strain increment
			dgamma, hatback, normal, theta, thetatilde = self.parametersPreCalc3D(stress_trial[:, plsInd], 
															back_n0[:, :, plsInd], plseq_n0[:, plsInd])

			# Update internal hardening variable
			plseq_n1[:, plsInd] += dgamma

			# Update stress
			stress_n1[:, plsInd] -= 2*self.lame_mu*np.sqrt(3.0/2.0)*dgamma*normal
			
			# Update plastic strain
			plasticstrain_n1[:, plsInd] += np.sqrt(3.0/2.0)*dgamma*normal

			# Update backstress
			for i in range(self._chabocheNBparameters):
				[ci, di] = self._chabocheTable[i, :]
				for j, ind in enumerate(plsInd):
					back_n1[i, :, ind] = (back_n0[i, :, ind] + np.sqrt(3.0/2.0)*ci*dgamma[:, j]*normal[:, j])/(1. + di*dgamma[:, j])
			
			normaltmp = np.copy(normal); hatbacktmp = np.copy(hatback)
			if nvoigtreal == 3:
				normaltmp = np.zeros((nvoigtreal, len(plsInd))); hatbacktmp = np.zeros((nvoigtreal, len(plsInd)))
				normaltmp[0:2, :] = normal[0:2, :]; normaltmp[-1, :] = normal[3, :]
				hatbacktmp[0:2, :] = hatback[0:2, :]; hatbacktmp[-1, :] = hatback[3, :]
			plsVars = {'plsInd': plsInd, 'normal': normaltmp, 'hatback': hatbacktmp, 'theta': theta, 'thetatilde': thetatilde}

		mechArgs = self.consistentTangentAlgorithm3D(nnz, isElasticMatrix+isElasticLoad, plsVars)
		output = {'stress': stress_n1, 'plastic': plasticstrain_n1, 'plseq': plseq_n1, 'back': back_n1, 'mechArgs': mechArgs}

		return output, isElasticLoad
	
	# 1D
	def __parametersPreCalc1D(self, stress_trial, back_n0, plseq_n0, nbIter=50, threshold=1e-8):
		
		dgamma = np.zeros(shape=np.shape(plseq_n0)); plseq_n1 = np.copy(plseq_n0); theta = np.zeros(shape=np.shape(plseq_n0))
		for k in range(nbIter):
			meanback, hatback, const1, const2 = self.__sumOverChabocheTable(dgamma, back_n0)
			hatshifted = stress_trial - meanback
			funyield = np.abs(hatshifted) - (self.elasticModulus + const1)*dgamma - self._isoHardening._isohardfun(plseq_n1)
			resNLj = np.sqrt(np.dot(np.ravel(funyield), np.ravel(funyield)))
			if k == 0: resNL0 = resNLj
			if resNLj <= max([threshold*resNL0, 1e-12]): break
			normal = np.sign(hatshifted)
			dersfunyield = normal*hatback - (self.elasticModulus + const2) - self._isoHardening._isohardfunders(plseq_n1)			
			dgamma -= funyield/dersfunyield; plseq_n1 = plseq_n0 + dgamma
			theta = -self.elasticModulus/dersfunyield
		return dgamma, hatback, normal, theta
	
	def __consistentTangentAlgorithm1D(self, nnz, isElastic, plsVars={}):
		mechArgs = np.zeros((1, nnz))
		mechArgs[0, :] = self.elasticModulus
		if not isElastic:
			plsInd = plsVars['plsInd']; theta = plsVars['theta']
			mechArgs[0, plsInd] = self.elasticModulus*(1 - theta)
		return mechArgs
	
	def J2returnMappingAlgorithm1D(self, strain_n1, plasticstrain_n0, plseq_n0, back_n0, isElasticMatrix=False, threshold=1e-8):
		""" Return mapping algorithm for multidimensional rate-independent plasticity. 
		"""		
		isElasticLoad = True; output = {}

		# Compute trial stress
		strain_trial = strain_n1 - plasticstrain_n0
		stress_trial = self.elasticModulus*strain_trial

		# Compute shifted stress
		shifted_trial = stress_trial - np.sum(back_n0, axis=0)

		# Check yield status
		norm_shifted_trial = np.abs(shifted_trial)
		yield_trial = norm_shifted_trial - self._isoHardening._isohardfun(plseq_n0)
		stress_n1 = np.copy(stress_trial); plasticstrain_n1= np.copy(plasticstrain_n0)
		plseq_n1 = np.copy(plseq_n0); back_n1 = np.copy(back_n0); plsVars = {}

		plsInd = np.nonzero(np.ravel(yield_trial)>threshold*self.elasticLimit)[0]
		if np.size(plsInd) > 0:

			isElasticLoad = False

			# Compute plastic-strain increment
			dgamma, hatback, normal, theta = self.parametersPreCalc1D(stress_trial[:, plsInd], 
															back_n0[:, :, plsInd], plseq_n0[:, plsInd])

			# Update internal hardening variable
			plseq_n1[:, plsInd] += dgamma

			# Update stress
			stress_n1[:, plsInd] -= self.elasticModulus*dgamma*normal
			
			# Update plastic strain
			plasticstrain_n1[:, plsInd] += dgamma*normal

			# Update backstress
			for i in range(self._chabocheNBparameters):
				[ci, di] = self._chabocheTable[i, :]
				for j, ind in enumerate(plsInd):
					back_n1[i, :, ind] = (back_n0[i, :, ind] + ci*dgamma[:, j]*normal[:, j])/(1. + di*dgamma[:, j])

			plsVars = {'plsInd': plsInd, 'normal': normal, 'hatback': hatback, 'theta': theta}

		mechArgs = self.consistentTangentAlgorithm1D(np.size(plseq_n0, axis=1), isElasticMatrix+isElasticLoad, plsVars)
		output = {'stress': stress_n1, 'plastic': plasticstrain_n1, 'plseq': plseq_n1, 'back': back_n1, 'mechArgs': mechArgs}

		return output, isElasticLoad

def computeTrace4All(arrays):
	dim = 3
	nvoigt, nnz = np.shape(arrays)
	assert nvoigt == 6, 'Try another method'
	trace = np.zeros(nnz)
	for i in range(dim):
		trace += arrays[i, :]
	return trace

def computeDeviatoric4All(arrays):
	dim = 3
	devarray = np.copy(arrays)
	trace = computeTrace4All(arrays)/3.0
	for i in range(dim): devarray[i, :] -= trace
	return devarray

def computeSymDoubleContraction4All(arrays1, arrays2):
	dim, nvoigt = 3, 6
	nnz = np.size(arrays1, axis=1) # or arrays2
	s = np.zeros(nnz)
	for i in range(dim):
		s += arrays1[i, :]*arrays2[i, :]
	for i in range(dim, nvoigt):
		s += 2*arrays1[i, :]*arrays2[i, :]
	return s

def computeSymTensorNorm4All(arrays):
	s = computeSymDoubleContraction4All(arrays, arrays)
	return np.sqrt(s)

