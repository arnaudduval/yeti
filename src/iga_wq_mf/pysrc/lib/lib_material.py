from .__init__ import *

class material():

	def __init__(self):
		return		
	
	def setScalarProperty(self, inpt, isIsotropic=False):
		if isIsotropic:
			# Isotropic material 
			if np.isscalar(inpt): 		
				prop = lambda x: inpt*np.ones(np.max(np.shape(x)))
			else: 
				raise Warning('Not possible')
		elif callable(inpt):
			# Anisotropic material using a continuous function
			prop = lambda x: inpt(x)
		else:
			raise Warning('Not implemented')
		return prop
	
	def setTensorProperty(self, inpt, shape=(3, 3), isIsotropic=False):

		def create3ArrayFrom2Array(inpt, x):
			lenx = np.max(np.shape(x))
			y = np.zeros((*np.shape(inpt), lenx))
			for i in range(np.shape(inpt)[0]):
				for j in range(np.shape(inpt)[1]):
					y[i, j, :] = inpt[i, j]
			return y

		if isIsotropic:
			# Isotropic material 
			if np.isscalar(inpt):
				prop = lambda x: inpt*np.eye((*shape, np.max(np.shape(x))))
			else:
				prop = lambda x: create3ArrayFrom2Array(inpt, x)
		elif callable(inpt):
			# Anisotropic material using a continuous function
			prop = lambda x: inpt(x)
		else:
			raise Warning('Not implemented')
		return prop

class thermomat(material):
	def __init__(self):
		super().__init__()
		self.density      = None
		self.capacity     = None
		self.conductivity = None

		self._isDensityIsotropic      = False
		self._isCapacityIsotropic     = False
		self._isConductivityIsotropic = False
		return
	
	def addDensity(self, inpt, isIsotropic):
		if isIsotropic: self._isDensityIsotropic = True
		self.density     = super().setScalarProperty(inpt, isIsotropic=isIsotropic)
		return
	
	def addCapacity(self, inpt, isIsotropic):
		if isIsotropic: self._isCapacityIsotropic = True
		self.capacity    = super().setScalarProperty(inpt, isIsotropic=isIsotropic)
		return
	
	def addConductivity(self, inpt, isIsotropic, shape=(3, 3)):
		if isIsotropic: self._isConductivityIsotropic = True
		self.conductivity = super().setTensorProperty(inpt, shape=shape, isIsotropic=isIsotropic)
		return

class plasticLaw():
	def __init__(self, elasticmodulus, elasticlimit, plasticArgs:dict):
		self.__elasticlimit   = elasticlimit
		self.__elasticmodulus = elasticmodulus
		self._Kfun = None; self._Kderfun = None
		self._Hfun = None; self._Hderfun = None
		
		self._plasticArgs = plasticArgs
		if   plasticArgs['name'] == 'linear': self.__setLinearModel(plasticArgs)
		elif plasticArgs['name'] == 'swift' : self.__setSwiftModel(plasticArgs)
		elif plasticArgs['name'] == 'voce'  : self.__setVoceModel(plasticArgs)
		else: raise Warning('Unknown method')
		funlist = [self._Hfun, self._Hderfun, self._Kfun, self._Kderfun]
		if any(fun is None for fun in funlist): raise Warning('Something went wrong')
		return	
	
	def __setLinearModel(self, plasticArgs:dict):
		theta = plasticArgs.get('theta', None)
		Hbar  = plasticArgs.get('Hbar', None)
		self._Kfun = lambda a: self.__elasticlimit + theta*Hbar*a 
		self._Hfun = lambda a: (1 - theta)*Hbar*a
		self._Kderfun = lambda a: theta*Hbar*np.ones(np.size(a))
		self._Hderfun = lambda a: (1 - theta)*Hbar*np.ones(np.size(a))
		return
	
	def __setSwiftModel(self, plasticArgs:dict):
		K = plasticArgs.get('K', None)
		n = plasticArgs.get('exp', None)
		def localKderfun(a):
			with np.errstate(divide='ignore'):
				y = np.where(a==0.0, 1e10*self.__elasticmodulus, (self.__elasticmodulus/K)*n*(a/K)**(n-1))
			return y
		self._Kfun = lambda a: self.__elasticlimit + self.__elasticmodulus*(a/K)**n
		self._Hfun = lambda a: np.zeros(np.size(a))
		self._Kderfun = lambda a: localKderfun(a)
		self._Hderfun = lambda a: np.zeros(np.size(a))
		return
	
	def __setVoceModel(self, plasticArgs:dict):
		theta = plasticArgs.get('theta', None)
		Hbar  = plasticArgs.get('Hbar', None)
		Kinf  = plasticArgs.get('Kinf', None)
		delta = plasticArgs.get('delta', None)
		self._Kfun = lambda a: self.__elasticlimit + theta*Hbar*a + Kinf*(1.0 - np.exp(-delta*a))
		self._Hfun = lambda a: (1 - theta)*Hbar*a
		self._Kderfun = lambda a: theta*Hbar + Kinf*delta*np.exp(-delta*a) 
		self._Hderfun = lambda a: (1 - theta)*Hbar*np.ones(np.size(a))
		return
	
	def __setExponential(self, plasticArgs:dict):
		theta = plasticArgs.get('theta', None)
		Hbar  = plasticArgs.get('Hbar', None)
		Kinf  = plasticArgs.get('Kinf', None)
		Kzero = plasticArgs.get('Kzero', None)
		delta = plasticArgs.get('delta', None)
		
		def h(a):
			return  Kinf - (Kinf - Kzero)*np.exp(-delta*a) + Hbar*a
		def dh(a):
			return (Kinf - Kzero)*delta*np.exp(-delta*a) + Hbar

		self._Kfun = lambda a: theta*h(a)
		self._Hfun = lambda a: (1 - theta)*h(a)
		self._Kderfun = lambda a: theta*dh(a)
		self._Hderfun = lambda a: (1 - theta)*dh(a)
		return
	
class mechamat(material):
	# Eventually this class should be similar to thermomat, but for the moment let's say it works !
	def __init__(self, matArgs:dict):
		super().__init__()
		self.density        = matArgs.get('density', None)
		self.elasticmodulus = matArgs.get('elastic_modulus', None)
		self.poissonratio   = matArgs.get('poisson_ratio', None)
		self.elasticlimit   = matArgs.get('elastic_limit', None)
		if any(prop is None for prop in [self.elasticmodulus, self.elasticlimit, self.poissonratio]): 
			raise Warning('Mechanics not well defined')

		self.plasticLaw            = None
		self._isPlasticityPossible = False
		tmp = matArgs.get('plasticLaw', None)
		if isinstance(tmp, dict): 
			self._isPlasticityPossible = True
			self.plasticLaw = plasticLaw(self.elasticmodulus, self.elasticlimit, tmp)
		self.__setExtraMechanicalProperties()
		return
	
	def __setExtraMechanicalProperties(self):
		E  = self.elasticmodulus
		nu = self.poissonratio
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
		Tstrain = array2symtensor4All(strain, 2)
		traceStrain = evalTrace4All(strain, 2)
		Tstress = 2*self.lame_mu*Tstrain
		for i in range(2): Tstress[i, i, :] += self.lame_lambda*traceStrain
		stress  = symtensor2array4All(Tstress, 2)
		return stress
	
	def returnMappingAlgorithm(self, strain, pls, a, b, threshold=1e-9):
		""" Return mapping algorithm for multidimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		def computeDeltaGamma(law:plasticLaw, lame_mu, a_n0, eta_trial, nbIter=20, threshold=1e-9):
			if law._plasticArgs['name'] == 'linear':
				f_trial = np.linalg.norm(eta_trial, axis=(0, 1)) - np.sqrt(2.0/3.0)*law._Kfun(a_n0)
				dgamma  = f_trial/(2.0*lame_mu + 2.0*law._plasticArgs['Hbar']/3.0)
			else:
				dgamma = np.zeros(np.size(a_n0))
				a_n1 = a_n0
				for i in range(nbIter):
					dH = law._Hfun(a_n1) - law._Hfun(a_n0) 
					G  = (-np.sqrt(2.0/3.0)*law._Kfun(a_n1) + np.linalg.norm(eta_trial, axis=(0, 1)) 
						- (2.0*lame_mu*dgamma + np.sqrt(2.0/3.0)*dH))
					if np.all(np.abs(G)<=threshold): break
					dG = - 2.0*(lame_mu + (law._Hderfun(a_n1) + law._Kderfun(a_n1))/3.0)
					dgamma -= G/dG
					a_n1 = a_n0 + np.sqrt(2.0/3.0)*dgamma
			return dgamma

		nvoigt, nnz = np.shape(strain)
		if nvoigt   == 3: dim = 2
		elif nvoigt == 6: dim = 3

		output  = np.zeros((4*(nvoigt+1), nnz))
		Cep     = np.zeros((3, nnz))
		Tstrain = array2symtensor4All(strain, dim)
		Tpls    = array2symtensor4All(pls, dim)
		Tb      = array2symtensor4All(b, dim)
		traceStrain = evalTrace4All(strain, dim)
		devStrain   = Tstrain
		for i in range(dim): devStrain[i, i, :] -= 1.0/3.0*traceStrain

		# Compute trial stress
		s_trial = 2*self.lame_mu*(devStrain - Tpls)

		# Compute shifted stress
		eta_trial = s_trial - Tb

		# Check yield status
		norm_trial = np.linalg.norm(eta_trial, axis=(0, 1))
		f_trial = norm_trial - np.sqrt(2.0/3.0)*self.plasticLaw._Kfun(a)
		sigma   = s_trial
		for i in range(dim): sigma[i, i, :] += self.lame_bulk*traceStrain
		Cep[0, :] = self.lame_lambda; Cep[1, :] = self.lame_mu
		pls_new = pls; a_new = a; b_new = b
		stress  = symtensor2array4All(sigma, dim)

		plsInd = np.nonzero(f_trial>threshold)[0]
		if np.size(plsInd)>0:

			# Compute plastic-strain increment
			dgamma_plsInd = computeDeltaGamma(self.plasticLaw, self.lame_mu, a[plsInd], eta_trial[:, :, plsInd])

			# Update internal hardening variable
			a_new[plsInd] = a[plsInd] + np.sqrt(2.0/3.0)*dgamma_plsInd

			# Compute df/dsigma
			Normal_plsInd = eta_trial[:, :, plsInd]/norm_trial[plsInd]
			output[3*nvoigt+4:, plsInd] = symtensor2array4All(Normal_plsInd, dim)

			# Update stress
			sigma[:, :, plsInd] -= 2*self.lame_mu*dgamma_plsInd*Normal_plsInd
			stress[:, plsInd] = symtensor2array4All(sigma[:, :, plsInd], dim)

			# Update plastic strain
			Tpls[:, :, plsInd] += dgamma_plsInd*Normal_plsInd
			pls_new[:, plsInd] = symtensor2array4All(Tpls[:, :, plsInd], dim)

			# Update backstress
			Tb[:, :, plsInd] += np.sqrt(2.0/3.0)*(self.plasticLaw._Hfun(a_new[plsInd]) - self.plasticLaw._Hfun(a[plsInd]))*Normal_plsInd
			b_new[:, plsInd] = symtensor2array4All(Tb[:, :, plsInd], dim)

			# Update tangent coefficients
			somme = self.plasticLaw._Kderfun(a_new[plsInd]) + self.plasticLaw._Hderfun(a_new[plsInd])
			c1 = 2*self.lame_mu*dgamma_plsInd/norm_trial[plsInd]
			c2 = 1.0/(1+somme/(3*self.lame_mu)) - c1
			Cep[0, plsInd] = self.lame_lambda + 2.0/3.0*self.lame_mu*c1
			Cep[1, plsInd] = self.lame_mu*(1.0 - c1)
			Cep[2, plsInd] = -2.0*self.lame_mu*c2

		output[0:nvoigt, :] = stress; output[nvoigt:2*nvoigt, :] = pls_new; output[2*nvoigt, :] = a_new
		output[2*nvoigt+1:3*nvoigt+1, :] = b_new; output[3*nvoigt+1:3*nvoigt+4, :] = Cep
		return output
	
def clean_dirichlet(A, dod):
	""" Set to 0 (Dirichlet condition) the values of an array using the indices in each dimension
		A is actually a vector arranged following each dimension [Au, Av, Aw]
	"""
	dim = np.size(A, axis=0)
	for i in range(dim): A[i, dod[i]] = 0.0
	return

def block_dot_product(d, A, B):
	""" Computes dot product of A and B. 
		Both are actually vectors arranged following each dimension
		A = [Au, Av, Aw] and B = [Bu, Bv, Bw]. Then A.B = Au.Bu + Av.Bv + Aw.Bw
	"""
	result = 0.0
	for i in range(d): result += A[i, :] @ B[i, :]
	return result

def symtensor2array4All(tensors, dim):
	nvoigt = int(dim*(dim+1)/2); nnz = np.size(tensors, axis=2)
	array  = np.zeros((nvoigt, nnz))
	k = 0
	for i in range(dim):
		array[k, :] = tensors[i, i, :]
		k += 1
	for i in range(dim-1):
		for j in range(i+1, dim):
			array[k, :] = tensors[i, j, :]
			k += 1
	return array

def array2symtensor4All(arrays, dim):
	nnz = np.size(arrays, axis=1)
	tensor = np.zeros((dim, dim, nnz))
	k = 0
	for i in range(dim):
		tensor[i, i, :] = arrays[k, :]
		k += 1
	for i in range(dim-1):
		for j in range(i+1, dim):
			tensor[i, j, :] = arrays[k, :] 
			tensor[j, i, :] = arrays[k, :] 
			k += 1
	return tensor

def evalTrace4All(arrays, dim):
	nnz = np.size(arrays, axis=1)
	trace = np.zeros(nnz)
	for i in range(dim):
		trace += arrays[i, :]
	return trace

def computeVMStress4All(arrays, dim):
	nnz = np.size(arrays, axis=1)
	nvgt  = int(dim*(dim+1)/2)
	trace = evalTrace4All(arrays, dim)
	dev   = np.copy(arrays)
	for i in range(dim):
		dev[i, :] -= 1.0/3.0*trace
	s = np.zeros(nnz)
	for i in range(dim):
		s += arrays[i, :]*arrays[i, :]
	for i in range(dim, nvgt):
		s += 2*arrays[i, :]*arrays[i, :]
	vm = np.sqrt(3.0/2.0*s)
	return vm