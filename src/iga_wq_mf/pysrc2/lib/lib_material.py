from lib.__init__ import *

class material():

	def __init__(self):
		self._threshold = 1e-8
		return		
	
	def verifyTable(self, table, isTensor=False):
		table = np.atleast_2d(table)
		x = np.diff(table[:, 0])
		if np.any(x<self._threshold): raise Warning('Table is not well defined')
		if isTensor and np.size(table, axis=1) <= 2: raise Warning('Table not well defined')
		return
	
	def interpolateScalarProperty(table, x):
		lenx = np.max(np.shape(x))
		if np.size(table, axis=0) == 1: y = table[:, 1]*np.ones(lenx)
		else:  							y = np.interp(x, table[:, 0], table[:, 1])
		return y
	
	def interpolateTensorProperty(table, x, shape=(3, 3)):
		if np.size(table, axis=1) != np.prod(shape)+1: raise Warning('Not possible')
		lenx = np.max(np.shape(x))
		y = np.zeros((*shape, lenx))
		for i in range(shape[0]):
			for j in range(shape[1]):
				k = 1 + j + shape[1]*i
				y[i, j, :] = np.interp(x, table[:, 0], table[:, k])
		return y
	
	def setScalarProperty(self, inpt, isIsotropic=False):
		if isIsotropic:
			# Isotropic material (position and temperature independent)
			if np.isscalar(inpt): 		
				prop = lambda x: inpt*np.ones(np.max(np.shape(x)))
			else: 
				raise Warning('Not possible')
		elif callable(inpt):
			# Anisotropic material (temperature independent but position dependent)
			prop = lambda x: inpt(x)
		elif len(inpt.shape) > 1:
			# Anisotropic material (position independent but temperature dependent)
			self.verifyTable(inpt, False)
			prop = lambda x: self.interpolateScalarProperty(inpt, x)
			print('It will be deprecated')
		else:
			# Anisotropic material (position dependent and temperature dependent)
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
			# Isotropic material (position and temperature independent)
			if np.isscalar(inpt):
				prop = lambda x: inpt*np.eye((*shape, np.max(np.shape(x))))
			else:
				prop = lambda x: create3ArrayFrom2Array(inpt, x)
		elif callable(inpt):
			# Anisotropic material (temperature independent but position dependent)
			prop = lambda x: inpt(x)
		elif len(inpt.shape) > 1:
			# Anisotropic material (position independent but temperature dependent)
			self.verifyTable(inpt, True)
			prop = lambda x: self.interpolateTensorProperty(inpt, x, shape=shape)
			print('It will be deprecated')
		else:
			# Anisotropic material (position dependent and temperature dependent)
			raise Warning('Not implemented')
		return prop

class thermomat(material):
	def __init__(self):
		super().__init__()
		self._capacity     = None
		self._conductivity = None
		self._density      = None
		self._isCapacityIsotropic     = False
		self._isDensityIsotropic      = False
		self._isConductivityIsotropic = False
		return
	
	def addDensity(self, inpt, isIsotropic):
		if isIsotropic: self._isDensityIsotropic = True
		self._density     = super().setScalarProperty(inpt, isIsotropic=isIsotropic)
		return
	
	def addCapacity(self, inpt, isIsotropic):
		if isIsotropic: self._isCapacityIsotropic = True
		self._capacity    = super().setScalarProperty(inpt, isIsotropic=isIsotropic)
		return
	
	def addConductivity(self, inpt, isIsotropic, shape=(3, 3)):
		if isIsotropic: self._isConductivityIsotropic = True
		self._conductivity = super().setTensorProperty(inpt, shape=shape, isIsotropic=isIsotropic)
		return
	
	def eval_capacityCoefficients(self, detJ, inpt): 
		prop = self._capacity(inpt)
		coefs, info = geophy.eval_capacity_coefficient(detJ, prop)
		if info == 0: raise Warning('It is not possible to compute coefficients')
		return coefs
	
	def eval_conductivityCoefficients(self, invJ, detJ, inpt):
		prop = self._conductivity(inpt)
		coefs, info = geophy.eval_conductivity_coefficient(invJ, detJ, prop)
		if info == 0: raise Warning('It is not possible to compute coefficients')
		return coefs
	
	def eval_heatForceCoefficients(self, fun, detJ, qp): 
		coefs = fun(qp)*detJ
		return coefs

class plasticLaw():
	def __init__(self, elasticmodulus, elasticlimit, kwargs:dict):
		self._elasticlimit   = elasticlimit
		self._elasticmodulus = elasticmodulus
		self._Kfun = None; self._Kderfun = None
		self._Hfun = None; self._Hderfun = None
		lawName = kwargs.get('name', 'linear').lower()
		if lawName == 'linear': self._setLinearModel(kwargs)
		if lawName == 'swift' : self._setSwiftModel(kwargs)
		if lawName == 'voce'  : self._setVoceModel(kwargs)
		funlist = [self._Hfun, self._Hderfun, self._Kfun, self._Kderfun]
		if any([fun is None for fun in funlist]): raise Warning('Something went wrong')
		return	
	
	def _setLinearModel(self, kwargs:dict):
		theta	  = kwargs.get('theta', None)
		Hbar      = kwargs.get('Hbar', None)
		self._Kfun = lambda a: self._elasticlimit + theta*Hbar*a 
		self._Hfun = lambda a: (1-theta)*Hbar*a
		self._Kderfun = lambda a: theta*Hbar
		self._Hderfun = lambda a: (1-theta)*Hbar
		return
	
	def _setSwiftModel(self, kwargs:dict):
		K = kwargs.get('K', None)
		n = kwargs.get('exp', None)
		self._Kfun = lambda a: self._elasticlimit + self._elasticmodulus*(a/K)**n
		self._Hfun = lambda a: 0.0
		self._Kderfun = lambda a: (self._elasticmodulus/K)*n*(a/K)**(n-1) if a!=0 else 1e10*self._elasticmodulus
		self._Hderfun = lambda a: 0.0
		return
	
	def _setVoceModel(self, kwargs:dict):
		theta = kwargs.get('theta', None)
		Hbar  = kwargs.get('Hbar', None)
		Kinf  = kwargs.get('Kinf', None)
		delta = kwargs.get('delta', None)
		self._Kfun = lambda a: self._elasticlimit + theta*Hbar*a + Kinf*(1.0 - np.exp(-delta*a))
		self._Hfun = lambda a: (1-theta)*Hbar*a
		self._Kderfun = lambda a: theta*Hbar + Kinf*delta*np.exp(-delta*a) 
		self._Hderfun = lambda a: (1-theta)*Hbar
		return
class mechamat(material):
	# Eventually this class should be similar to thermomat, but for the moment let's say it works !
	def __init__(self, kwargs:dict):
		super().__init__()
		self.getInfo(kwargs)
		return
	
	def _setExtraMechanicalProperties(self):
		E  = self._elasticmodulus
		nu = self._poissonratio
		self._lame_lambda, self._lame_mu, self._lame_bulk = None, None, None
		if E is not None and nu is not None:
			lamb = nu*E/((1+nu)*(1-2*nu))
			mu = E/(2*(1+nu))
			bulk = lamb + 2.0/3.0*mu
			self._lame_lambda = lamb
			self._lame_mu = mu
			self._lame_bulk = bulk
		return
	
	def getInfo(self, kwargs:dict):		
		self._elasticmodulus = kwargs.get('elastic_modulus', None)
		self._poissonratio   = kwargs.get('poisson_ratio', None)
		self._elasticlimit   = kwargs.get('elastic_limit', None)
		self._density        = kwargs.get('density', None)
		plasticKwargs 		 = kwargs.get('law', None)
		self._isPlasticityPossible = False
		if plasticKwargs is not None: 
			self._isPlasticityPossible = True
			self._mechaBehavLaw = plasticLaw(self._elasticmodulus, self._elasticlimit, plasticKwargs)
		self._setExtraMechanicalProperties()
		return
	
	def verifyMechanicalProperties(self):
		" Verifies if mechanical properties exits "
		proplist = [self._elasticmodulus, self._elasticlimit, self._poissonratio]
		if any([prop is None for prop in proplist]): raise Warning('Mechanics not well defined')
		return
	
	def eval_volForceCoefficients(self, fun, detJ, qp):
		qp = np.atleast_2d(qp)
		coefs = fun(qp)*detJ*self._density
		return coefs
	
	def returnMappingAlgorithm(self, law:plasticLaw, strain, pls, a, b):
		""" Return mapping algorithm for multidimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		def computeDeltaGamma(law: plasticLaw, lame_mu, a_n0, eta_trial, nbIter=20, threshold=1e-8):
			dgamma = 0.0
			a_n1 = a_n0
			for i in range(nbIter):
				dH = law._Hfun(a_n1) - law._Hfun(a_n0) 
				G  = (-np.sqrt(2.0/3.0)*law._Kfun(a_n1) + np.linalg.norm(eta_trial, axis=(1, 2)) 
					- (2.0*lame_mu*dgamma + np.sqrt(2.0/3.0)*dH))
				if np.abs(G) <=threshold: break
				dG = - 2.0*(lame_mu + (law._Hderfun(a_n1) + law._Kderfun(a_n1))/3.0)
				dgamma -= G/dG
				a_n1 = a_n0 + np.sqrt(2.0/3.0)*dgamma

			return dgamma

		ddl, nnz = np.shape(strain)
		if ddl == 3: dim = 2
		elif ddl == 6: dim = 3

		Cep     = np.zeros((ddl+3, nnz))
		Tstrain = array2symtensor(strain)
		Tpls    = array2symtensor(pls)
		Tb      = array2symtensor(b)
		traceStrain = evalTrace(Tstrain)
		devStrain   = Tstrain
		for i in range(dim): devStrain[i, i, :] -= 1.0/3.0*traceStrain

		# Compute trial stress
		s_trial = 2*self._lame_mu*(devStrain - Tpls)

		# Compute shifted stress
		eta_trial = s_trial - Tb

		# Check yield condition
		norm_trial = np.linalg.norm(eta_trial, axis=(1, 2))
		f_trial = norm_trial - np.sqrt(2.0/3.0)*law._Kfun(a)
		sigma   = s_trial
		for i in range(dim): sigma[i, i, :] += self._lame_bulk*traceStrain
		Cep[0, :] = self._lame_lambda; Cep[1, :] = self._lame_mu

		if f_trial<=0.0:
			pls_new = pls; a_new = a; b_new = b
			stress  = symtensor2array(sigma)
			return

		# Compute plastic-strain increment
		dgamma = computeDeltaGamma(law, self._lame_mu, a, eta_trial)

		# Update internal hardening variable
		a_new = a + np.sqrt(2.0/3.0)*dgamma

		# Compute df/dsigma
		Normal = eta_trial/norm_trial
		Cep[3:, :] = symtensor2array(Normal)

		# Update stress
		sigma -= 2*self._lame_mu*dgamma*Normal
		stress = symtensor2array(sigma)

		# Update plastic strain
		Tpls += dgamma*Normal
		pls_new = symtensor2array(Tpls)

		# Update backstress
		Tb += np.sqrt(2.0/3.0)*(law._Hfun(a_new) - law._Hfun(a))*Normal
		b_new = symtensor2array(Tb)

		# Update new coefficients
		somme = law._Kderfun(a_new) + law._Hderfun(a_new)
		c1 = 2*self._lame_mu*dgamma/norm_trial
		c2 = 1.0/(1+somme/(3*self._lame_mu)) - c1
		Cep[0, :] = self._lame_lambda + 2.0/3.0*self._lame_mu*c1
		Cep[1, :] = self._lame_mu*(1.0 - c1)
		Cep[2, :] = -2.0*self._lame_mu*c2

		return [stress, pls_new, a_new, b_new, Cep]
	
def symtensor2array(tensor):
	tensor = np.atleast_3d(tensor)
	dim, _, nnz = np.shape(tensor)
	ddl = int(dim*(dim+1)/2)

	array = np.zeros((ddl, nnz))
	k = 0
	for i in range(dim):
		array[k, :] = tensor[i, i, :]
		k += 1

	for i in range(dim-1):
		for j in range(i+1, dim):
			array[k, :] = tensor[i, j, :]

	return array

def array2symtensor(array):
	array = np.atleast_2d(array)
	ddl, nnz = np.shape(array)
	if ddl == 3: dim = 2
	elif ddl == 6: dim = 3
	tensor = np.zeros((dim, dim, nnz))
	k = 0
	for i in range(dim):
		tensor[i, i, :] = array[k, :]
		k += 1

	for i in range(dim-1):
		for j in range(i+1, dim):
			tensor[i, j, :] = array[k, :] 
			tensor[j, i, :] = array[k, :] 

	return tensor

def evalTrace(tensor):
	tensor = np.atleast_3d(tensor)
	dim, _, nnz = np.shape(tensor)
	trace = np.zeros(nnz)
	for i in range(dim):
		trace += tensor[i, i, :]
	return trace

def computeVonMisesStress(tensor):
	tensor = np.atleast_3d(tensor)
	dim, _, nnz = np.shape(tensor)
	trace = evalTrace(tensor)
	dev   = np.copy(tensor)
	for i in range(dim):
		dev[i, i, :] -= 1.0/3.0*trace
	vm = np.linalg.norm(dev, axis=(1, 2))
	vm = np.sqrt(3.0/2.0)*vm
	return vm

def clean_dirichlet(A, dod):
	""" Set to 0 (Dirichlet condition) the values of an array using the indices in each dimension
		A is actually a vector arranged following each dimension [Au, Av, Aw]
	"""
	dim = np.size(A, axis=1)
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