from lib.__init__ import *

class material():

	def __init__(self):
		return		
	
	def __verifyTable(self, table, isTensor=False, threshold=1.e-8):
		table = np.atleast_2d(table)
		x = np.diff(table[:, 0])
		if np.any(x<threshold): raise Warning('Table is not well defined')
		if isTensor and np.size(table, axis=1) <= 2: raise Warning('Table not well defined')
		return
	
	def __interpolateScalarProperty(table, x):
		lenx = np.max(np.shape(x))
		if np.size(table, axis=0) == 1: y = table[:, 1]*np.ones(lenx)
		else:  							y = np.interp(x, table[:, 0], table[:, 1])
		return y
	
	def __interpolateTensorProperty(table, x, shape=(3, 3)):
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
			# Isotropic material 
			if np.isscalar(inpt): 		
				prop = lambda x: inpt*np.ones(np.max(np.shape(x)))
			else: 
				raise Warning('Not possible')
		elif callable(inpt):
			# Anisotropic material a continuous function
			prop = lambda x: inpt(x)
		elif len(inpt.shape) > 1:
			# Anisotropic material a discrete table
			self.__verifyTable(inpt, False)
			prop = lambda x: self.__interpolateScalarProperty(inpt, x)
			print('It will be deprecated')
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
		elif len(inpt.shape) > 1:
			# Anisotropic material using a discrete table
			self.__verifyTable(inpt, True)
			prop = lambda x: self.__interpolateTensorProperty(inpt, x, shape=shape)
			print('It will be deprecated')
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
	
	def eval_capacityCoefficients(self, detJ, inpt): 
		prop = self.capacity(inpt)
		coefs, info = geophy.eval_capacity_coefficient(detJ, prop)
		if info == 0: raise Warning('It is not possible to compute coefficients')
		return coefs
	
	def eval_conductivityCoefficients(self, invJ, detJ, inpt):
		prop = self.conductivity(inpt)
		coefs, info = geophy.eval_conductivity_coefficient(invJ, detJ, prop)
		if info == 0: raise Warning('It is not possible to compute coefficients')
		return coefs
	
	def eval_heatForceCoefficients(self, fun, detJ, qp): 
		coefs = fun(qp)*detJ
		return coefs

class plasticLaw():
	def __init__(self, elasticmodulus, elasticlimit, plasticVars:dict):
		self.__elasticlimit   = elasticlimit
		self.__elasticmodulus = elasticmodulus
		self._Kfun = None; self._Kderfun = None
		self._Hfun = None; self._Hderfun = None
		lawName = plasticVars.get('name', '').lower()
		if   lawName == 'linear': self.__setLinearModel(plasticVars)
		elif lawName == 'swift' : self.__setSwiftModel(plasticVars)
		elif lawName == 'voce'  : self.__setVoceModel(plasticVars)
		else: raise Warning('Unknown method')
		funlist = [self._Hfun, self._Hderfun, self._Kfun, self._Kderfun]
		if any(fun is None for fun in funlist): raise Warning('Something went wrong')
		return	
	
	def __setLinearModel(self, kwargs:dict):
		theta	  = kwargs.get('theta', None)
		Hbar      = kwargs.get('Hbar', None)
		self._Kfun = lambda a: self.__elasticlimit + theta*Hbar*a 
		self._Hfun = lambda a: (1 - theta)*Hbar*a
		self._Kderfun = lambda a: theta*Hbar
		self._Hderfun = lambda a: (1 - theta)*Hbar
		return
	
	def __setSwiftModel(self, kwargs:dict):
		K = kwargs.get('K', None)
		n = kwargs.get('exp', None)
		self._Kfun = lambda a: self.__elasticlimit + self.__elasticmodulus*(a/K)**n
		self._Hfun = lambda a: 0.0
		self._Kderfun = lambda a: (self.__elasticmodulus/K)*n*(a/K)**(n-1) if a!=0 else 1e10*self.__elasticmodulus
		self._Hderfun = lambda a: 0.0
		return
	
	def __setVoceModel(self, kwargs:dict):
		theta = kwargs.get('theta', None)
		Hbar  = kwargs.get('Hbar', None)
		Kinf  = kwargs.get('Kinf', None)
		delta = kwargs.get('delta', None)
		self._Kfun = lambda a: self.__elasticlimit + theta*Hbar*a + Kinf*(1.0 - np.exp(-delta*a))
		self._Hfun = lambda a: (1 - theta)*Hbar*a
		self._Kderfun = lambda a: theta*Hbar + Kinf*delta*np.exp(-delta*a) 
		self._Hderfun = lambda a: (1 - theta)*Hbar
		return
	
class mechamat(material):
	# Eventually this class should be similar to thermomat, but for the moment let's say it works !
	def __init__(self, kwargs:dict):
		super().__init__()
		self.density        = kwargs.get('density', None)
		self.elasticmodulus = kwargs.get('elastic_modulus', None)
		self.poissonratio   = kwargs.get('poisson_ratio', None)
		self.elasticlimit   = kwargs.get('elastic_limit', None)
		if any(prop is None for prop in [self.elasticmodulus, self.elasticlimit, self.poissonratio]): 
			raise Warning('Mechanics not well defined')

		self.plasticLaw            = None
		self._isPlasticityPossible = False
		tmp = kwargs.get('plasticlaw', None)
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
		
	def eval_volForceCoefficients(self, fun, detJ, qp):
		qp = np.atleast_2d(qp)
		coefs = fun(qp)*detJ*self.density
		return coefs
	
	def returnMappingAlgorithm(self, strain, pls, a, b, threshold=1e-9):
		""" Return mapping algorithm for multidimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		def computeDeltaGamma(law:plasticLaw, lame_mu, a_n0, eta_trial, nbIter=20, threshold=1e-9):
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

		nvoigt = np.size(strain)
		if nvoigt   == 3: dim = 2
		elif nvoigt == 6: dim = 3

		Cep     = np.zeros(nvoigt+3)
		Tstrain = array2symtensor(strain, dim)
		Tpls    = array2symtensor(pls, dim)
		Tb      = array2symtensor(b, dim)
		traceStrain = evalTrace(Tstrain, dim)
		devStrain   = Tstrain
		for i in range(dim): devStrain[i, i] -= 1.0/3.0*traceStrain

		# Compute trial stress
		s_trial = 2*self.lame_mu*(devStrain - Tpls)

		# Compute shifted stress
		eta_trial = s_trial - Tb

		# Check yield condition
		norm_trial = np.linalg.norm(eta_trial)
		f_trial = norm_trial - np.sqrt(2.0/3.0)*self.plasticLaw._Kfun(a)
		sigma   = s_trial
		for i in range(dim): sigma[i, i] += self.lame_bulk*traceStrain
		Cep[0] = self.lame_lambda; Cep[1] = self.lame_mu
		pls_new = pls; a_new = a; b_new = b
		stress  = symtensor2array(sigma, dim)

		if f_trial>threshold:

			# Compute plastic-strain increment
			dgamma = computeDeltaGamma(self.plasticLaw, self.lame_mu, a, eta_trial)

			# Update internal hardening variable
			a_new = a + np.sqrt(2.0/3.0)*dgamma

			# Compute df/dsigma
			Normal = eta_trial/norm_trial
			Cep[3:] = symtensor2array(Normal, dim)

			# Update stress
			sigma -= 2*self.lame_mu*dgamma*Normal
			stress = symtensor2array(sigma, dim)

			# Update plastic strain
			Tpls += dgamma*Normal
			pls_new = symtensor2array(Tpls, dim)

			# Update backstress
			Tb += np.sqrt(2.0/3.0)*(self.plasticLaw._Hfun(a_new) - self.plasticLaw._Hfun(a))*Normal
			b_new = symtensor2array(Tb, dim)

			# Update new coefficients
			somme = self.plasticLaw._Kderfun(a_new) + self.plasticLaw._Hderfun(a_new)
			c1 = 2*self.lame_mu*dgamma/norm_trial
			c2 = 1.0/(1+somme/(3*self.lame_mu)) - c1
			Cep[0] = self.lame_lambda + 2.0/3.0*self.lame_mu*c1
			Cep[1] = self.lame_mu*(1.0 - c1)
			Cep[2] = -2.0*self.lame_mu*c2

		return stress, pls_new, a_new, b_new, Cep
	
def symtensor2array(tensor, dim):
	nvoigt = int(dim*(dim+1)/2)
	array  = np.zeros(nvoigt)
	k = 0
	for i in range(dim):
		array[k] = tensor[i, i]
		k += 1

	for i in range(dim-1):
		for j in range(i+1, dim):
			array[k] = tensor[i, j]
			k += 1

	return array

def array2symtensor(array, dim):
	tensor = np.zeros((dim, dim))
	k = 0
	for i in range(dim):
		tensor[i, i] = array[k]
		k += 1

	for i in range(dim-1):
		for j in range(i+1, dim):
			tensor[i, j] = array[k] 
			tensor[j, i] = array[k] 
			k += 1

	return tensor

def evalTrace(tensor, dim):
	trace = 0.0
	for i in range(dim):
		trace += tensor[i, i]
	return trace

def computeVonMisesStress(tensor, dim):
	trace  = evalTrace(tensor, dim)
	dev    = np.copy(tensor)
	for i in range(dim):
		dev[i, i] -= 1.0/3.0*trace
	vm = np.linalg.norm(dev)
	vm = np.sqrt(3.0/2.0)*vm
	return vm

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