from .__init__ import *

class material():

	def __init__(self):
		self.density = None
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
	
	def addDensity(self, inpt, isIsotropic):
		self.density = self.setScalarProperty(inpt, isIsotropic=isIsotropic)
		return

class heatmat(material):
	def __init__(self, matArgs=dict()):
		super().__init__()
		self.capacity     = None
		self.conductivity = None
		self._isCapacityIsotropic     = False
		self._isConductivityIsotropic = False
		self.refTemperature = matArgs.get('RefTemp', 300)
		return
	
	def addCapacity(self, inpt, isIsotropic):
		if isIsotropic: self._isCapacityIsotropic = True
		self.capacity = super().setScalarProperty(inpt, isIsotropic=isIsotropic)
		return
	
	def addConductivity(self, inpt, isIsotropic, shape=(3, 3)):
		if isIsotropic: self._isConductivityIsotropic = True
		self.conductivity = super().setTensorProperty(inpt, shape=shape, isIsotropic=isIsotropic)
		return

class plasticLaw():
	def __init__(self, elasticlimit, plasticArgs:dict):
		self.__elasticlimit = elasticlimit
		self._IsotropicHard = None; self._IsotropicHardDer = None
		self._KinematicHard = None; self._KinematicHardDer = None
		self._plasticArgs = plasticArgs

		# Define Isotropic hardening
		isoname = plasticArgs.get('Isoname', None)
		if   isoname == 'linear': self.__setIsoLinearModel(plasticArgs)
		elif isoname == 'swift' : self.__setIsoSwiftModel(plasticArgs)
		elif isoname == 'voce'  : self.__setIsoVoceModel(plasticArgs)
		elif isoname == 'none' : self.__setIsoNoneModel()
		else: raise Warning('Unknown method')

		# Define kinematic hardening
		kinename = plasticArgs.get('Kinename', None)
		if   kinename == 'linear': self.__setKineLinearModel(plasticArgs)
		else: 
			print('By default, we do not consider kinematic hardening')
			self.__setKineNoneModel()

		funlist = [self._KinematicHard, self._KinematicHardDer, self._IsotropicHard, self._IsotropicHardDer]
		if any(fun is None for fun in funlist): raise Warning('Something went wrong')
		return	
	
	def __setIsoNoneModel(self):
		self._IsotropicHard = lambda a: 1e6*self.__elasticlimit*np.ones(np.size(a))
		self._IsotropicHardDer = lambda a: np.zeros(np.size(a))
		return
	
	def __setIsoLinearModel(self, plasticArgs:dict):
		Eiso = plasticArgs.get('Eiso', None)
		self._IsotropicHard = lambda a: self.__elasticlimit + Eiso*a 
		self._IsotropicHardDer = lambda a: Eiso*np.ones(np.size(a))
		return
	
	def __setIsoSwiftModel(self, plasticArgs:dict):
		e0 = plasticArgs.get('e0', None)
		n = plasticArgs.get('n', None)
		self._IsotropicHard = lambda a: self.__elasticlimit*(1 + a/e0)**n
		self._IsotropicHardDer = lambda a: self.__elasticlimit*n/e0*(1 + a/e0)**(n - 1)
		return
	
	def __setIsoVoceModel(self, plasticArgs:dict):
		ssat = plasticArgs.get('ssat', None)
		beta = plasticArgs.get('beta', None)
		self._IsotropicHard = lambda a: self.__elasticlimit +  ssat*(1.0 - np.exp(-beta*a))
		self._IsotropicHardDer = lambda a: ssat*beta*np.exp(-beta*a) 
		return
	
	def __setKineNoneModel(self):
		self._KinematicHard = lambda a: np.zeros(np.size(a))
		self._KinematicHardDer = lambda a: np.zeros(np.size(a))
		return
	
	def __setKineLinearModel(self, plasticArgs:dict):
		Ek = plasticArgs.get('Ekine', None)
		print('If working in 2 or 3 dimensions, please be sure of scailing the coefficients by 2/3')
		self._KinematicHard = lambda a: Ek*np.ones(np.size(a))
		self._KinematicHardDer = lambda a: np.zeros(np.size(a))
		return
	
class mechamat(material):
	# Eventually this class should be similar to thermomat, but for the moment let's say it works !
	def __init__(self, matArgs:dict):
		super().__init__()
		self.elasticmodulus = matArgs.get('elastic_modulus', None)
		self.poissonratio   = matArgs.get('poisson_ratio', None)
		self.elasticlimit   = matArgs.get('elastic_limit', None)
		self.thexpansion    = matArgs.get('thermal_expansion', 1.0)
		if any(prop is None for prop in [self.elasticmodulus, self.elasticlimit, self.poissonratio]): 
			raise Warning('Mechanics not well defined')

		self.plasticLaw            = None
		self._isPlasticityPossible = False
		tmp = matArgs.get('plasticLaw', None)
		if isinstance(tmp, dict): 
			self._isPlasticityPossible = True
			self.plasticLaw = plasticLaw(self.elasticlimit, tmp)
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
		
	def evalElasticStress(self, strain, dim):
		Tstrain = array2symtensor4All(strain, dim)
		traceStrain = evalTrace4All(strain, dim)
		Tstress = 2*self.lame_mu*Tstrain
		for i in range(dim): Tstress[i, i, :] += self.lame_lambda*traceStrain
		stress  = symtensor2array4All(Tstress, dim)
		return stress
	
	def returnMappingAlgorithm(self, strain, pls, a, b, threshold=1e-9):
		""" Return mapping algorithm for multidimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		def computeDeltaGamma(law:plasticLaw, lame_mu, a_n0, eta_trial, nbIter=50, threshold=1e-9):
			dgamma = np.zeros(np.size(a_n0))
			a_n1 = a_n0
			for i in range(nbIter):
				G  = (np.linalg.norm(eta_trial, axis=(0, 1)) - (2.0*lame_mu + law._KinematicHard(a_n1))*dgamma
				- np.sqrt(2.0/3.0)*law._IsotropicHard(a_n1))
				if np.all(np.abs(G)<=threshold): break
				dG = -(2.0*lame_mu + law._KinematicHard(a_n1) 
					+ np.sqrt(2.0/3.0)*dgamma*law._KinematicHardDer(a_n1) 
					+ 2.0/3.0*law._IsotropicHardDer(a_n1))
				dgamma -= G/dG
				a_n1 = a_n0 + np.sqrt(2.0/3.0)*dgamma
			return dgamma

		nvoigt, nnz = np.shape(strain)
		if nvoigt   == 3: dim = 2
		elif nvoigt == 6: dim = 3
		isElasticLoad = True

		output  = np.zeros((4*(nvoigt+1), nnz)); Cep = np.zeros((3, nnz))
		Tstrain = array2symtensor4All(strain, dim)
		Tpls    = array2symtensor4All(pls, dim)
		Tb      = array2symtensor4All(b, dim)

		# Compute strain deviator
		traceStrain = evalTrace4All(strain, dim)
		devStrain   = np.copy(Tstrain)
		for i in range(dim): devStrain[i, i, :] -= 1.0/3.0*traceStrain

		# Compute trial stress deviator
		s_trial = 2*self.lame_mu*(devStrain - Tpls)

		# Compute shifted stress
		eta_trial = s_trial - Tb

		# Check yield status
		norm_trial = np.linalg.norm(eta_trial, axis=(0, 1))
		f_trial = norm_trial - np.sqrt(2.0/3.0)*self.plasticLaw._IsotropicHard(a)
		sigma   = np.copy(s_trial)
		for i in range(dim): sigma[i, i, :] += self.lame_bulk*traceStrain
		Cep[0, :] = self.lame_lambda; Cep[1, :] = self.lame_mu
		pls_new = np.copy(pls); a_new = np.copy(a); b_new = np.copy(b)
		stress  = symtensor2array4All(sigma, dim)

		plsInd = np.nonzero(f_trial>threshold)[0]
		if np.size(plsInd) > 0:

			isElasticLoad = False

			# Compute plastic-strain increment
			dgamma_plsInd = computeDeltaGamma(self.plasticLaw, self.lame_mu, a[plsInd], eta_trial[:, :, plsInd])

			# Update internal hardening variable
			a_new[plsInd] = a[plsInd] + np.sqrt(2.0/3.0)*dgamma_plsInd

			# Compute d f_trial/d eta_trial
			Normal_plsInd = eta_trial[:, :, plsInd]/norm_trial[plsInd]
			output[3*nvoigt+4:, plsInd] = symtensor2array4All(Normal_plsInd, dim)

			# Update stress
			sigma[:, :, plsInd] -= 2*self.lame_mu*dgamma_plsInd*Normal_plsInd
			stress[:, plsInd] = symtensor2array4All(sigma[:, :, plsInd], dim)

			# Update plastic strain
			Tpls[:, :, plsInd] += dgamma_plsInd*Normal_plsInd
			pls_new[:, plsInd] = symtensor2array4All(Tpls[:, :, plsInd], dim)

			# Update backstress
			Tb[:, :, plsInd] += self.plasticLaw._KinematicHard(a_new[plsInd])*dgamma_plsInd*Normal_plsInd
			b_new[:, plsInd] = symtensor2array4All(Tb[:, :, plsInd], dim)

			# Update tangent coefficients
			sumofterms = (3.0*self.plasticLaw._KinematicHard(a_new[plsInd]) 
						+ 2.0*self.plasticLaw._IsotropicHardDer(a_new[plsInd]) 
						+ np.sqrt(6.0)*dgamma_plsInd*self.plasticLaw._KinematicHardDer(a_new[plsInd]))
			c1 = 1.0 - 2.0*self.lame_mu*dgamma_plsInd/norm_trial[plsInd]
			c2 = 1.0/(1.0 + sumofterms/(6.0*self.lame_mu)) + c1 - 1.0
			Cep[0, plsInd] = self.lame_lambda + 2.0/3.0*self.lame_mu*(1.0 - c1)
			Cep[1, plsInd] = self.lame_mu*c1
			Cep[2, plsInd] = -2.0*self.lame_mu*c2

		output[0:nvoigt, :] = stress; output[nvoigt:2*nvoigt, :] = pls_new; output[2*nvoigt, :] = a_new
		output[2*nvoigt+1:3*nvoigt+1, :] = b_new; output[3*nvoigt+1:3*nvoigt+4, :] = Cep
		return output, isElasticLoad
	
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