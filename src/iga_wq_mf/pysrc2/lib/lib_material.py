from lib.__init__ import *
	
class thermomat():
	def __init__(self, **kwargs):
		self._density = None
		self._heatconductivity = None
		self._heatcapacity     = None
		self.getInfo(kwargs)
		return
	
	def __setConductivity(self, fun):
		self._heatconductivity = fun
		return
	
	def __setCapacity(self, fun):
		self._heatcapacity = fun
		return
	
	def __setDensity(self, val):
		self._density = val
		return

	def getInfo(self, kwargs:dict):
		conductivity = kwargs.get('conductivity', None)
		capacity     = kwargs.get('capacity', None)
		density      = kwargs.get('density', None)
		self.__setCapacity(capacity)
		self.__setConductivity(conductivity)
		self.__setDensity(density)
		return
	
class mechamat():
	def __init__(self, **kwargs):
		self._density = None
		self._elasticmodulus = None
		self._poissonratio   = None
		self._elasticlimit   = None
		self.getInfo(kwargs)
		return
	
	def __setElasticModulus(self, val):
		self._elasticmodulus = val
		return
	
	def __setPoissonRatio(self, val):
		self._poissonratio = val
		return
	
	def __setElasticLimit(self, val):
		self._elasticlimit = val
		return
	
	def __setDensity(self, val):
		self.__setDensity(val)
		return
	
	def getInfo(self, kwargs:dict):
		elasticmodulus = kwargs.get('elastic_modulus', None)
		poissonratio   = kwargs.get('poisson_ratio', None)
		elasticlimit   = kwargs.get('elastic_limit', None)
		density        = kwargs.get('density', None)
		self.__setElasticModulus(elasticmodulus)
		self.__setPoissonRatio(poissonratio)
		self.__setElasticLimit(elasticlimit)
		self.__setDensity(density)
		return

class mechamat1D(mechamat):
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.Hfun, self.Hderfun = None, None
		self.Kfun, self.Kderfun = None, None
		name = kwargs.get('model', 'linear')
		if name.lower() == 'linear': self.__setLinearModel(kwargs)
		if name.lower() == 'swift' : self.__setSwiftModel(kwargs)
		if name.lower() == 'voce'  : self.__setVoceModel(kwargs)
		funlist = [self.Hfun, self.Hderfun, self.Kfun, self.Kderfun]
		if any([fun is None for fun in funlist]): raise Warning('Something went wrong')
		return
	
	def __setLinearModel(self, kwargs:dict):
		theta	  = kwargs.get('theta', None)
		Hbar      = kwargs.get('Hbar', None)
		self.Kfun = lambda a: self._elasticlimit + theta*Hbar*a 
		self.Hfun = lambda a: (1-theta)*Hbar*a
		self.Kderfun = lambda a: theta*Hbar
		self.Hderfun = lambda a: (1-theta)*Hbar
		return
	
	def __setSwiftModel(self, kwargs:dict):
		K = kwargs.get('K', None)
		n = kwargs.get('exp', None)
		self.Kfun = lambda a: self._elasticlimit + self._elasticmodulus*(a/K)**n
		self.Hfun = lambda a: 0.0
		self.Kderfun = lambda a: (self._elasticmodulus/K)*n*(a/K)**(n-1) if a!=0 else 1e10*self._elasticmodulus
		self.Hderfun = lambda a: 0.0
		return
	
	def __setVoceModel(self, kwargs:dict):
		theta = kwargs.get('theta', None)
		Hbar  = kwargs.get('Hbar', None)
		Kinf  = kwargs.get('Kinf', None)
		delta = kwargs.get('delta', None)
		self.Kfun = lambda a: self._elasticlimit + theta*Hbar*a + Kinf*(1.0 - np.exp(-delta*a))
		self.Hfun = lambda a: (1-theta)*Hbar*a
		self.Kderfun = lambda a: theta*Hbar + Kinf*delta*np.exp(-delta*a) 
		self.Hderfun = lambda a: (1-theta)*Hbar
		return

	def returnMappingAlgorithm(self, strain, pls, a, b, threshold=1e-8):
		""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		def computeDeltaGamma(self, eta_trial, a_n0, nbIter=20, threshold=1e-8):
			dgamma = 0.0
			a_n1   = a_n0 
			for i in range(nbIter):
				dH = self.Hfun(a_n1) - self.Hfun(a_n0) 
				G  = -self.Kfun(a_n1) + np.abs(eta_trial) - (self._elasticmodulus*dgamma + dH)
				if G <=threshold: break
				dG = - (self._elasticmodulus + self.Hderfun(a_n1) + self.Kderfun(a_n1))
				dgamma = dgamma - G/dG
				a_n1   = a_n0 + dgamma 
			return dgamma

		# Elastic predictor
		sigma_trial = self._elasticmodulus*(strain - pls)
		eta_trial   = sigma_trial - b

		# Check yield status
		f_trial = np.abs(eta_trial) - self.Kfun(a)

		if f_trial <= threshold: # Elastic
			stress = sigma_trial
			pls_new = pls
			a_new = a
			b_new = b
			Cep = self._elasticmodulus

		else: # Plastic
			N = np.sign(eta_trial)
			dgamma = computeDeltaGamma(eta_trial, a)
			stress = sigma_trial - dgamma*self._elasticmodulus*N
			pls_new = pls + dgamma*N
			a_new = a + dgamma
			b_new = b + (self.Hfun(a_new)-self.Hfun(a))*N
			somme = self.Kderfun(a_new) + self.Hderfun(a_new)
			Cep = self._elasticmodulus*somme/(self._elasticmodulus + somme)

		return [stress, pls_new, a_new, b_new, Cep]
	
class material:
	def __init__(self):
		self._thermomat = None
		self._mechamat  = None
		return
	
	def addThermalMaterial(self, mat:thermomat):
		self._thermomat = mat
		return
	
	def addMechanicalMaterial(self, mat:mechamat):
		self._mechamat = mat
		return
	
	def eraseMaterial(self):
		self._thermomat = None
		self._mechamat  = None
		return