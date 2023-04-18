from lib.__init__ import *
	
class thermomat():
	def __init__(self, **kwargs):
		self._kwargs = kwargs
		self.getInfo()
		return
	
	def getInfo(self):
		kwargs = self._kwargs
		self._heatconductivity = kwargs.get('conductivity', None)
		self._heatcapacity     = kwargs.get('capacity', None)
		self._density          = kwargs.get('density', None)
		return
	
	def verifyThermalProperties(self):
		" Verifies if thermal properties exits "
		if self._heatconductivity is None: raise Warning('Conductivity not defined')
		if self._heatcapacity is None: raise Warning('Capacity not defined')
		return
	
	def eval_conductivityCoefficients(self, invJ, detJ, prop):
		prop = np.atleast_3d(prop)
		coefs = np.zeros(np.shape(invJ))
		coefs, info = assembly.eval_conductivity_coefficient(invJ, detJ, prop)
		if info == 0: raise Warning('It is not possible to compute coefficients')
		return coefs
	
	def eval_capacityCoefficients(self, detJ, prop): 
		prop = np.atleast_1d(prop)
		coefs = np.zeros(np.shape(detJ))
		coefs, info = assembly.eval_capacity_coefficient(detJ, prop)
		if info == 0: raise Warning('It is not possible to compute coefficients')
		return coefs
	
	def eval_heatForceCoefficients(self, fun, detJ, qp): 
		qp = np.atleast_2d(qp)
		coefs = fun(qp)*detJ
		return coefs
	
class mechamat():
	def __init__(self, **kwargs):
		self._kwargs = kwargs
		self.getInfo()
		return
	
	def __setExtraMechanicalProperties(self):
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
	
	def getInfo(self):
		kwargs = self._kwargs
		self._elasticmodulus = kwargs.get('elastic_modulus', None)
		self._poissonratio   = kwargs.get('poisson_ratio', None)
		self._elasticlimit   = kwargs.get('elastic_limit', None)
		self._density        = kwargs.get('density', None)
		self.__setExtraMechanicalProperties()
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