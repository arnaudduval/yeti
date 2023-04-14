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
	
	def verifyThermalProperties(self):
		" Verifies if thermal properties exits "
		if self._heatconductivity is None: raise Warning('Conductivity not defined')
		if self._heatcapacity is None: raise Warning('Capacity not defined')
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
	
	def getInfo(self, kwargs:dict):
		elasticmodulus = kwargs.get('elastic_modulus', None)
		poissonratio   = kwargs.get('poisson_ratio', None)
		elasticlimit   = kwargs.get('elastic_limit', None)
		density        = kwargs.get('density', None)
		self.__setElasticModulus(elasticmodulus)
		self.__setPoissonRatio(poissonratio)
		self.__setElasticLimit(elasticlimit)
		self.__setDensity(density)
		self.__setExtraMechanicalProperties()
		return
	
	def verifyMechanicalProperties(self):
		" Verifies if mechanical properties exits "
		proplist = [self._elasticmodulus, self._elasticlimit, self._poissonratio]
		if any([prop is None for prop in proplist]): raise Warning('Mechanics not well defined')
		
		return
	
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