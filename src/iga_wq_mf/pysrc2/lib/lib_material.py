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