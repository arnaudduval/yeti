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
		y = np.interp(x, table[:, 0], table[:, 1])
		return y
	
	def interpolateTensorProperty(table, x, shape=(3, 3)):
		if np.size(table, axis=1) != np.prod(shape)+1: raise Warning('Not possible')
		y = np.zeros((*shape, len(x)))
		for i in range(shape[0]):
			for j in range(shape[1]):
				k = 1 + j + shape[1]*i
				y[i, j, :] = np.interp(x, table[:, 0], table[:, k])
		return y
	
	def setScalarProperty(self, input):
		if np.isscalar(input):
			# Isotropic material (position and temperature independent)
			prop = lambda x: input*np.ones(len(x))
		elif callable(input):
			# Anisotropic material (temperature independent but position dependent)
			prop = lambda x: input(x)
		elif len(input.shape) > 1:
			# Anisotropic material (position independent but temperature dependent)
			self.verifyTable(input, False)
			prop = lambda x: self.interpolateScalarProperty(input, x)
		else:
			# Anisotropic material (position dependent and temperature dependent)
			raise Warning('Not implemented')
		return prop
	
	def setTensorProperty(self, input, shape=(3, 3)):
		if np.isscalar(input):
			# Isotropic material (position and temperature independent)
			prop = lambda x: input*np.ones((*shape, len(x)))
		elif callable(input):
			# Anisotropic material (temperature independent but position dependent)
			prop = lambda x: input(x)
		elif len(input.shape) > 1:
			# Anisotropic material (position independent but temperature dependent)
			self.verifyTable(input, True)
			prop = lambda x: self.interpolateScalarProperty(input, x, shape=shape)
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
		return
	
	def addDensity(self, input):
		self._density     = super().setScalarProperty(input)
		return
	
	def addCapacity(self, input):
		self._capacity    = super().setScalarProperty(input)
		return
	
	def addConductivity(self, input, shape=(3, 3)):
		self._conductivity = super().setTensorProperty(input, shape=shape)
		return
	
	def eval_capacityCoefficients(self, detJ, input): 
		prop = self._capacity(input)
		coefs, info = assembly.eval_capacity_coefficient(detJ, prop)
		if info == 0: raise Warning('It is not possible to compute coefficients')
		return coefs
	
	def eval_conductivityCoefficients(self, invJ, detJ, input):
		prop = self._conductivity(input)
		coefs, info = assembly.eval_conductivity_coefficient(invJ, detJ, prop)
		if info == 0: raise Warning('It is not possible to compute coefficients')
		return coefs
	
	def eval_heatForceCoefficients(self, fun, detJ, qp): 
		coefs = fun(qp)*detJ
		return coefs
	
class mechamat(material):
	# Eventually this class should be similar to thermomat, but for the moment let's say it works !
	def __init__(self, **kwargs):
		super().__init__()
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