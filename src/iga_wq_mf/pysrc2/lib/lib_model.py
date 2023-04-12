"""
.. This module helps to read information from a given geometry (YETI format)
.. and to create thermo-mechanical model (material, boundary conditions, etc.)
.. Joaquin Cornejo
"""

# My libraries
from .lib_quadrules import *
from .lib_material import *

class part(): 

	def __init__(self, modelIGA, **kwargs):
		print('\nInitializing thermo-mechanical model')
		self._sampleSize = kwargs.get('sample_size', 101)
		self._name       = self.read_name(modelIGA)
		self._dim        = self.read_dimension(modelIGA)
		self._degree     = self.read_degree(modelIGA)
		self._ctrlpts    = self.read_controlPoints(modelIGA)
		self._kwargs     = kwargs
		self._knotvector, self._size_kv = self.read_knotvector(modelIGA)

		self._nbqp, self._nbqp_total = [], None
		self._qpPar, self._basis, self._weights, self._indices = [], [], [], []
		self.__setQuadratureRules()

		self._Jqp, self._detJ, self._invJ, self._qpPhy = None, None, None, None
		self.__setJacobienPhysicalPoints()

		return
	
	# ----------
	# READ FILE
	# ----------
	
	def read_name(self, modelIGA: IGAparametrization): 
		" Reads name from model "
		try: name = modelIGA._name
		except: name = 'IGAparametrization'
		return name

	def read_degree(self, modelIGA: IGAparametrization): 
		" Reads degree from model "
		degree = modelIGA._Jpqr.flatten()
		if any(p == 1 for p in degree[:self._dim]): 
			raise Warning('Model must have at least degree p = 2')
		return degree

	def read_dimension(self, modelIGA: IGAparametrization):
		" Reads dimensions from model "
		dim = modelIGA._dim[0]
		if dim != 3: print("WARNING: Some functions have not been well-implemented for 2D geometries")
		return dim

	def read_knotvector(self, modelIGA: IGAparametrization):
		" Reads knot-vector from model "
		knotvector = modelIGA._Ukv[0]
		size_kv = modelIGA._Nkv.flatten()
		return knotvector, size_kv

	def read_controlPoints(self, modelIGA: IGAparametrization): 
		" Reads control points from model "
		ctrlpts = modelIGA._COORDS[:self._dim, :]
		return ctrlpts
	
	# -----------

	def __setQuadratureRules(self):

		def eval_basisWeights(quadRule: WeightedQuadrature):
			quadPtsPos, dersIndices, dersBasis, dersWeights = quadRule.getQuadratureRulesInfo()
			nbqp = len(quadPtsPos)
			return nbqp, quadPtsPos, dersIndices, dersBasis, dersWeights
		
		kwargs        = self._kwargs
		quadRuleName  = kwargs.get('quadrule', 'wq').lower()
		quadRule_list = []
		
		print('Evaluating basis and weights')
		start = time.process_time()
		if quadRuleName == 'iga':
			for i in range(self._dim):
				quadRule_list.append(GaussQuadrature(self._degree[i], self._knotvector[i], **kwargs))
		elif quadRuleName == 'wq':
			for i in range(self._dim):
				quadRule_list.append(WeightedQuadrature(self._degree[i], self._knotvector[i], **kwargs))
		else:
			raise Warning('Not found')
		
		for quadRule in quadRule_list:
			nbqp, quadPtsPos, dersIndices, dersBasis, dersWeights = eval_basisWeights(quadRule)
			self._nbqp.append(nbqp); self._qpPar.append(quadPtsPos); self._indices.append(dersIndices)
			self._basis.append(dersBasis), self._weights(dersWeights)
		stop = time.process_time()
		print('\tBasis and weights in : %.5f s' %(stop-start))

		# Update number of quadrature points
		self._nbqp_total = np.prod(self._nbqp)

		return
	
	def __setJacobienPhysicalPoints(self):
		" Computes jacobien and physical position "

		print('Evaluating jacobien and physical position')
		start = time.process_time()
		
		inputs = [*self._nbqp, *self._indices, *self._basis, self._ctrlpts]
		if self._dim == 2:
			self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_2d(*inputs)
			self._qp_PS = assembly.interpolate_fieldphy_2d(*inputs)
		if self._dim == 3:
			self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_3d(*inputs)
			self._qp_PS = assembly.interpolate_fieldphy_3d(*inputs)
		stop = time.process_time()
		print('\t Time jacobien: %.5f s' %(stop-start))
		return
	
	# ----------------
	# POST-PROCESSING 
	# ----------------