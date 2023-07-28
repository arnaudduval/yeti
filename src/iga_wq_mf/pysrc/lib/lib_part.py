"""
.. This module helps to read information from a given geometry (YETI format)
.. and to create thermo-mechanical model (material, boundary conditions, etc.)
.. Joaquin Cornejo
"""

# My libraries
from lib.__init__ import *
from lib.lib_base import evalDersBasisFortran
from lib.lib_quadrules import WeightedQuadrature, GaussQuadrature

class part(): 

	def __init__(self, modelIGA, quadArgs:dict):
		print('\nInitializing thermo-mechanical model')
		self.name       = self.__read_name(modelIGA)
		self.dim        = self.__read_dimension(modelIGA)
		self.degree     = self.__read_degree(modelIGA)
		self.ctrlpts    = self.__read_controlPoints(modelIGA)
		self.nbctrlpts  = np.ones(3, dtype=int)
		self.knotvector, self.size_kv = self.__read_knotvector(modelIGA)
		for i in range(self.dim): self.nbctrlpts[i] = len(self.knotvector[i]) - self.degree[i] - 1
		self.nbctrlpts_total = np.product(self.nbctrlpts)

		self.nbqp, self.nbqp_total = np.ones(3, dtype=int), None
		self.basis, self.weights, self.indices = [], [], []
		self.__setQuadratureRules(quadArgs)

		self.Jqp, self.detJ, self.invJ, self.qpPhy = None, None, None, None
		self.__setJacobienPhysicalPoints()
		if np.any(self.detJ<0.0): raise Warning('Geometry problem. See control points positions')

		return
	
	# ----------
	# READ FILE
	# ----------
	
	def __read_name(self, modelIGA:IGAparametrization): 
		" Reads name from model "
		try: name = modelIGA._name
		except: name = 'IGAparametrization'
		return name

	def __read_degree(self, modelIGA:IGAparametrization): 
		" Reads degree from model "
		degree = modelIGA._Jpqr.flatten()
		if any(p == 1 for p in degree[:self.dim]): 
			raise Warning('Model must have at least degree p = 2')
		return degree

	def __read_dimension(self, modelIGA:IGAparametrization):
		" Reads dimensions from model "
		dim = modelIGA._dim[0]
		if dim != 3: print("WARNING: Some functions may have not been well-implemented for 2D geometries")
		return dim

	def __read_knotvector(self, modelIGA:IGAparametrization):
		" Reads knot-vector from model "
		knotvector = modelIGA._Ukv[0]
		size_kv = modelIGA._Nkv.flatten()
		return knotvector, size_kv

	def __read_controlPoints(self, modelIGA:IGAparametrization): 
		" Reads control points from model "
		ctrlpts = modelIGA._COORDS[:self.dim, :]
		return ctrlpts
	
	def __setQuadratureRules(self, quadArgs:dict):

		name = quadArgs.get('quadrule', 'wq').lower()
		quadRuleByDirection = []
		
		print('Evaluating basis and weights')
		start = time.process_time()
		if name == 'iga':
			for i in range(self.dim):
				quadRuleByDirection.append(GaussQuadrature(self.degree[i], self.knotvector[i], quadArgs=quadArgs))
		if name == 'wq':
			for i in range(self.dim):
				quadRuleByDirection.append(WeightedQuadrature(self.degree[i], self.knotvector[i], quadArgs=quadArgs))
		
		for i, quadRule in enumerate(quadRuleByDirection):
			quadPtsPos, dersIndices, dersBasis, dersWeights = quadRule.getQuadratureRulesInfo()
			nbqp = len(quadPtsPos); indi, indj = dersIndices
			self.nbqp[i] = nbqp; self.indices.append(indi); self.indices.append(indj)
			self.basis.append(dersBasis); self.weights.append(dersWeights)
		stop = time.process_time()
		print('\tBasis and weights in : %.5f s' %(stop-start))

		# Update number of quadrature points
		self.nbqp_total = np.prod(self.nbqp)

		return
	
	def addQuadratureRule(self, direction:int, quadArgs:dict):
		if direction<0 or direction>=self.dim: raise Warning('Direction not valid')
		name = quadArgs.get('quadrule', None)
		if name == 'iga':
				quadRule = GaussQuadrature(self.degree[direction], self.knotvector[direction], quadArgs=quadArgs)
		elif name == 'wq':
				quadRule = WeightedQuadrature(self.degree[direction], self.knotvector[direction], quadArgs=quadArgs)
		else: raise Warning('Insert a valid quadrature rule')

		# Update quadrature rule
		quadPtsPos, dersIndices, dersBasis, dersWeights = quadRule.getQuadratureRulesInfo()
		nbqp = len(quadPtsPos); indi, indj = dersIndices
		self.nbqp[direction] = nbqp; self.indices[2*direction] = indi; self.indices[2*direction+1] = indj
		self.basis[direction] = dersBasis; self.weights[direction] = dersWeights
		self.nbqp_total = np.prod(self.nbqp)
		return
	
	def __setJacobienPhysicalPoints(self):
		" Computes jacobien and physical position "

		print('Evaluating jacobien and physical position')
		start = time.process_time()
		inputs = [*self.nbqp[:self.dim], *self.indices, *self.basis, self.ctrlpts]
		if self.dim == 2:
			self.Jqp = geophy.eval_jacobien_2d(*inputs)
			self.detJ, self.invJ = geophy.eval_inverse_det(self.Jqp)
			self.qpPhy = geophy.interpolate_meshgrid_2d(*inputs)
		if self.dim == 3:
			self.Jqp = geophy.eval_jacobien_3d(*inputs)
			self.detJ, self.invJ = geophy.eval_inverse_det(self.Jqp)
			self.qpPhy = geophy.interpolate_meshgrid_3d(*inputs)
		stop = time.process_time()
		print('\t Time jacobien: %.5f s' %(stop-start))
		return
	
	# ----------------
	# POST-PROCESSING 
	# ----------------

	def interpolateMeshgridField(self, u_ctrlpts=None, sampleSize=101, isAll=True):
		# Initialize all outputs
		meshinterp, Jinterp, detJinterp, uinterp = None, None, None, None

		# Get basis using fortran
		knots = np.linspace(0, 1, sampleSize)
		basis, indices = [], []
		for i in range(self.dim):
			dersb, indi, indj = evalDersBasisFortran(self.degree[i], self.knotvector[i], knots)
			basis.append(dersb); indices.append(indi); indices.append(indj)

		# Get position and determinant 
		if isAll:
			inputs = [*self.dim*[sampleSize], *indices, *basis, self.ctrlpts]
			if self.dim == 2:
				Jinterp = geophy.eval_jacobien_2d(*inputs)
				detJinterp = geophy.eval_inverse_det(Jinterp)[0]
				meshinterp = geophy.interpolate_meshgrid_2d(*inputs)
			elif self.dim == 3: 
				Jinterp = geophy.eval_jacobien_3d(*inputs)
				detJinterp = geophy.eval_inverse_det(Jinterp)[0]
				meshinterp = geophy.interpolate_meshgrid_3d(*inputs)

		if u_ctrlpts is not None: 
			u_temp = np.atleast_2d(u_ctrlpts)
			nr = np.size(u_temp, axis=0)
			inputs = [*self.dim*[sampleSize], *indices, *basis, u_temp]

			if self.dim == 2:   uinterp = geophy.interpolate_meshgrid_2d(*inputs)    
			elif self.dim == 3: uinterp = geophy.interpolate_meshgrid_3d(*inputs)
			if nr == 1: uinterp = np.ravel(uinterp)

		return meshinterp, Jinterp, detJinterp, uinterp

	def exportResultsCP(self, fields={}, folder=None, sampleSize=101): 
		""" Export solution in VTK format. 
			It is possible to use Paraview to visualize data
		"""

		if folder == None: 
			full_path = os.path.realpath(__file__)
			dirname = os.path.dirname
			folder = dirname(dirname(full_path)) + '/results/'
		if not os.path.isdir(folder): os.mkdir(folder)
		print("File saved in %s" %folder)

		# Only interpolate meshgrid
		meshinterp, _, detJinterp = self.interpolateMeshgridField(sampleSize=sampleSize)[:-1]
		detJinterp /= statistics.mean(detJinterp)

		shape = [1, 1, 1]
		for dim in range(self.dim): shape[dim] = sampleSize
		shape = tuple(shape)
		X = [np.zeros(shape) for i in range(3)]
		for i in range(2): X[i] = np.reshape(np.ravel(meshinterp[i, :]), shape, order='F')
		if self.dim == 3:  X[2] = np.reshape(np.ravel(meshinterp[-1, :]), shape, order='F')
		
		# Create point data 
		pointData = {}
		for fieldname, fieldvalue in fields.items():
			fieldvalue = np.atleast_2d(fieldvalue); nr = np.size(fieldvalue, axis=0)
			fieldinterp = self.interpolateMeshgridField(u_ctrlpts=fieldvalue, sampleSize=sampleSize, isAll=False)[-1]
			if nr>1:
				for l in range(nr):
					newfieldname = fieldname + str(l+1)
					pointData[newfieldname] = np.reshape(np.ravel(fieldinterp[l, :]), shape, order='F')
			else:
				pointData[fieldname] = np.reshape(np.ravel(fieldinterp), shape, order='F')
		pointData['detJ'] = np.reshape(np.ravel(detJinterp), shape, order='F')

		# Export geometry
		try: name    = self.name
		except: name = 'ExportFile'
		name = folder + name
		gridToVTK(name, X[0], X[1], X[2], pointData=pointData)

		return

