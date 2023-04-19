"""
.. This module helps to read information from a given geometry (YETI format)
.. and to create thermo-mechanical model (material, boundary conditions, etc.)
.. Joaquin Cornejo
"""

# My libraries
from lib.lib_base import evalDersBasisFortran
from lib.lib_quadrules import *
from lib.lib_material import *

class part(): 

	def __init__(self, modelIGA, **kwargs):
		print('\nInitializing thermo-mechanical model')
		self._sampleSize = kwargs.get('sample_size', 101)
		self._name       = self.read_name(modelIGA)
		self._dim        = self.read_dimension(modelIGA)
		self._degree     = self.read_degree(modelIGA)
		self._ctrlpts    = self.read_controlPoints(modelIGA)
		self._nbctrlpts  = np.ones(3, dtype=int)
		self._kwargs     = kwargs
		self._knotvector, self._size_kv = self.read_knotvector(modelIGA)
		for i in range(self._dim): self._nbctrlpts[i] = len(self._knotvector[i]) - self._degree[i] - 1
		self._nbctrlpts_total = np.product(self._nbctrlpts)

		self._nbqp, self._nbqp_total = [], None
		self._qpPar, self._basis, self._weights, self._indices = [], [], [], []
		self.setQuadratureRules()

		self._Jqp, self._detJ, self._invJ, self._qpPhy = None, None, None, None
		self.setJacobienPhysicalPoints()

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

	def setQuadratureRules(self):

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
				quadRule_list.append(GaussQuadrature(self._degree[i], self._knotvector[i], kwargs=kwargs))
		elif quadRuleName == 'wq':
			for i in range(self._dim):
				quadRule_list.append(WeightedQuadrature(self._degree[i], self._knotvector[i], kwargs=kwargs))
		else:
			raise Warning('Not found')
		
		for quadRule in quadRule_list:
			nbqp, quadPtsPos, dersIndices, dersBasis, dersWeights = eval_basisWeights(quadRule)
			indi, indj = dersIndices
			self._nbqp.append(nbqp); self._qpPar.append(quadPtsPos); self._indices.append(indi); 
			self._indices.append(indj); self._basis.append(dersBasis); self._weights.append(dersWeights)
		stop = time.process_time()
		print('\tBasis and weights in : %.5f s' %(stop-start))

		# Update number of quadrature points
		self._nbqp_total = np.prod(self._nbqp)

		return
	
	def setJacobienPhysicalPoints(self):
		" Computes jacobien and physical position "

		print('Evaluating jacobien and physical position')
		start = time.process_time()
		inputs = [*self._nbqp, *self._indices, *self._basis, self._ctrlpts]
		if self._dim == 2:
			self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_2d(*inputs)
			self._qpPhy = assembly.interpolate_fieldphy_2d(*inputs)
		if self._dim == 3:
			self._Jqp, self._detJ, self._invJ = assembly.eval_jacobien_3d(*inputs)
			self._qpPhy = assembly.interpolate_fieldphy_3d(*inputs)
		stop = time.process_time()
		print('\t Time jacobien: %.5f s' %(stop-start))
		return
	
	# ----------------
	# POST-PROCESSING 
	# ----------------

	def interpolateField(self, samplesize=None, u_ctrlpts=None, nbDOF=3):

		# Get basis using fortran
		if samplesize == None: samplesize = self._sampleSize
		knots = np.linspace(0, 1, samplesize)
		basis, indices = [], []
		for i in range(self._dim):
			dersb, indi, indj = evalDersBasisFortran(self._degree[i], self._knotvector[i], knots)
			basis.append(dersb); indices.append(indi); indices.append(indj)


		# Get position and determinant 
		inputs = [*self._dim*[samplesize], *indices, *basis, self._ctrlpts]
		if self._dim == 2:
			Jinterp, detJinterp = assembly.eval_jacobien_2d(*inputs)[:2]
			posinterp = assembly.interpolate_fieldphy_2d(*inputs)
		elif self._dim == 3: 
			Jinterp, detJinterp = assembly.eval_jacobien_3d(*inputs)[:2]
			posinterp = assembly.interpolate_fieldphy_3d(*inputs)

		# Get interpolation
		if u_ctrlpts is not None:
			u_temp = np.atleast_2d(u_ctrlpts)
			inputs = [*self._dim*[samplesize], *indices, *basis, u_temp]

			if self._dim == 2:   uinterp = assembly.interpolate_fieldphy_2d(*inputs)    
			elif self._dim == 3: uinterp = assembly.interpolate_fieldphy_3d(*inputs)
			if nbDOF == 1: uinterp = np.ravel(uinterp)
	
		else: uinterp = None

		return Jinterp, posinterp, detJinterp, uinterp

	def exportResults(self, u_ctrlpts=None, folder=None, name=None, nbDOF=3): 
		""" Export solution in VTK format. 
			It is possible to use Paraview to visualize data
		"""

		if folder == None: 
			full_path = os.path.realpath(__file__)
			dirname = os.path.dirname
			folder = dirname(dirname(full_path)) + '/results/'
		if not os.path.isdir(folder): os.mkdir(folder)
		print("File saved in %s" %folder)

		if u_ctrlpts is None: pass
		elif isinstance(u_ctrlpts, np.ndarray): 
			if np.size(u_ctrlpts)%self._nb_ctrlpts_total != 0: 
				raise Warning('Not enough control points')
		else: raise Warning('Solution must be ndarray type')

		# ------------------
		# Get interpolation
		# ------------------
		qpPhy, detJ, u_interp = self.interpolateField(u_ctrlpts=u_ctrlpts, nbDOF=nbDOF)[1:]
		mean_detJ = statistics.mean(detJ)
		detJ /= mean_detJ

		# ------------------
		# Export results
		# ------------------
		shape_pts = [1, 1, 1]
		for dim in range(self._dim): shape_pts[dim] = self._sample_size
		shape_pts  = tuple(shape_pts)
		X1, X2, X3 = np.zeros(shape_pts), np.zeros(shape_pts), np.zeros(shape_pts)
		U   = np.zeros((nbDOF, *shape_pts))
		DET = np.zeros(shape_pts)

		for k in range(shape_pts[2]):
			for j in range(shape_pts[1]):
				for i in range(shape_pts[0]):
					pos = i + j * self._sample_size + k * self._sample_size**2
					X1[i,j,k] = qpPhy[0, pos]
					X2[i,j,k] = qpPhy[1, pos]
					DET[i,j,k] = detJ[pos]
					if self._dim == 3: X3[i,j,k] = qpPhy[2, pos]
					if u_interp is not None: 
						u_interp = np.atleast_2d(u_interp)
						for l in range(nbDOF):
							U[l,i,j,k] = u_interp[l, pos]
		
		# Create point data 
		pointData = {}
		if u_interp is not None: 
			for l in range(nbDOF):
				varname = 'U' + str(l+1)
				pointData[varname] = U[l,:,:,:]
		pointData['detJ'] = DET

		# Export geometry
		if name is None: name = self._name
		name = folder + name
		gridToVTK(name, X1, X2, X3, pointData=pointData)
		
		return
