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
			self.Jqp, self.detJ, self.invJ = geophy.eval_jacobien_2d(*inputs)
			self.qpPhy = geophy.interpolate_fieldphy_2d(*inputs)
		if self.dim == 3:
			self.Jqp, self.detJ, self.invJ = geophy.eval_jacobien_3d(*inputs)
			self.qpPhy = geophy.interpolate_fieldphy_3d(*inputs)
		stop = time.process_time()
		print('\t Time jacobien: %.5f s' %(stop-start))
		return
	
	# ----------------
	# POST-PROCESSING 
	# ----------------

	def interpolateField(self, u_ctrlpts=None, nbDOF=3, sampleSize=101):

		# Get basis using fortran
		knots = np.linspace(0, 1, sampleSize)
		basis, indices = [], []
		for i in range(self.dim):
			dersb, indi, indj = evalDersBasisFortran(self.degree[i], self.knotvector[i], knots)
			basis.append(dersb); indices.append(indi); indices.append(indj)


		# Get position and determinant 
		inputs = [*self.dim*[sampleSize], *indices, *basis, self.ctrlpts]
		if self.dim == 2:
			Jinterp, detJinterp = geophy.eval_jacobien_2d(*inputs)[:2]
			posinterp = geophy.interpolate_fieldphy_2d(*inputs)
		elif self.dim == 3: 
			Jinterp, detJinterp = geophy.eval_jacobien_3d(*inputs)[:2]
			posinterp = geophy.interpolate_fieldphy_3d(*inputs)

		# Get interpolation
		if u_ctrlpts is not None:
			u_temp = np.atleast_2d(u_ctrlpts)
			inputs = [*self.dim*[sampleSize], *indices, *basis, u_temp]

			if self.dim == 2:   uinterp = geophy.interpolate_fieldphy_2d(*inputs)    
			elif self.dim == 3: uinterp = geophy.interpolate_fieldphy_3d(*inputs)
			if nbDOF == 1: uinterp = np.ravel(uinterp)
	
		else: uinterp = None

		return Jinterp, posinterp, detJinterp, uinterp

	def exportResults(self, u_ctrlpts=None, folder=None, name=None, nbDOF=3, sampleSize=101): 
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
			if np.size(u_ctrlpts)%self.nbctrlpts_total != 0: 
				raise Warning('Not enough control points')
		else: raise Warning('Solution must be ndarray type')

		# ------------------
		# Get interpolation
		# ------------------
		qpPhy, detJ, u_interp = self.interpolateField(u_ctrlpts=u_ctrlpts, nbDOF=nbDOF, sampleSize=sampleSize)[1:]
		mean_detJ = statistics.mean(detJ)
		detJ /= mean_detJ

		# ------------------
		# Export results
		# ------------------
		shape_pts = [1, 1, 1]
		for dim in range(self.dim): shape_pts[dim] = sampleSize
		shape_pts  = tuple(shape_pts)
		X1, X2, X3 = np.zeros(shape_pts), np.zeros(shape_pts), np.zeros(shape_pts)
		U   = np.zeros((nbDOF, *shape_pts))
		DET = np.zeros(shape_pts)

		for k in range(shape_pts[2]):
			for j in range(shape_pts[1]):
				for i in range(shape_pts[0]):
					pos = i + j * sampleSize + k * sampleSize**2
					X1[i,j,k] = qpPhy[0, pos]
					X2[i,j,k] = qpPhy[1, pos]
					DET[i,j,k] = detJ[pos]
					if self.dim == 3: X3[i,j,k] = qpPhy[2, pos]
					if u_interp is not None: 
						u_interp = np.atleast_2d(u_interp)
						for l in range(nbDOF):
							U[l,i,j,k] = u_interp[l, pos]
		
		# Create point data 
		pointData = {}
		if u_interp is not None: 
			for l in range(nbDOF):
				varname = 'U' + str(l+1)
				pointData[varname] = U[l, :, :, :]
		pointData['detJ'] = DET

		# Export geometry
		if name is None: name = self.name
		name = folder + name
		gridToVTK(name, X1, X2, X3, pointData=pointData)
		
		return
