"""
.. This module helps to read information from a given geometry (YETI format)
.. and to create thermo-mechanical model (material, boundary conditions, etc.)
.. Joaquin Cornejo
"""

# My libraries
from .__init__ import *
from .lib_base import evalDersBasisCSRFortran, array2csr_matrix, findInterpolationSpan
from .lib_quadrules import WeightedQuadrature, GaussQuadrature

class part1D():

	def __init__(self, crv:BSpline.Curve, kwargs:dict):
		self.dim 	    = 1
		self.degree     = crv.degree
		self.knotvector = np.array(crv.knotvector)
		self.ctrlpts    = np.array(crv.ctrlpts)[:, 0]
		self.nbctrlpts  = np.array([len(self.ctrlpts), 1, 1])
		self.nbctrlpts_total = np.product(self.nbctrlpts)

		self.__setQuadratureRules(kwargs.get('quadArgs', {}))
		self.__setJacobienPhysicalPoints()
		assert np.all(self.detJ>0.0), 'Geometry problem. See control points positions'
		return
	
	def __setQuadratureRules(self, quadArgs:dict):
		quadRuleName = quadArgs.get('quadrule', '').lower()
		if quadRuleName == 'iga':
			quadRule = GaussQuadrature(self.degree, self.knotvector, quadArgs=quadArgs)
		elif quadRuleName == 'wq':
			quadRule = WeightedQuadrature(self.degree, self.knotvector, quadArgs=quadArgs)
		else: raise Warning('Not found')
		quadRule.getQuadratureRulesInfo()
		self.quadRule = quadRule
		self._densebasis, self._denseweights = quadRule.getDenseQuadRules()	
		self.nbqp = quadRule.nbqp
		return
	
	def __setJacobienPhysicalPoints(self):
		self.Jqp  = self._densebasis[1].T @ self.ctrlpts
		self.detJ = self.Jqp
		self.invJ = 1.0/self.Jqp
		self.qpPhy = self._densebasis[0].T @ self.ctrlpts
		return
	
	def _compute_discrete_mesh_parameter(self, knots):

		def compute_mesh_parameter_element(part:part1D):
			kvunique = np.unique(part.knotvector)
			dersBasis, indi, indj = evalDersBasisCSRFortran(part.degree, part.knotvector, kvunique)
			B0 = array2csr_matrix(dersBasis[:, 0], indi, indj, isfortran=True)
			nodesPhy = B0.T @ part.ctrlpts
			return np.diff(nodesPhy)

		distPhy = compute_mesh_parameter_element(self); mesh_par = []
		for knot in knots:
			span = findInterpolationSpan(np.unique(self.knotvector), knot)
			mesh_par.append(distPhy[span])
		return np.array(mesh_par)

class part(): 

	def __init__(self, modelIGA, quadArgs:dict):
		print('\nInitializing thermo-mechanical model')
		self.name       = self.__read_name(modelIGA)
		self.dim        = self.__read_dimension(modelIGA)
		self.degree     = self.__read_degree(modelIGA)
		self.NURBSwgt   = self.__read_NURBS_weights(modelIGA)
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
		assert np.all(self.detJ>0.0), 'Geometry problem. See control points positions'

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
		if any(p == 1 for p in degree[:self.dim]): print('ATTENTION: Model must have at least degree 2 when working with weighted quadrature')
		return degree

	def __read_dimension(self, modelIGA:IGAparametrization):
		" Reads dimensions from model "
		dim = modelIGA._dim[0]
		if dim != 3: print("WARNING: Maybe there are functions not well-implemented for 2D geometries")
		return dim

	def __read_knotvector(self, modelIGA:IGAparametrization):
		" Reads knot-vector from model "
		knotvector = modelIGA._Ukv[0]
		size_kv = modelIGA._Nkv.flatten()
		return knotvector, size_kv

	def __read_controlPoints(self, modelIGA:IGAparametrization): 
		" Reads control points from model "
		ctrlpts = modelIGA._COORDS[:self.dim, :]
		for i in range(self.dim):
			ctrlpts[i, :] *= self.NURBSwgt
		return ctrlpts
	
	def __read_NURBS_weights(self, modelIGA:IGAparametrization):
		try: NURBS_weights = modelIGA._vectWeight
		except: NURBS_weights = np.ones(shape=np.shape(modelIGA._COORDS[:self.dim, :]))
		return NURBS_weights
	
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
			_, dersIndices, dersBasis, dersWeights = quadRule.getQuadratureRulesInfo()
			nbqp = quadRule.nbqp; indi, indj = dersIndices
			self.nbqp[i] = nbqp; self.indices.append(indi); self.indices.append(indj)
			self.basis.append(dersBasis); self.weights.append(dersWeights)
		stop = time.process_time()
		print('\tBasis and weights in : %.5f s' %(stop-start))

		# Update number of quadrature points
		self.nbqp_total = np.prod(self.nbqp)

		# Save quadrature rules	
		self._quadraturerules = quadRuleByDirection

		return
	
	def __setJacobienPhysicalPoints(self):
		" Computes jacobien and physical position "

		print('Evaluating jacobien and physical position')
		start = time.process_time()
		inpts = [*self.nbqp[:self.dim], *self.indices, *self.basis, np.atleast_2d(self.NURBSwgt)]

		if self.dim == 2:
			self.NURBScorrection = geophy.interpolate_meshgrid_2d(*inpts)[0, :]
			dersNURBScorrection = geophy.eval_jacobien_2d(*inpts)[0, :, :]

		elif self.dim == 3: 
			self.NURBScorrection = geophy.interpolate_meshgrid_3d(*inpts)[0, :]
			dersNURBScorrection = geophy.eval_jacobien_3d(*inpts)[0, :, :]

		inpts = [*self.nbqp[:self.dim], *self.indices, *self.basis, self.ctrlpts]
		if self.dim == 2:
			self.qpPhy = geophy.interpolate_meshgrid_2d(*inpts)
			pseudojac = geophy.eval_jacobien_2d(*inpts)
		if self.dim == 3:
			self.qpPhy = geophy.interpolate_meshgrid_3d(*inpts)
			pseudojac = geophy.eval_jacobien_3d(*inpts)
		for i in range(self.dim):
			self.qpPhy[i, :] /= self.NURBScorrection
			
		self.Jqp = np.zeros((self.dim, self.dim, len(self.NURBScorrection)))
		for i in range(self.dim):
			for j in range(self.dim):
				self.Jqp[i, j, :] = (pseudojac[i, j, :] - self.qpPhy[i, :]*dersNURBScorrection[j, :])/self.NURBScorrection
		
		self.detJ, self.invJ = geophy.eval_inverse_det(self.Jqp)
		
		stop = time.process_time()
		print('\t Time jacobien: %.5f s' %(stop-start))
		return	

	# ----------------
	# POST-PROCESSING 
	# ----------------
	def __compute_mesh_parameter_element(self, meantype=None):
		if meantype is None: meantype = 'max'
		if meantype.lower() == 'max': func = lambda x: np.max(x)
		elif meantype.lower() == 'arithmetic': func = lambda x: np.mean(x)
		elif meantype.lower() == 'geometric': func = lambda x: np.exp(np.mean(np.log(x)))

		nbqp, indices, basis = [], [], []
		for i in range(self.dim):
			knots = np.unique(self.knotvector[i])
			dersBasis, indi, indj = evalDersBasisCSRFortran(self.degree[i], self.knotvector[i], knots)
			nbqp.append(len(knots)); indices.append(indi); indices.append(indj); basis.append(dersBasis)
		
		inpts = [*nbqp[:self.dim], *indices, *basis, self.ctrlpts]; distmean = []
		if self.dim == 2:
			qpPhy = geophy.interpolate_meshgrid_2d(*inpts)
			for j in range(nbqp[1]-1):
				for i in range(nbqp[0]-1):
					parGridIndices = [i + j*nbqp[0], 
									i + (j + 1)*nbqp[0],
									i + 1 + j*nbqp[0], 
									i + 1 + (j + 1)*nbqp[0]]
					dist1 = np.linalg.norm(qpPhy[:, parGridIndices[0]] - qpPhy[:, parGridIndices[3]])
					dist2 = np.linalg.norm(qpPhy[:, parGridIndices[1]] - qpPhy[:, parGridIndices[2]])
					distmean.append(func([dist1, dist2]))

		if self.dim == 3:
			qpPhy = geophy.interpolate_meshgrid_3d(*inpts)
			for k in range(nbqp[2]-1):
				for j in range(nbqp[1]-1):
					for i in range(nbqp[0]-1):
						parGridIndices = [i + j*nbqp[0] + k*nbqp[0]*nbqp[1], 
										i + j*nbqp[0] + (k + 1)*nbqp[0]*nbqp[1],
										i + (j + 1)*nbqp[0] + k*nbqp[0]*nbqp[1],
										i + (j + 1)*nbqp[0] + (k + 1)*nbqp[0]*nbqp[1], 
										i + 1 + j*nbqp[0] + k*nbqp[0]*nbqp[1], 
										i + 1 + j*nbqp[0] + (k + 1)*nbqp[0]*nbqp[1], 
										i + 1 + (j + 1)*nbqp[0] + k*nbqp[0]*nbqp[1], 
										i + 1 + (j + 1)*nbqp[0] + (k + 1)*nbqp[0]*nbqp[1]]
						dist1 = np.linalg.norm(qpPhy[:, parGridIndices[0]] - qpPhy[:, parGridIndices[7]])
						dist2 = np.linalg.norm(qpPhy[:, parGridIndices[1]] - qpPhy[:, parGridIndices[6]])
						dist3 = np.linalg.norm(qpPhy[:, parGridIndices[2]] - qpPhy[:, parGridIndices[5]])
						dist4 = np.linalg.norm(qpPhy[:, parGridIndices[3]] - qpPhy[:, parGridIndices[4]])
						distmean.append(func([dist1, dist2]))
		return np.array(distmean)

	def compute_global_mesh_parameter(self):
		return np.max(self.__compute_mesh_parameter_element(meantype='max'))

	def interpolateMeshgridField(self, u_ctrlpts=None, sampleSize=None, isAll=True):
		
		if sampleSize is None: sampleSize = np.max(self.nbqp)*np.ones(self.dim, dtype=int)
		elif np.isscalar(sampleSize): sampleSize = np.copy(sampleSize)*np.ones(self.dim, dtype=int)

		# Initialize all outputs
		pts, Jpts, detJ, uinterp = None, None, None, None

		# Get basis using fortran
		basis, indices = [], []
		for i in range(self.dim):
			knots = np.linspace(0, 1, sampleSize[i])
			dersb, indi, indj = evalDersBasisCSRFortran(self.degree[i], self.knotvector[i], knots)
			basis.append(dersb); indices.append(indi); indices.append(indj)

		# Get position and determinant 
		if isAll:
			inpts = [*sampleSize[:self.dim], *indices, *basis, np.atleast_2d(self.NURBSwgt)]
			if self.dim == 2:
				NURBScorrection = geophy.interpolate_meshgrid_2d(*inpts)[0, :]
				dersNURBScorrection = geophy.eval_jacobien_2d(*inpts)[0, :, :]
			elif self.dim == 3: 
				NURBScorrection = geophy.interpolate_meshgrid_3d(*inpts)[0, :]
				dersNURBScorrection = geophy.eval_jacobien_3d(*inpts)[0, :, :]

			inpts = [*sampleSize[:self.dim], *indices, *basis, self.ctrlpts]
			if self.dim == 2:
				pts = geophy.interpolate_meshgrid_2d(*inpts)
				pseudojac = geophy.eval_jacobien_2d(*inpts)
			if self.dim == 3:
				pts = geophy.interpolate_meshgrid_3d(*inpts)
				pseudojac = geophy.eval_jacobien_3d(*inpts)
			for i in range(self.dim):
				pts[i, :] /= NURBScorrection
				
			Jpts = np.zeros((self.dim, self.dim, len(NURBScorrection)))
			for i in range(self.dim):
				for j in range(self.dim):
					Jpts[i, j, :] = (pseudojac[i, j, :] - pts[i, :]*dersNURBScorrection[j, :])/NURBScorrection
			
			detJ = geophy.eval_inverse_det(Jpts)[0]

		if u_ctrlpts is not None: 
			inpts = [*sampleSize[:self.dim], *indices, *basis, np.atleast_2d(self.NURBSwgt)]
			if self.dim == 2:
				NURBScorrection = geophy.interpolate_meshgrid_2d(*inpts)[0, :]
				dersNURBScorrection = geophy.eval_jacobien_2d(*inpts)[0, :, :]
			elif self.dim == 3: 
				NURBScorrection = geophy.interpolate_meshgrid_3d(*inpts)[0, :]
				dersNURBScorrection = geophy.eval_jacobien_3d(*inpts)[0, :, :]

			inpts = [*sampleSize[:self.dim], *indices, *basis, np.atleast_2d(u_ctrlpts)]
			if self.dim == 2:   uinterp = geophy.interpolate_meshgrid_2d(*inpts)    
			elif self.dim == 3: uinterp = geophy.interpolate_meshgrid_3d(*inpts)
			for i in range(np.size(uinterp, axis=0)):
				uinterp[i, :] /= NURBScorrection

		return pts, Jpts, detJ, uinterp

	def postProcessingPrimal(self, fields={}, folder=None, sampleSize=None, addDetJ=False, name='primal', extraArgs={}): 
		""" Export solution in VTK format. 
			It is possible to use Paraview to visualize data
		"""
		from pyevtk.hl import gridToVTK
		if folder is None: 
			full_path = os.path.realpath(__file__)
			dirname = os.path.dirname
			folder = dirname(dirname(full_path)) + '/results/'
		if not os.path.isdir(folder): os.mkdir(folder)
		print("File saved in %s" %folder)

		if sampleSize is None: sampleSize = np.max(self.nbqp)*np.ones(self.dim, dtype=int)
		elif np.isscalar(sampleSize): sampleSize = np.copy(sampleSize)*np.ones(self.dim, dtype=int)

		# Only interpolate meshgrid
		pts, _, detJ = self.interpolateMeshgridField(sampleSize=sampleSize)[:-1]

		ptsshape = [1, 1, 1]
		for dim in range(self.dim): ptsshape[dim] = sampleSize[dim]
		ptsshape = tuple(ptsshape)
		X = [np.zeros(ptsshape) for i in range(3)]
		for i in range(2): X[i] = np.reshape(np.ravel(pts[i, :]), ptsshape, order='F')
		if self.dim == 3:  X[2] = np.reshape(np.ravel(pts[-1, :]), ptsshape, order='F')
		
		# Create point data 
		pointData = {}
		for fieldname, fieldvalue in fields.items():
			if fieldvalue is None: continue
			if isinstance(fieldvalue, np.ndarray):
				fieldvalue = np.atleast_2d(fieldvalue)
				fieldinterp = self.interpolateMeshgridField(u_ctrlpts=fieldvalue, sampleSize=sampleSize, isAll=False)[-1]
			elif callable(fieldvalue):
				if not 'position' in extraArgs.keys(): extraArgs['position'] = pts
				fieldinterp = fieldvalue(extraArgs)
				fieldinterp = np.atleast_2d(fieldinterp)
			else: continue
		
			nr = np.size(fieldinterp, axis=0)
			if nr>1:
				for l in range(nr):
					newfieldname = fieldname + str(l+1)
					pointData[newfieldname] = np.reshape(np.ravel(fieldinterp[l, :]), ptsshape, order='F')
			else:
				pointData[fieldname] = np.reshape(np.ravel(fieldinterp), ptsshape, order='F')

		if addDetJ: pointData['detJ'] = np.reshape(np.ravel(detJ), ptsshape, order='F')

		gridToVTK(folder + name, X[0], X[1], X[2], pointData=pointData)

		return

	def postProcessingDual(self, fields={}, folder=None, addDetJ=False, name='dual'):
		from pyevtk.hl import gridToVTK
		if folder == None: 
			full_path = os.path.realpath(__file__)
			dirname = os.path.dirname
			folder = dirname(dirname(full_path)) + '/results/'
		if not os.path.isdir(folder): os.mkdir(folder)
		print("File saved in %s" %folder)

		ptsshape = tuple(self.nbqp)
		X = [np.zeros(ptsshape) for i in range(3)]
		for i in range(2): X[i] = np.reshape(np.ravel(self.qpPhy[i, :]), ptsshape, order='F')
		if self.dim == 3:  X[2] = np.reshape(np.ravel(self.qpPhy[-1, :]), ptsshape, order='F')

		# Create point data 
		pointData = {}
		for fieldname, fieldvalue in fields.items():
			if fieldvalue is None or not isinstance(fieldvalue, np.ndarray): continue
			fieldinterp = np.atleast_2d(fieldvalue)
			nr = np.size(fieldinterp, axis=0)
			if nr>1:
				for l in range(nr):
					newfieldname = fieldname + str(l+1)
					pointData[newfieldname] = np.reshape(np.ravel(fieldinterp[l, :]), ptsshape, order='F')
			else:
				pointData[fieldname] = np.reshape(np.ravel(fieldinterp), ptsshape, order='F')

		if addDetJ: pointData['detJ'] = np.reshape(np.ravel(self.detJ), ptsshape, order='F')

		# Export geometry
		gridToVTK(folder + name, X[0], X[1], X[2], pointData=pointData)

		return

