from .__init__ import *
from .lib_base import (legendreTable, lobattoTable, findMultiplicity, 
					evalDersBasisCSRFortran, array2csr_matrix, evalDersBasisDensePy, 
					evalDersBasisCSRPy, increaseMultiplicity
)

class QuadratureRules:
	def __init__(self, degree, knotvector):
		self.degree       = degree
		self.knotvector   = knotvector
		self.__getInfo()

		self.nbqp         = None
		self.quadPtsPos   = None
		self.dersIndices  = None
		self.dersBasis    = None
		self.dersWeights  = None
		self._denseBasis   = None
		self._denseWeights = None
		return

	def __getInfo(self):
		self.nbctrlpts = len(self.knotvector) - self.degree - 1
		self._uniqueKV = np.unique(self.knotvector)
		self._nbel     = len(self._uniqueKV) - 1
		return
	
	def getQuadratureRulesInfo(self):
		return self.quadPtsPos, self.dersIndices, self.dersBasis, self.dersWeights
	
	def getDenseQuadRules(self, isfortran=True):
		denseBasis, denseWeights = [], []
		indi, indj = self.dersIndices
		for i in range(np.size(self.dersBasis, axis=1)):
			tmp = array2csr_matrix(self.dersBasis[:, i], indi, indj, isfortran=isfortran)
			denseBasis.append(tmp)

		for i in range(np.size(self.dersWeights, axis=1)):
			tmp = array2csr_matrix(self.dersWeights[:, i], indi, indj, isfortran=isfortran)
			denseWeights.append(tmp)
		self._denseBasis, self._denseWeights = denseBasis, denseWeights
		return denseBasis, denseWeights
	
	def getSampleBasis(self, sampleSize=101):
		basis = []
		knots = np.linspace(0, 1, sampleSize)
		dersBasis, indi, indj = evalDersBasisCSRFortran(self.degree, self.knotvector, knots)
		for i in range(np.size(dersBasis, axis=1)):
			tmp = sp.csr_matrix((dersBasis[:, i], indj-1, indi-1)) # -1 because indj is in fortran notation
			basis.append(tmp)
		return basis, knots
	
class GaussQuadrature(QuadratureRules):
	def __init__(self, degree, knotvector, quadArgs:dict):
		QuadratureRules.__init__(self, degree, knotvector)
		self._gaussType = quadArgs.get('type', 'leg').lower()
		extraArgs = quadArgs.get('extra', {'nb_qp_el': 1})
		if self._gaussType == 'leg': # Legendre 
			self._order = self.degree + 1
			self._tablefunction = legendreTable
		elif self._gaussType == 'lob': # Lobatto
			self._order = self.degree + 2
			self._tablefunction = lobattoTable
		elif self._gaussType == 'legextra':
			self._order = self.degree + extraArgs.get('nb_qp_el')
			self._tablefunction = legendreTable
		else: 
			raise Warning('Not found')
		return
	
	def __getGaussInfo(self):
		" Gets the position of Gauss quadrature points in isoparametric space and its weights using known tables "
		self._isoPositions, self._isoWeights = self._tablefunction(self._order)
		return

	def __findQuadraturePositions(self):
		" Gets the position of Gauss quadrature points in parametric space "
		tmp   = np.array([])
		knots = self._uniqueKV
		for i in range(self._nbel):
			xg  = 0.5*((knots[i+1] - knots[i])*self._isoPositions + knots[i] + knots[i+1])
			tmp = np.append(tmp, xg)
		self.quadPtsPos = np.atleast_1d(tmp)
		self.nbqp = np.size(self.quadPtsPos)
		return
	
	def __findParametricWeights(self):
		" Gets the weight of Gauss quadrature points in parametric space "
		tmp   = np.array([])
		knots = self._uniqueKV
		for i in range(self._nbel):
			wg  = 0.5*(knots[i+1] - knots[i])*self._isoWeights
			tmp = np.append(tmp, wg)
		self._parametricWeights = np.atleast_1d(tmp)
		return
	
	def evalDersBasisWeights(self):
		" Gets the basis and weights evaluated at the Gauss quadrature points "
		basis, indi, indj = basisweights.get_genbasis_csr(self.degree, self.knotvector, self.quadPtsPos)
		self.dersBasis   = basis
		self.dersIndices = [indi, indj]
		nnz = np.shape(basis)[0]
		weights = np.zeros((nnz, 4))
		for i in range(nnz):
			weights[i, 0]  = basis[i, 0]*self._parametricWeights[indj[i]-1] # -1 because indj is in fortran notation
			weights[i, 3]  = basis[i, 1]*self._parametricWeights[indj[i]-1]
		weights[:, 1] = weights[:, 0]; weights[:, 2] = weights[:, 3]
		self.dersWeights = weights
		return
	
	def getQuadratureRulesInfo(self):
		self.__getGaussInfo()
		self.__findQuadraturePositions()
		self.__findParametricWeights()
		self.evalDersBasisWeights()
		self.getDenseQuadRules()
		return self.quadPtsPos, self.dersIndices, self.dersBasis, self.dersWeights

class WeightedQuadrature(QuadratureRules):
	def __init__(self, degree, knotvector, quadArgs:dict):
		QuadratureRules.__init__(self, degree, knotvector)
		self._wqType  = quadArgs.get('type', 1)
		if   self._wqType == 1: 
			self._posRule = 'midpoint'
			extraArgsDefault = {'s': 1, 'r': 3}
		elif self._wqType == 2: 
			self._posRule = 'midpoint'
			extraArgsDefault = {'s': 2, 'r': 3}
		elif self._wqType == 3:
			self._posRule = 'internal'
			extraArgsDefault = {'s': 2, 'r': 4}
		self._extraArgs = quadArgs.get('extra', extraArgsDefault)
		return
	
	def __findQuadraturePositions(self):
		if self._posRule == 'midpoint':
			s = self._extraArgs.get('s', 1); r = self._extraArgs.get('r', 2)
			self.quadPtsPos = self.__QuadPosMidPointRule(s=s, r=r)
		elif self._posRule == 'internal':
			s = self._extraArgs.get('s', 2); r = self._extraArgs.get('r', 4)
			self.quadPtsPos = self.__QuadPosInternalRule(s=s, r=r)
		# Add more methods 
		self._Bshape = self.__getBShape(self.quadPtsPos)
		self.nbqp    = np.size(self.quadPtsPos)
		return

	# Add rules to get the quadrature points
	def __QuadPosMidPointRule(self, s=2, r=3):	
		""" Find quadrature points using mid-point rule
			r is the 'extra' points on the knotspans at the boundaries
			s is the 'extra' points on the inner knotspans. The number must respect the max rule (see bibliography)
		"""
		quadPtsPos = np.array([])
		knots = self._uniqueKV
		
		# First span
		tmp = np.linspace(knots[0], knots[1], self.degree+r)
		quadPtsPos = np.append(quadPtsPos, tmp)

		# Last span
		tmp = np.linspace(knots[-2], knots[-1], self.degree+r)
		quadPtsPos = np.append(quadPtsPos, tmp)

		# Inner span
		for i in range(1, self._nbel-1):
			tmp = np.linspace(knots[i], knots[i+1], 2+s)
			quadPtsPos = np.append(quadPtsPos, tmp)
		
		quadPtsPos = np.unique(quadPtsPos)
		quadPtsPos.sort()
		return quadPtsPos
	
	def __QuadPosInternalRule(self, s=2, r=4):	
		""" Find quadrature points using mid-point rule
			r is the 'extra' points on the knotspans at the boundaries
			s is the 'extra' points on the inner knotspans. The number must respect the max rule (see bibliography)
		"""
		quadPtsPos = np.array([])
		knots = self._uniqueKV
		
		# First span
		tmp = np.linspace(knots[0], knots[1], self.degree+r)[1:-1]
		quadPtsPos = np.append(quadPtsPos, tmp)

		# Last span
		tmp = np.linspace(knots[-2], knots[-1], self.degree+r)[1:-1]
		quadPtsPos = np.append(quadPtsPos, tmp)

		# Inner span
		for i in range(1, self._nbel-1):
			tmp = np.linspace(knots[i], knots[i+1], 2+s)[1:-1]
			quadPtsPos = np.append(quadPtsPos, tmp)
		
		quadPtsPos = np.unique(quadPtsPos)
		quadPtsPos.sort()
		return quadPtsPos
	
	def __compute_weights(self, Ishape, Bshape, II, BB):
		
		def solveLinearSystem(B, I):
			sol = np.linalg.lstsq(B, I, rcond=None)[0]
			w = sol.reshape((1, -1)).tolist()[0]
			return w
		
		nr_obj, _ = np.shape(Bshape)
		nr_test, _ = np.shape(Ishape)
		_, nc = np.shape(BB)
		weights = np.zeros((nr_obj, nc))
		for i in range(nr_obj):
			Pmin, Pmax = Bshape[i, :]
			Fmin = 0
			for j in range(nr_test):
				if Ishape[j, i]>0:
					Fmin = j
					break
			Fmax = nr_test - 1
			for j in range(nr_test-1, -1, -1):
				if Ishape[j, i]>0:
					Fmax = j
					break
			tmp = solveLinearSystem(BB[Fmin:Fmax+1, Pmin:Pmax+1], II[Fmin:Fmax+1, i])
			weights[i, Pmin:Pmax+1] = tmp
		return weights
	
	def __getBShape(self, knots, isfortran=True):
		" Return the shape of basis in WQ approach. "
		B0shape = np.zeros((self.nbctrlpts, 2), dtype=int)
		B1shape = np.zeros((self.nbctrlpts, 2), dtype=int)
		knots   = np.atleast_1d(knots)

		tableKVSpan = np.zeros((self._nbel, 2))
		for i in range(0, self._nbel):
			tableKVSpan[i, :] = [self._uniqueKV[i], self._uniqueKV[i+1]]
		tablePointsOverKVSpan = np.zeros((self._nbel, 2), dtype=int)
		for i in range(0, self._nbel):
			left = self._uniqueKV[i]; right = self._uniqueKV[i+1]
			boolean = (knots>=left)*(knots<=right)
			tablePointsOverKVSpan[i, 0] = np.nonzero(boolean)[0][0]
			tablePointsOverKVSpan[i, 1] = np.nonzero(boolean)[0][-1]

		tableFunctionsOverKVSpan = np.zeros((self._nbel, self.degree+1), dtype=int)
		for j in range(0, self.degree+1): tableFunctionsOverKVSpan[0, j] = j
		for i in range(1, self._nbel): 
			multiplicity = findMultiplicity(self.knotvector, self._uniqueKV[i])
			tableFunctionsOverKVSpan[i, 0] = tableFunctionsOverKVSpan[i-1, 0] + multiplicity
			for j in range(1, self.degree+1): 
				tableFunctionsOverKVSpan[i, j] = tableFunctionsOverKVSpan[i, 0] + j

		tableFunctionSpans = np.zeros((self.nbctrlpts, 2), dtype=int)
		for i in range(0, self.nbctrlpts):
			minFuncSpan = 1
			for j in range(0, self._nbel):
				if np.any(tableFunctionsOverKVSpan[j, :]==i):
					minFuncSpan = j; break
			maxFuncSpan = self._nbel
			for j in range(self._nbel-1, -1, -1):
				if np.any(tableFunctionsOverKVSpan[j, :]==i):
					maxFuncSpan = j; break
			tableFunctionSpans[i, :] = [minFuncSpan, maxFuncSpan]

		minFuncSpan, maxFuncSpan = tableFunctionSpans[0, :]
		minKnotOverFuncSpan = tableKVSpan[minFuncSpan, 0]
		maxKnotOverFuncSpan = tableKVSpan[maxFuncSpan, 1]
		boolean = (knots>=minKnotOverFuncSpan)*(knots<maxKnotOverFuncSpan)
		B0shape[0, 0] = np.nonzero(boolean)[0][0]
		B0shape[0, 1] = np.nonzero(boolean)[0][-1]

		for i in range(1, self.nbctrlpts-1):
			minFuncSpan, maxFuncSpan = tableFunctionSpans[i, :]
			minKnotOverFuncSpan = tableKVSpan[minFuncSpan, 0]
			maxKnotOverFuncSpan = tableKVSpan[maxFuncSpan, 1]
			boolean = (knots>minKnotOverFuncSpan)*(knots<maxKnotOverFuncSpan)
			B0shape[i, 0] = np.nonzero(boolean)[0][0]
			B0shape[i, 1] = np.nonzero(boolean)[0][-1]

		minFuncSpan, maxFuncSpan = tableFunctionSpans[self.nbctrlpts-1, :]
		minKnotOverFuncSpan = tableKVSpan[minFuncSpan, 0]
		maxKnotOverFuncSpan = tableKVSpan[maxFuncSpan, 1]
		boolean = (knots>minKnotOverFuncSpan)*(knots<=maxKnotOverFuncSpan)
		B0shape[self.nbctrlpts-1, 0] = np.nonzero(boolean)[0][0]
		B0shape[self.nbctrlpts-1, 1] = np.nonzero(boolean)[0][-1]
	
		B1shape = np.copy(B0shape) # Eventually it could be different
		if isfortran: B0shape += 1
		if isfortran: B1shape += 1 # change from 0 to 1 index
		return [B0shape, B1shape]
	
	def __getWeightsM1(self, isfortran=True):
		" Computes the weights at quadrature points in WQ approach using trial space S^[p-1]_[r-1]. "

		# Space S^p_r
		gauss_p0 = GaussQuadrature(self.degree, self.knotvector, {'type':'leg'})
		gauss_p0.getQuadratureRulesInfo()
		B0cgg_p0, B1cgg_p0 = gauss_p0._denseBasis
		basis_csr, indi_csr, indj_csr = evalDersBasisCSRPy(self.degree, self.knotvector, self.quadPtsPos, isfortran=isfortran)
		B0wq_p0 = array2csr_matrix(basis_csr[:, 0], indi_csr, indj_csr, isfortran=isfortran)

		# Space S^[p-1]_[r-1]
		degree_p1 = self.degree - 1
		knotvector_p1 = self.knotvector[1:-1]
		gauss_p1 = QuadratureRules(degree_p1, knotvector_p1)
		B0cgg_p1 = evalDersBasisDensePy(degree_p1, knotvector_p1, gauss_p0.quadPtsPos, isfortran=isfortran)[0]
		B0wq_p1 = evalDersBasisDensePy(degree_p1, knotvector_p1, self.quadPtsPos, isfortran=isfortran)[0]  

		# Compute Integrals
		weights_csr = np.zeros((np.size(basis_csr, axis=0), 4))
		Bcgg_p0_int = np.zeros((gauss_p0.nbctrlpts, gauss_p0.nbqp))
		Bcgg_p1_int = np.zeros((gauss_p1.nbctrlpts, gauss_p0.nbqp))
		Bcgg_p0_int[np.where(np.abs(B0cgg_p0.todense())>1e-8)] = 1
		Bcgg_p1_int[np.where(np.abs(B0cgg_p1.todense())>1e-8)] = 1
		
		tmp = []
		# Computation of W00
		Ishape = Bcgg_p0_int @ Bcgg_p0_int.T 
		II = B0cgg_p0 @ sp.diags(gauss_p0._parametricWeights) @ B0cgg_p0.T
		tmp.append(self.__compute_weights(Ishape, self._Bshape[0]-1, II.todense(), B0wq_p0.todense()))

		# Computation of W01
		Ishape = Bcgg_p1_int @ Bcgg_p0_int.T 
		II = B0cgg_p1 @ sp.diags(gauss_p0._parametricWeights) @ B0cgg_p0.T 
		tmp.append(self.__compute_weights(Ishape, self._Bshape[0]-1, II.todense(), B0wq_p1.todense()))
		
		# Computation of W10
		Ishape = Bcgg_p0_int @ Bcgg_p0_int.T 
		II = B0cgg_p0 @ sp.diags(gauss_p0._parametricWeights) @ B1cgg_p0.T
		tmp.append(self.__compute_weights(Ishape, self._Bshape[1]-1, II.todense(), B0wq_p0.todense()))
			
		# Computation of W11
		Ishape = Bcgg_p1_int @ Bcgg_p0_int.T
		II = B0cgg_p1 @ sp.diags(gauss_p0._parametricWeights) @ B1cgg_p0.T    
		tmp.append(self.__compute_weights(Ishape, self._Bshape[1]-1, II.todense(), B0wq_p1.todense()))

		c = 0
		indi_csr_copy = np.copy(indi_csr) 
		indj_csr_copy = np.copy(indj_csr) 
		if isfortran: indi_csr_copy -= 1
		if isfortran: indj_csr_copy -= 1
		for i in range(self.nbctrlpts):
			for j in np.array(indj_csr_copy[indi_csr_copy[i]:indi_csr_copy[i+1]]):
				weights_csr[c, :] = [tmp[0][i, j], tmp[1][i, j], tmp[2][i, j], tmp[3][i, j]]
				c += 1

		return basis_csr, weights_csr, indi_csr, indj_csr
	
	# def __getWeightsM2(self, isfortran=True):
	# 	" Computes the weights at quadrature points in WQ approach using trial space S^[p]_[r-1]. "

	# 	# Space S^p_r
	# 	gauss_p0 = GaussQuadrature(self.degree, self.knotvector, {'type':'leg'})
	# 	gauss_p0.getQuadratureRulesInfo()
	# 	B0cgg_p0, B1cgg_p0 = gauss_p0._denseBasis
	# 	basis_csr, indi_csr, indj_csr = evalDersBasisCSRPy(self.degree, self.knotvector, self.quadPtsPos, isfortran=isfortran)
	# 	B0wq_p0 = array2csr_matrix(basis_csr[:, 0], indi_csr, indj_csr, isfortran=isfortran)

	# 	# Space S^[p]_[r-1]
	# 	degree_p1 = self.degree
	# 	knotvector_p1 = increaseMultiplicity(1, degree_p1, self.knotvector)
	# 	gauss_p1 = QuadratureRules(degree_p1, knotvector_p1)
	# 	B0cgg_p1 = evalDersBasisDensePy(degree_p1, knotvector_p1, gauss_p0.quadPtsPos, isfortran=isfortran)[0]
	# 	B0wq_p1 = evalDersBasisDensePy(degree_p1, knotvector_p1, self.quadPtsPos, isfortran=isfortran)[0]  

	# 	# Compute Integrals
	# 	weights_csr = np.zeros((np.size(basis_csr, axis=0), 4))
	# 	Bcgg_p0_int = np.zeros((gauss_p0.nbctrlpts, gauss_p0.nbqp))
	# 	Bcgg_p1_int = np.zeros((gauss_p1.nbctrlpts, gauss_p0.nbqp))
	# 	Bcgg_p0_int[np.where(np.abs(B0cgg_p0.todense())>1e-8)] = 1
	# 	Bcgg_p1_int[np.where(np.abs(B0cgg_p1.todense())>1e-8)] = 1
		
	# 	tmp = []
	# 	# Computation of W00
	# 	Ishape = Bcgg_p1_int @ Bcgg_p0_int.T 
	# 	II = B0cgg_p0 @ sp.diags(gauss_p0._parametricWeights) @ B0cgg_p0.T
	# 	tmp.append(self.__compute_weights(Ishape, self._Bshape[0]-1, II.todense(), B0wq_p0.todense()))

	# 	# Computation of W11
	# 	Ishape = Bcgg_p1_int @ Bcgg_p0_int.T
	# 	II = B0cgg_p1 @ sp.diags(gauss_p0._parametricWeights) @ B1cgg_p0.T    
	# 	tmp.append(self.__compute_weights(Ishape, self._Bshape[1]-1, II.todense(), B0wq_p1.todense()))

	# 	c = 0
	# 	indi_csr_copy = np.copy(indi_csr) 
	# 	indj_csr_copy = np.copy(indj_csr) 
	# 	if isfortran: indi_csr_copy -= 1
	# 	if isfortran: indj_csr_copy -= 1
	# 	for i in range(self.nbctrlpts):
	# 		for j in np.array(indj_csr_copy[indi_csr_copy[i]:indi_csr_copy[i+1]]):
	# 			weights_csr[c, :] = [tmp[0][i, j], tmp[0][i, j], tmp[1][i, j], tmp[1][i, j]]
	# 			c += 1

	# 	return basis_csr, weights_csr, indi_csr, indj_csr
	
	def evalDersBasisWeights(self):
		assert self._wqType in [1, 2, 3], 'Method unknown'
		size_data = (self.degree + 1)*self.nbqp
		basis   = np.zeros((size_data, 2))
		weights = np.zeros((size_data, 4))
		indj = np.zeros(size_data, dtype=int)
		indi = np.zeros(self.nbctrlpts+1, dtype=int)
		if (self._wqType == 1 or self._wqType == 3) and self.degree == 1: 
			basis, weights, indi, indj = self.__getWeightsM1()
		else: 
			basis, weights, indi, indj = basisweights.wq_getbasisweights_csr(self.degree, self.knotvector, self.quadPtsPos, 
													self._Bshape[0], self._Bshape[1], size_data, self._wqType)
		self.dersBasis   = basis
		self.dersWeights = weights
		self.dersIndices = [indi, indj]

		return
	
	def getQuadratureRulesInfo(self):
		self.__findQuadraturePositions()
		self.evalDersBasisWeights()
		self.getDenseQuadRules()
		return self.quadPtsPos, self.dersIndices, self.dersBasis, self.dersWeights
	