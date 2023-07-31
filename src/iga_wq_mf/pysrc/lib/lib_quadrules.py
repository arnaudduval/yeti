from .__init__ import *
from .lib_base import (legendreTable, lobattoTable, 
						findMultiplicity, createUniformMaxregularKnotvector, 
						insertRowCSR, evalDersBasisFortran, array2csr_matrix
)

class QuadratureRules:
	def __init__(self, degree, knotvector):
		self.degree       = degree
		self.knotvector   = knotvector
		self.__getInfo()
		self.__verifyUniformityRegularity()

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
	
	def __verifyUniformityRegularity(self, threshold=1e-8):
		self._isMaxReg = False
		if (self._nbel+self.degree==self.nbctrlpts): self._isMaxReg = True

		self._isUniform = True
		diffknot = np.diff(np.diff(self._uniqueKV))
		if np.any(np.abs(diffknot)>=threshold): self._isUniform = False
		return
	
	def getQuadratureRulesInfo(self):
		return self.quadPtsPos, self.dersIndices, self.dersBasis, self.dersWeights
	
	def getDenseQuadRules(self, isFortran=True):
		denseBasis, denseWeights = [], []
		indi, indj = self.dersIndices
		for i in range(np.size(self.dersBasis, axis=1)):
			tmp = array2csr_matrix(self.dersBasis[:, i], indi, indj, isfortran=isFortran)
			denseBasis.append(tmp)

		for i in range(np.size(self.dersWeights, axis=1)):
			tmp = array2csr_matrix(self.dersWeights[:, i], indi, indj, isfortran=isFortran)
			denseWeights.append(tmp)
		self._denseBasis, self._denseWeights = denseBasis, denseWeights
		return denseBasis, denseWeights
	
	def getSampleBasis(self, sampleSize=101):
		basis = []
		knots = np.linspace(0, 1, sampleSize)
		dersBasis, indi, indj = evalDersBasisFortran(self.degree, self.knotvector, knots)
		for i in range(np.size(dersBasis, axis=1)):
			tmp = sp.csr_matrix((dersBasis[:, i], indj-1, indi-1)) # -1 because indj is in fortran notation
			basis.append(tmp)
		return basis, knots
	
class GaussQuadrature(QuadratureRules):
	def __init__(self, degree, knotvector, quadArgs:dict):
		super().__init__(degree, knotvector)
		self._gaussType = quadArgs.get('type', 'leg').lower()
		extraArgs = quadArgs.get('extra', {'nbQPEL': 1})
		if self._gaussType == 'leg': # Legendre 
			self._order = self.degree + 1
			self._tablefunction = legendreTable
		elif self._gaussType == 'lob': # Lobatto
			self._order = self.degree + 2
			self._tablefunction = lobattoTable
		elif self._gaussType == 'legextra':
			self._order = self.degree + extraArgs.get('nbQPEL')
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
		return self.quadPtsPos, self.dersIndices, self.dersBasis, self.dersWeights

class WeightedQuadrature(QuadratureRules):
	def __init__(self, degree, knotvector, quadArgs:dict):
		super().__init__(degree, knotvector)
		self._wqType  = quadArgs.get('type', 1)
		self._posRule = 'midpoint' # By default
		if   self._wqType == 1: extraArgsDefault = {'s': 1, 'r': 2}
		elif self._wqType == 2: extraArgsDefault = {'s': 2, 'r': 2}
		self._extraArgs = quadArgs.get('extra', extraArgsDefault)
		return
	
	def __findQuadraturePositions(self):
		if self._posRule == 'midpoint':
			s = self._extraArgs.get('s', 1); r = self._extraArgs.get('r', 2)
			self.quadPtsPos = self.__QuadPosMidPointRule(s=s, r=r)
		# Add more methods 
		self._Bshape = self.__getBShape(self.quadPtsPos)
		self.nbqp    = np.size(self.quadPtsPos)
		return

	# Add rules to get the quadrature points
	def __QuadPosMidPointRule(self, s=1, r=2):	
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
	
	def __getBShape(self, knots, isfortran=True):
		" Return the shape of basis in WQ approach. "
		B0shape = np.zeros((self.nbctrlpts, 2), dtype=int)
		B1shape = np.zeros((self.nbctrlpts, 2), dtype=int)
		knots   = np.atleast_1d(knots)

		tablePointsOverSpan = np.zeros((self._nbel, 2), dtype=int)
		for i in range(0, self._nbel):
			left = self._uniqueKV[i]; right = self._uniqueKV[i+1]
			boolean = (knots>=left)*(knots<=right)
			tablePointsOverSpan[i, 0] = np.nonzero(boolean)[0][0]
			tablePointsOverSpan[i, 1] = np.nonzero(boolean)[0][-1]

		tableFunctionsOverSpan = np.zeros((self._nbel, self.degree+1), dtype=int)
		for j in range(0, self.degree+1): tableFunctionsOverSpan[0, j] = j
		for i in range(1, self._nbel): 
			multiplicity = findMultiplicity(self.knotvector, self._uniqueKV[i])
			tableFunctionsOverSpan[i, 0] = tableFunctionsOverSpan[i-1, 0] + multiplicity
			for j in range(1, self.degree+1): 
				tableFunctionsOverSpan[i, j] = tableFunctionsOverSpan[i, 0] + j

		tableSpansOverFunction = np.zeros((self.nbctrlpts, 2), dtype=int)
		for i in range(0, self.nbctrlpts):
			minspan = 1
			for j in range(0, self._nbel):
				if np.any(tableFunctionsOverSpan[j, :]==i):
					minspan = j; break
			maxspan = self._nbel
			for j in range(self._nbel-1, -1, -1):
				if np.any(tableFunctionsOverSpan[j, :]==i):
					maxspan = j; break
			tableSpansOverFunction[i, :] = [minspan, maxspan]

		for i in range(0, self.nbctrlpts):
			minspan, maxspan = tableSpansOverFunction[i, :]
			minknot = tablePointsOverSpan[minspan, 0] + 1
			maxknot = tablePointsOverSpan[maxspan, 1] - 1
			if i == 0: minknot -= 1
			if i == self.nbctrlpts-1: maxknot += 1
			B0shape[i, :] = [minknot, maxknot]

		for i in range(0, self.nbctrlpts):
			minspan, maxspan = tableSpansOverFunction[i, :]
			minknot = tablePointsOverSpan[minspan, 0] + 1
			maxknot = tablePointsOverSpan[maxspan, 1] - 1
			if (i == 0) or (i == 1): minknot -= 1
			if (i == self.nbctrlpts-1) or (i == self.nbctrlpts-2): maxknot += 1
			B1shape[i, :] = [minknot, maxknot]
		
		if isfortran: 
			B0shape += 1; B1shape += 1 # change from 0 to 1 index
		
		return [B0shape, B1shape]
	
	def evalDersBasisWeights(self):	
		if self._wqType not in [1, 2]: raise Warning('Method unknown')
		size_data = (self.degree + 1)*self.nbqp
		basis   = np.zeros((size_data, 2))
		weights = np.zeros((size_data, 4))
		indj = np.zeros(size_data, dtype=int)
		indi = np.zeros(self.nbctrlpts+1, dtype=int)

		if (self._posRule == 'midpoint') and (self._isUniform) and (self._nbel > self.degree+3):
			s = self._extraArgs.get('s', 1); r = self._extraArgs.get('r', 2)
			# Create model
			degree_model = self.degree
			kv_model     = createUniformMaxregularKnotvector(degree_model, degree_model + 3)
			kwargs       = {'type': self._wqType, 'extra': self._extraArgs}
			WQmodel      = WeightedQuadrature(degree_model, kv_model, quadArgs=kwargs)
			WQmodel.__findQuadraturePositions()
			size_data_model = (degree_model + 1)*WQmodel.nbqp
			Bm, Wm, indim, indjm = basisweights.wq_getbasisweights_csr(WQmodel.degree, WQmodel.knotvector, WQmodel.quadPtsPos, 
													WQmodel._Bshape[0], WQmodel._Bshape[1], size_data_model, WQmodel._wqType)
			# Scale the results
			Wm[:, :2] = Wm[:, :2]*WQmodel._nbel/self._nbel
			Bm[:, -1] = Bm[:, -1]*self._nbel/WQmodel._nbel

			#Copy model data
			times = self._nbel - (self.degree + 3); row2copy = self.degree + 2
			left  = indim[row2copy-1] - 1; right = indim[row2copy] - 1
			indj2copy = indjm[left:right];  data2copy = np.zeros((len(indj2copy), 6)) 
			data2copy[:, :2] = Bm[left:right, :]; data2copy[:, 2:] = Wm[left:right, :]

			indi = np.copy(indim); indj = np.copy(indjm)
			data = np.hstack((Bm, Wm))
			for i in range(1, times+1):
				row2copyt  = row2copy + i
				indj2copyt = indj2copy + i*(s+1)

				indit = indi; indjt = indj; datat = data
				inditt, indjtt, datatt = insertRowCSR(row2copyt, data2copy, indj2copyt, indit, indjt, datat)
				indi = np.copy(inditt); indj = np.copy(indjtt); data = np.copy(datatt)

			# Offset of last p+1 rows
			nbcols = 2*(self.degree + r) + self._nbel*(s + 1) - 2*s - 3  
			offset = nbcols - np.max(indjm)

			for i in range(self.nbctrlpts-self.degree-1, self.nbctrlpts):
				indj[indi[i]:indi[i+1]] += offset

			basis = data[:, :2]; weights = data[:, 2:]

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
		return self.quadPtsPos, self.dersIndices, self.dersBasis, self.dersWeights
	