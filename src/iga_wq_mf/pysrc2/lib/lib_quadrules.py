from lib.__init__ import *
from lib.lib_base import gaussTable, lobattoTable, findMultiplicity, createKnotVector, insertRowCSR

class QuadratureRules:
	def __init__(self, degree, knotvector):
		self._degree      = degree
		self._knotvector  = knotvector
		self._nbqp        = None
		self._quadPtsPos  = None
		self._dersIndices = None
		self._dersBasis   = None
		self._dersWeights = None
		self.getInfoFromKnotvector()
		self.verifyUniformityRegularity()
		return

	def getInfoFromKnotvector(self):
		self._nbctrlpts = len(self._knotvector) - self._degree - 1
		self._uniqueKV  = np.unique(self._knotvector)
		self._nbel      = len(self._uniqueKV) - 1
		return
	
	def verifyUniformityRegularity(self, threshold=1e-8):
		self._isMaxReg = False
		if (self._nbel+self._degree==self._nbctrlpts): self._isMaxReg = True

		self._isUniform = True
		diffknot = np.diff(np.diff(self._uniqueKV))
		if np.any(np.abs(diffknot)>=threshold): self._isUniform = False
		return
	
	def getQuadratureRulesInfo(self):
		return self._quadPtsPos, self._dersIndices, self._dersBasis, self._dersWeights
	
	def getDenseQuadRules(self, isFortran=True):
		self._denseBasis   = []
		self._denseWeights = []
		indi, indj = self._dersIndices
		for i in range(np.size(self._dersBasis, axis=1)):
			if isFortran: tmp = sp.csr_matrix((self._dersBasis[:, i], indj-1, indi-1))
			else: tmp = sp.csr_matrix((self._dersBasis[:, i], indj, indi))
			self._denseBasis.append(tmp)

		for i in range(np.size(self._dersWeights, axis=1)):
			if isFortran: tmp = sp.csr_matrix((self._dersWeights[:, i], indj-1, indi-1))
			else: tmp = sp.csr_matrix((self._dersWeights[:, i], indj, indi))
			self._denseWeights.append(tmp)
		return
	
class GaussQuadrature(QuadratureRules):
	def __init__(self, degree, knotvector, kwargs={}):
		super().__init__(degree, knotvector)
		self._kwargs = kwargs
		self._quadMethod = kwargs.get('quadmethod', 'leg').lower()
		if self._quadMethod == 'leg': # Legendre 
			self._order = self._degree + 1
			self._table = gaussTable
		elif self._quadMethod == 'lob': # Lobatto
			self._order = self._degree + 2
			self._table = lobattoTable
		else: 
			raise Warning('Not found')
		return
	
	def getGaussInfo(self):
		" Gets Gauss quadrature points position and weights from tables "
		gp, gw = self._table(self._order)
		self._isoPositions = gp
		self._isoWeights   = gw
		return

	def findQuadraturePositions(self):
		quadPtsInfo = np.array([])
		knots = self._uniqueKV
		for i in range(self._nbel):
			xg = 0.5*((knots[i+1] - knots[i])*self._isoPositions + knots[i] + knots[i+1])
			quadPtsInfo = np.append(quadPtsInfo, xg)
		self._quadPtsPos = np.atleast_1d(quadPtsInfo)
		self._nbqp = np.size(self._quadPtsPos)
		return
	
	def findParametricWeights(self):
		quadPtsInfo = np.array([])
		knots = self._uniqueKV
		for i in range(self._nbel):
			wg = 0.5*(knots[i+1] - knots[i])*self._isoWeights
			quadPtsInfo = np.append(quadPtsInfo, wg)
		self._parametricWeights = np.atleast_1d(quadPtsInfo)
		return
	
	def evalDersBasisWeights(self):
		B, indi, indj = basis_weights.get_genbasis_csr(self._degree, self._knotvector, self._quadPtsPos)
		self._dersBasis   = B
		self._dersIndices = [indi, indj]
		nnz = np.shape(B)[0]
		W = np.zeros((nnz, 4))
		for i in range(nnz):
			W[i, 0]  = B[i, 0]*self._parametricWeights[indj[i]-1] # -1 because indj is in fortran notation
			W[i, 3]  = B[i, 1]*self._parametricWeights[indj[i]-1]
		W[:, 1] = W[:, 0]; W[:, 2] = W[:, 3]
		self._dersWeights = W
		return
	
	def getQuadratureRulesInfo(self):
		self.getGaussInfo()
		self.findQuadraturePositions()
		self.findParametricWeights()
		self.evalDersBasisWeights()
		return self._quadPtsPos, self._dersIndices, self._dersBasis, self._dersWeights

class WeightedQuadrature(QuadratureRules):
	def __init__(self, degree, knotvector, kwargs={}):
		super().__init__(degree, knotvector)
		self._kwargs     = kwargs
		self._quadMethod = kwargs.get('quadmethod', 1)
		return
	
	def findQuadraturePositions(self):
		rule = self._kwargs.get('rule', 'midpoint')
		self._posRule = rule.lower()
		if self._posRule == 'midpoint':
			s = self._kwargs.get('s', 1); r = self._kwargs.get('r', 2)
			self._quadPtsPos = self.QuadPosMidPointRule(s=s, r=r)
		# Add more methods 
		self._Bshape = self.getBShape(self._quadPtsPos)
		self._nbqp   = np.size(self._quadPtsPos)
		return

	# Add rules to get the quadrature points
	def QuadPosMidPointRule(self, s=1, r=2):	
		""" Find quadrature points using mid-point rule
			r is the 'extra' points on the knotspans at the boundaries
			s is the 'extra' points on the inner knotspans. The number must respect the max rule (see bibliography)
		"""
		quadPtsInfo = np.array([])
		knots = self._uniqueKV
		# First span
		tmp = np.linspace(knots[0], knots[1], self._degree+r)
		quadPtsInfo = np.append(quadPtsInfo, tmp)

		# Last span
		tmp = np.linspace(knots[-2], knots[-1], self._degree+r)
		quadPtsInfo = np.append(quadPtsInfo, tmp)

		# Inner span
		for i in range(1, self._nbel-1):
			tmp = np.linspace(knots[i], knots[i+1], 2+s)
			quadPtsInfo = np.append(quadPtsInfo, tmp)
		
		quadPtsInfo = np.unique(quadPtsInfo)
		quadPtsInfo.sort()
		return quadPtsInfo
	
	# ----

	def getBShape(self, knots, isfortran=True):
		""" Return the shape of basis in WQ approach.
		"""
		B0shape = np.zeros((self._nbctrlpts, 2), dtype=int)
		B1shape = np.zeros((self._nbctrlpts, 2), dtype=int)
		knots   = np.atleast_1d(knots)

		tablePointsOverSpan = np.zeros((self._nbel, 2), dtype=int)
		for i in range(0, self._nbel):
			left = self._uniqueKV[i]; right = self._uniqueKV[i+1]
			boolean = (knots>=left)*(knots<=right)
			tablePointsOverSpan[i, 0] = np.nonzero(boolean)[0][0]
			tablePointsOverSpan[i, 1] = np.nonzero(boolean)[0][-1]

		tableFunctionsOverSpan = np.zeros((self._nbel, self._degree+1), dtype=int)
		for j in range(0, self._degree+1): tableFunctionsOverSpan[0, j] = j
		for i in range(1, self._nbel): 
			multiplicity = findMultiplicity(self._degree, self._knotvector, self._uniqueKV[i])
			tableFunctionsOverSpan[i, 0] = tableFunctionsOverSpan[i-1, 0] + multiplicity
			for j in range(1, self._degree+1): 
				tableFunctionsOverSpan[i, j] = tableFunctionsOverSpan[i, 0] + j

		tableSpansOverFunction = np.zeros((self._nbctrlpts, 2), dtype=int)
		for i in range(0, self._nbctrlpts):
			minspan = 1
			for j in range(0, self._nbel):
				if np.any(tableFunctionsOverSpan[j, :]==i):
					minspan = j; break
			maxspan = self._nbel
			for j in range(self._nbel-1, -1, -1):
				if np.any(tableFunctionsOverSpan[j, :]==i):
					maxspan = j; break
			tableSpansOverFunction[i, :] = [minspan, maxspan]

		for i in range(0, self._nbctrlpts):
			minspan, maxspan = tableSpansOverFunction[i, :]
			minknot = tablePointsOverSpan[minspan, 0] + 1
			maxknot = tablePointsOverSpan[maxspan, 1] - 1
			if i == 0: minknot -= 1
			if i == self._nbctrlpts-1: maxknot += 1
			B0shape[i, :] = [minknot, maxknot]

		for i in range(0, self._nbctrlpts):
			minspan, maxspan = tableSpansOverFunction[i, :]
			minknot = tablePointsOverSpan[minspan, 0] + 1
			maxknot = tablePointsOverSpan[maxspan, 1] - 1
			if (i == 0) or (i == 1): minknot -= 1
			if (i == self._nbctrlpts-1) or (i == self._nbctrlpts-2): maxknot += 1
			B1shape[i, :] = [minknot, maxknot]
		
		if isfortran: 
			B0shape += 1; B1shape += 1 # change from 0 to 1 index
		
		return [B0shape, B1shape]
	
	def evalDersBasisWeights(self):
		if self._quadMethod not in [1, 2]: raise Warning('Method unknown')
		
		size_data = (self._degree + 1)*self._nbqp
		B = np.zeros((size_data, 2))
		W = np.zeros((size_data, 4))
		indj = np.zeros(size_data, dtype=int)
		indi = np.zeros(self._nbctrlpts+1, dtype=int)

		if (self._posRule == 'midpoint') and (self._isUniform) and (self._nbel > self._degree+3):
			s = self._kwargs.get('s', 1); r = self._kwargs.get('r', 2)
			# Create model
			degree_model = self._degree
			kv_model     = createKnotVector(degree_model, degree_model + 3)
			WQmodel      = WeightedQuadrature(degree_model, kv_model, kwargs=self._kwargs)
			WQmodel.findQuadraturePositions()
			size_data_model = (degree_model + 1)*WQmodel._nbqp
			Bm, Wm, indim, indjm = basis_weights.wq_getbasisweights_csr(WQmodel._degree, WQmodel._knotvector, WQmodel._quadPtsPos, 
													WQmodel._Bshape[0], WQmodel._Bshape[1], size_data_model, WQmodel._quadMethod)
			# Scale the results
			Wm[:, :2] = Wm[:, :2]*WQmodel._nbel/self._nbel
			Bm[:, -1] = Bm[:, -1]*self._nbel/WQmodel._nbel

			#Copy model data
			times = self._nbel - (self._degree + 3); row2copy = self._degree + 2
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
			nbcols = 2*(self._degree + r) + self._nbel*(s + 1) - 2*s - 3  
			offset = nbcols - np.max(indjm)

			for i in range(self._nbctrlpts-self._degree-1, self._nbctrlpts):
				indj[indi[i]:indi[i+1]] += offset

			B = data[:, :2]; W = data[:, 2:]

		else:
			B, W, indi, indj = basis_weights.wq_getbasisweights_csr(self._degree, self._knotvector, self._quadPtsPos, 
														self._Bshape[0], self._Bshape[1], size_data, self._quadMethod)
		
		self._dersBasis   = B
		self._dersWeights = W
		self._dersIndices = [indi, indj]
		return
	
	def getQuadratureRulesInfo(self):
		self.findQuadraturePositions()
		self.evalDersBasisWeights()
		return self._quadPtsPos, self._dersIndices, self._dersBasis, self._dersWeights
	