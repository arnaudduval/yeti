"""
.. This module contains functions that will be used in other modules. 
.. Some of them use fortran functions to accelerate computing time
.. Joaquin Cornejo
"""

from .__init__ import *

# ==========================
# GENERAL FUNCTIONS
# ==========================
def sigmoid(x, c1=1, c2=0):
	f = 1.0/(1.0 + np.exp(-c1*(x - c2)))
	return f

def cropImage(filename):
	from PIL import Image
	im = Image.open(filename).convert('RGB')
	na = np.array(im)
	colorY, colorX = np.where(np.all(na!=[255, 255, 255], axis=2))

	# Find first and last row containing colored pixels
	top, bottom = min(colorY), max(colorY)
	left, right = min(colorX), max(colorX)

	# Extract region of Interest
	ROI = na[top:bottom, left:right]
	Image.fromarray(ROI).save(filename)
	return 

def eraseRowsCSR(rows2er, indi_in, indj_in, data_in, isfortran=True):
	" Returns new data after erasing rows in CSR format "
	
	indi_int, indj_int = np.copy(indi_in), np.copy(indj_in) 
	if isfortran: indi_int -= 1; indj_int -= 1 
	indi_outt = np.delete(indi_int, rows2er)
	indj_out, data_out = [], []
	
	# Copy column indices
	for i in range(len(indi_outt)-1): 
		indj_out.extend(indj_int[indi_outt[i]:indi_outt[i+1]])
	indj_out = np.array(indj_out, dtype=int)

	# Copy data 
	for a_in in data_in: 
		a_in = np.atleast_2d(a_in); a_out = []
		for i in range(len(indi_outt)-1): 
			a_out.extend(a_in[indi_outt[i]:indi_outt[i+1], :])
		a_out = np.array(a_out)
		data_out.append(a_out)

	# Update row indices
	indi_out = [0]
	for i in range(len(indi_outt)-1): 
		nnz = indi_outt[i+1] - indi_outt[i]
		newvalue = indi_out[-1] + nnz
		indi_out.append(newvalue)
	indi_out = np.array(indi_out, dtype=int)

	if isfortran: indi_out += 1; indj_out += 1 

	return indi_out, indj_out, data_out

def insertRowCSR(row2in, data_in, indj_in, indi, indj, data):
	""" Returns the data and indices after inserting a new row (CSR format).
		Indices must start at 1, like fortran style.
	"""
	nnz = np.size(data, axis=0); ncols = np.size(data, axis=1)
	lenOfNewRow = np.size(data_in, axis=0) 
	assert lenOfNewRow == len(indj_in), 'Size problem'
	assert (row2in < len(indi)) and (row2in > 1), 'This is not a possible configuration'

	data_out = np.zeros((nnz+lenOfNewRow, ncols))
	indj_out = np.zeros(nnz+lenOfNewRow, dtype=int)
	indi_out = np.ones(len(indi)+1, dtype=int)

	left  = indi[row2in-1] - 1
	indj_out[0:left] = indj[0:left];           data_out[0:left, :] = data[0:left, :]
	indj_out[left:left+lenOfNewRow] = indj_in; data_out[left:left+lenOfNewRow, :] = data_in
	indj_out[left+lenOfNewRow:] = indj[left:]; data_out[left+lenOfNewRow:, :]   = data[left:, :]

	left  = row2in
	indi_out[0:left] = indi[0:left]
	for i in range(left, len(indi_out)): indi_out[i] = indi[i-1] + lenOfNewRow

	return indi_out, indj_out, data_out

def array2csr_matrix(data, indi, indj, isFortran=True):
	indjcopy = np.copy(indj)
	indicopy = np.copy(indi)
	if isFortran: 
		indjcopy = indjcopy - 1
		indicopy = indicopy - 1 
	sparse_matrix = sp.csr_matrix((data, indjcopy, indicopy))
	return sparse_matrix

def get_faceInfo(nb):
	dir  = int(np.floor(nb/2))
	side = 0
	if nb%2 == 1: side = 1
	return dir, side

def get_INCTable(nnzByDimension):
	" Sets topology table, also known as INC: NURBS coordinates. "
	# Create INC: NURBS coordinates
	nnz_total = np.prod(nnzByDimension)
	table = np.zeros((nnz_total, max([3, len(nnzByDimension)])), dtype= int)
	if len(nnzByDimension) == 3:
		for i3 in range(nnzByDimension[2]): 
			for i2 in range(nnzByDimension[1]): 
				for i1 in range(nnzByDimension[0]):
					genPos = i1 + i2*nnzByDimension[0] + i3*nnzByDimension[0]*nnzByDimension[1]
					table[genPos, :] = [i1, i2, i3]
	elif len(nnzByDimension) == 4:
		for i4 in range(nnzByDimension[3]): 
			for i3 in range(nnzByDimension[2]): 
				for i2 in range(nnzByDimension[1]): 
					for i1 in range(nnzByDimension[0]):
						genPos = i1 + i2*nnzByDimension[0] + i3*nnzByDimension[0]*nnzByDimension[1] + i4*nnzByDimension[0]*nnzByDimension[1]*nnzByDimension[2]
						table[genPos, :] = [i1, i2, i3, i4]

	else: raise Warning('Try other array')

	return table

# ==========================
# B-SPLINE FUNCTIONS
# ==========================

def createUniformKnotvector_Rmultiplicity(p, nbel, multiplicity=1):
	" Creates an uniform and open knot-vector with a given regularity "

	kv_unique = np.linspace(0., 1., nbel + 1)[1 : -1]

	knotvector = []
	for _ in range(p+1): 
		knotvector.append(0.0)

	for knot in kv_unique: 
		for _ in range(multiplicity): 
			knotvector.append(knot)

	for _ in range(p+1): 
		knotvector.append(1.0)

	knotvector = np.array(knotvector)
	
	return knotvector

def createAsymmetricalKnotvector_Rmultiplicity(p, nbel, xasym=0.25, multiplicity=1):
	assert nbel > 1, 'Not possible. Try higher number of elements'
	if nbel == 2: createUniformKnotvector_Rmultiplicity(p, nbel, multiplicity=multiplicity)
	else:
		nbel1 = int(np.floor((nbel+1)/2))
		kv_unique = np.array([])
		kv_unique = np.append(kv_unique, np.linspace(0.0, xasym, nbel1 + 1)[1:]) 
		kv_unique = np.append(kv_unique, np.linspace(xasym, 1.0, nbel - nbel1 + 1)[1:-1]) 
		knotvector = []
		for _ in range(p+1): 
			knotvector.append(0.0)

		for knot in kv_unique: 
			for _ in range(multiplicity): 
				knotvector.append(knot)

		for _ in range(p+1): 
			knotvector.append(1.0)

		knotvector = np.array(knotvector)
	return knotvector

def createUniformCurve(p, nbel, length):
	knotvector = createUniformKnotvector_Rmultiplicity(p, 1)
	ctrlpts    = [[i*length/p, 0.0] for i in range(p+1)]
	crv = BSpline.Curve()
	crv.degree  = p
	crv.ctrlpts = ctrlpts
	crv.knotvector = knotvector
	for knot in np.linspace(0, 1, nbel+1)[1:-1]:
		operations.insert_knot(crv, [knot], [1])
	return crv

def createAsymmetricalCurve(p, nbel, length, xasym=0.25):
	knotvector0 = createUniformKnotvector_Rmultiplicity(p, 1)
	ctrlpts    = [[i*length/p, 0.0] for i in range(p+1)]
	crv = BSpline.Curve()
	crv.degree  = p
	crv.ctrlpts = ctrlpts
	crv.knotvector = knotvector0
	knotvector1 = createAsymmetricalKnotvector_Rmultiplicity(p, nbel, xasym=xasym)
	for knot in knotvector1[p+1:-p-1]:
		operations.insert_knot(crv, [knot], [1])
	return crv

def findInterpolationSpan(array, x, threshold=1e-8):
	span = 1
	while (span<len(array)-1 and (array[span]-x<=threshold)):
		span += 1
	span -= 1
	return span

def findMultiplicity(knotvector, knot, threshold=1e-8):
	""" Finds the multiplicity of a given knot.
		Ex: Given the knot-vector {0, 0, 0, 0.5, 0.5, 1, 1, 1} and x = 0.5, the multiplicity is 2.
	"""
	multiplicity = 0
	for i in range(0, len(knotvector)):
		if (np.abs(knot-knotvector[i])<=threshold):
			multiplicity += 1

	return multiplicity

def evalDersBasisPy(degree, knotvector, knots): 
	""" Evaluates B-spline functions at given knots. 
		Knot-vector needs to be regular.
		THIS FUNCTION IS DEPRECATED
	"""

	nbknots   = len(knots)
	nbctrlpts = len(knotvector) - degree - 1
	uvk  = np.unique(knotvector)
	nbel = len(uvk) - 1 

	basis, indices = np.zeros(((degree+1)*nbknots, 2)), np.zeros(((degree+1)*nbknots, 2), dtype=int)
	# Set table of functions per element 
	table_functions_physpan = np.zeros((nbel, degree + 1), dtype=int); 
	table_functions_physpan[0, :] = np.arange(degree + 1) 

	for i in range(1, nbel):
		multiplicity = findMultiplicity(knotvector, uvk[i])
		table_functions_physpan[i, 0] = table_functions_physpan[i-1, 0] + multiplicity
		table_functions_physpan[i, 1:] = table_functions_physpan[i, 0] + np.arange(1, degree + 1) 

	k = 0
	for i, knot in enumerate(knots):
		knot_span = helpers.find_span_linear(degree, knotvector, nbctrlpts, knot)
		phy_span  = findInterpolationSpan(uvk, knot)
		functions_span = table_functions_physpan[phy_span, :]
		B0t, B1t = helpers.basis_function_ders(degree, knotvector, knot_span, knot, 1)

		for j in range(degree + 1):
			basis[k, :]   = [B0t[j], B1t[j]]
			indices[k, :] = [functions_span[j], i]
			k += 1

	B0 = sp.coo_matrix((basis[:, 0], (indices[:, 0], indices[:, 1])), shape=(nbctrlpts, nbknots))
	B1 = sp.coo_matrix((basis[:, 1], (indices[:, 0], indices[:, 1])), shape=(nbctrlpts, nbknots))
	B0, B1 = B0.tocsr(), B1.tocsr()

	return B0, B1

def evalDersBasisFortran(degree, knotvector, knots):
	" Evaluates B-spline functions at given knots using fortran libraries "
	B, indi, indj = basisweights.get_genbasis_csr(degree, knotvector, knots)
	return B, indi, indj

def legendreTable(order):
	" Computes Gauss-Legendre weights and positions in isoparametric space for a given degree "
	assert order <= 10, 'Only degrees below 10'
	pos, wgt = basisweights.gauss_quadrature_table(order)
	return pos, wgt

def lobattoTable(order):
	" Computes Gauss-Lobatto weights and positions in isoparametric space for a given degree "
	if order == 2:
		pos = [-1.0, 1.0]
		wgt = [1.0, 1.0]
	elif order == 3:
		pos = [-1.0, 0.0, 1.0]
		wgt = [1.0/3.0, 4.0/3.0, 1.0/3.0]
	elif order == 4:
		pos = [-1.0,
			-0.447_213_595_499_958,
			0.447_213_595_499_958,
			1.0]
		wgt = [1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0]
	elif order == 5:
		pos = [-1.0,
			-0.654_653_670_707_977,
			0.0,
			0.654_653_670_707_977,
			1.0]
		wgt = [0.1, 4.9/9.0, 6.4/9.0, 4.9/9.0, 0.1]
	elif order == 6:
		pos = [-1.0,
			-0.765_055_323_929_465,
			-0.285_231_516_480_645,
			0.285_231_516_480_645,
			0.765_055_323_929_465,
			1.0]
		wgt = [0.066_666_666_666_667,
			0.378_474_956_297_847,
			0.554_858_377_035_486,
			0.554_858_377_035_486,
			0.378_474_956_297_847,
			0.066_666_666_666_667]
	elif order == 7:
		pos = [-1.0, -0.830_223_896_278_5670, -0.468_848_793_470_7142, 0.0,
			0.468_848_793_470_7142, 0.830_223_896_278_5670, 1.0]
		wgt = [0.047_619_047_619_0476, 0.276_826_047_361_5659, 0.431_745_381_209_8627,
				0.487_619_047_619_0476, 0.431_745_381_209_8627, 0.276_826_047_361_5659,
				0.047_619_047_619_0476]
	elif order == 8:
		pos = [-1.0, -0.871_740_148_509_6066, -0.591_700_181_433_1423,
				-0.209_299_217_902_4789, 0.209_299_217_902_4789, 0.591_700_181_433_1423,
				0.871_740_148_509_6066, 1.0]
		wgt = [0.035_714_285_714_2857, 0.210_704_227_143_5061, 0.341_122_692_483_5044,
				0.412_458_794_658_7038, 0.412_458_794_658_7038, 0.341_122_692_483_5044,
				0.210_704_227_143_5061, 0.035_714_285_714_2857]
	elif order == 9:
		pos = [-1.0, -0.899_757_995_411_4602, -0.677_186_279_510_7377,
				-0.363_117_463_826_1782, 0.0, 0.363_117_463_826_1782,
				0.677_186_279_510_7377, 0.899_757_995_411_4602, 1.0]
		wgt = [0.027_777_777_777_7778, 0.165_495_361_560_8055, 0.274_538_712_500_1617,
				0.346_428_510_973_0463, 0.371_519_274_376_4172, 0.346_428_510_973_0463,
				0.274_538_712_500_1617, 0.165_495_361_560_8055, 0.027_777_777_777_7778]
	elif order == 10:
		pos = [-1.0, -0.919_533_908_166_4589, -0.738_773_865_105_5050,
				-0.477_924_949_810_4445, -0.165_278_957_666_3870, 0.165_278_957_666_3870,
				0.477_924_949_810_4445, 0.738_773_865_105_5050, 0.919_533_908_166_4589, 1.0]
		wgt = [0.022_222_222_222_2222, 0.133_305_990_851_0701, 0.224_889_342_063_1264,
				0.292_042_683_679_6838, 0.327_539_761_183_8976, 0.327_539_761_183_8976,
				0.292_042_683_679_6838, 0.224_889_342_063_1264, 0.133_305_990_851_0701,
				0.022_222_222_222_2222]
	elif order == 11:
		pos = [-1.0, -0.934_001_430_408_0592, -0.784_483_473_663_1444,-0.565_235_326_996_2050, 
				-0.295_758_135_586_9394, 0.0, 0.295_758_135_586_9394, 0.565_235_326_996_2050, 
				0.784_483_473_663_1444, 0.934_001_430_408_0592, 1.0]
		wgt = [0.018_181_818_181_8182, 0.109_612_273_266_9949, 0.187_169_881_780_3052,
				0.248_048_104_264_0284, 0.286_879_124_779_0080, 0.300_217_595_455_6907,
				0.286_879_124_779_0080, 0.248_048_104_264_0284, 0.187_169_881_780_3052,
				0.109_612_273_266_9949, 0.018_181_818_181_8182]
	elif order == 12:
		pos = [-1.0, -0.944_899_272_222_8822, -0.819_279_321_644_0067,
				-0.632_876_153_031_8606, -0.399_530_940_965_3489, -0.136_552_932_854_9276,
				0.136_552_932_854_9276, 0.399_530_940_965_3489, 0.632_876_153_031_8606,
				0.819_279_321_644_0067, 0.944_899_272_222_8822, 1.0]
		wgt = [0.015_151_515_151_5152, 0.091_684_517_413_1962, 0.157_974_705_564_3701,
				0.212_508_417_761_0211, 0.251_275_603_199_2013, 0.271_405_240_910_6962,
				0.271_405_240_910_6962, 0.251_275_603_199_2013, 0.212_508_417_761_0211,
				0.157_974_705_564_3701, 0.091_684_517_413_1962, 0.015_151_515_151_5152]
	else: raise Warning('Not defined')
	pos = np.array(pos); wgt = np.array(wgt)
	return pos, wgt

# =========================
# OTHERS
# =========================

class solver():
	
	def __init__(self):
		self._itersLin = 100
		self._thresLin = 1e-12
		return
	
	def eigs(self, N, Afun, Bfun=None, neigvals=100, which='SM', allowcomplex=False):
		""" 
		Computes the eigenvalues of the linear system A x = lambda x or A x = lambda B x 
		by an iterative solver (with a matrix free approach).
		By default, we compute the smallest eigevalues.
		"""
		ALinOp = sp.linalg.LinearOperator((N, N), matvec=Afun)
		if Bfun is not None: BLinOp = sp.linalg.LinearOperator((N, N), matvec=Bfun)
		if Bfun is None: eigvals, eigvecs = sp.linalg.eigs(A=ALinOp, k=neigvals, which=which)
		else: eigvals, eigvecs = sp.linalg.eigs(A=ALinOp, k=neigvals, B=BLinOp, which=which)
		if not allowcomplex: eigvals = np.absolute(eigvals)
		return eigvals, eigvecs
	
	def CG(self, Afun, b, Pfun=None, dotfun=None, cleanfun=None, dod=None):
		if Pfun is None: Pfun = lambda x: x
		if dotfun is None: dotfun = lambda x, y: np.dot(x, y)
		if cleanfun is None: cleanfun = lambda x, y: x

		x = np.zeros(np.shape(b))
		r = np.copy(b); cleanfun(r, dod)
		normb = np.sqrt(dotfun(r, r))
		resLin = [1.0]
		if normb <= self._thresLin: return
		z = Pfun(r); cleanfun(z, dod)
		p = np.copy(z)
		rsold = dotfun(r, z)

		for i in range(self._itersLin):
			Ap = Afun(p); cleanfun(Ap, dod)
			alpha = rsold/dotfun(p, Ap)
			x += alpha*p
			r -= alpha*Ap

			resLin.append(np.sqrt(dotfun(r, r))/normb)
			if (resLin[-1]<=self._thresLin): break

			z = Pfun(r); cleanfun(z, dod)
			rsnew = dotfun(r, z)
			p = z + rsnew/rsold*p
			rsold = np.copy(rsnew)
		output = {'sol':x, 'res':resLin}
		return output
	
	def BiCGSTAB(self, Afun, b, Pfun=None, dotfun=None, cleanfun=None, dod=None):
		if Pfun is None: Pfun = lambda x: x
		if dotfun is None: dotfun = lambda x, y: np.dot(x, y)
		if cleanfun is None: cleanfun = lambda x, y: x

		x = np.zeros(np.shape(b))
		r = np.copy(b); cleanfun(r, dod)
		normb = np.sqrt(dotfun(r, r))
		resLin = [1.0]
		if normb <= self._thresLin: return
		rhat = np.copy(r)
		p = np.copy(r)
		rsold = dotfun(r, rhat)

		for i in range(self._itersLin):
			ptilde = Pfun(p); cleanfun(ptilde, dod)
			Aptilde = Afun(ptilde); cleanfun(Aptilde, dod)
			alpha = rsold/dotfun(Aptilde, rhat)
			s = r - alpha*Aptilde

			stilde = Pfun(s); cleanfun(stilde, dod)
			Astilde = Afun(stilde); cleanfun(Astilde, dod)
			omega = dotfun(Astilde, s)/dotfun(Astilde, Astilde)
			x += alpha*ptilde + omega*stilde
			r = s - omega*Astilde

			resLin.append(np.sqrt(dotfun(r, r))/normb)
			if (resLin[-1]<=self._thresLin): break

			rsnew = dotfun(r, rhat)
			beta = (alpha/omega)*(rsnew/rsold)
			p = r + beta*(p - omega*Aptilde)
			rsold = np.copy(rsnew)

		output = {'sol':x, 'res':resLin}
		return output
	
	def GMRES(self, Afun, b, Pfun=None, dotfun=None, cleanfun=None, dod=None, n_restarts=1):
		if Pfun is None: Pfun = lambda x: x
		if dotfun is None: dotfun = lambda x, y: np.dot(x, y)
		if cleanfun is None: cleanfun = lambda x, y: x
    
		x = np.zeros(*np.shape(b))
		H = np.zeros((self._itersLin + 1, self._itersLin))
		V = np.zeros((self._itersLin + 1, *np.shape(b)))
		Z = np.zeros((self._itersLin + 1, *np.shape(b)))
		normb = np.zeros(n_restarts)
		
		AllresLin = []
		for m in range(n_restarts):
			r = b - Afun(x); cleanfun(r, dod)
			normb[m] = np.sqrt(dotfun(r, r))
			if normb[m] <= self._thresLin*normb[0]: break
			V[0] = r/normb[m]
			e1 = np.zeros(self._itersLin + 1); e1[0] = normb[m]
			
			resLin = [1.0]
			for k in range(self._itersLin):
				Z[k] = Pfun(V[k]); cleanfun(Z[k], dod)
				w = Afun(Z[k]); cleanfun(w, dod)
				for j in range(k+1):
					H[j, k] = dotfun(w, V[j])
					w -= H[j, k] * V[j]
				H[j+1, j] = np.sqrt(dotfun(w, w))
				if H[j+1, j] != 0: V[k+1] = w/H[j+1, j]
				y = np.linalg.lstsq(H[:k+2, :k+1], e1[:k+2])[0]
				zeta = np.linalg.norm(H[:k+2, :k+1] @ y - e1[:k+2])
				resLin.append(zeta/normb[0])
				if resLin[-1] <= self._thresLin: break
				
			for j in range(k+1): x = x + Z[j]*y[j]
			AllresLin.append(resLin)
		output = {'sol':x, 'res':resLin}
		return x, output