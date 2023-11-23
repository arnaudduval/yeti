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
		Knot-vector needs to be regular
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
		pos = [-1.0, -0.8997579954114602, -0.6771862795107377,
				-0.3631174638261782, 0.0, 0.3631174638261782,
				0.6771862795107377, 0.8997579954114602, 1.0]
		wgt = [0.0277777777777778, 0.1654953615608055, 0.2745387125001617,
				0.3464285109730463, 0.3715192743764172, 0.3464285109730463,
				0.2745387125001617, 0.1654953615608055, 0.0277777777777778]
	elif order == 10:
		pos = [-1.0, -0.9195339081664589, -0.7387738651055050,
				-0.4779249498104445, -0.1652789576663870, 0.1652789576663870,
				0.4779249498104445, 0.7387738651055050, 0.9195339081664589, 1.0]
		wgt = [0.0222222222222222, 0.1333059908510701, 0.2248893420631264,
				0.2920426836796838, 0.3275397611838976, 0.3275397611838976,
				0.2920426836796838, 0.2248893420631264, 0.1333059908510701,
				0.0222222222222222]
	elif order == 11:
		pos = [-1.0, -0.9340014304080592, -0.7844834736631444,
				-0.5652353269962050, -0.2957581355869394, 0.0,
				0.2957581355869394, 0.5652353269962050, 0.7844834736631444,
				0.9340014304080592, 1.0]
		wgt = [0.0181818181818182, 0.1096122732669949, 0.1871698817803052,
				0.2480481042640284, 0.2868791247790080, 0.3002175954556907,
				0.2868791247790080, 0.2480481042640284, 0.1871698817803052,
				0.1096122732669949, 0.0181818181818182]
	elif order == 12:
		pos = [-1.0, -0.9448992722228822, -0.8192793216440067,
				-0.6328761530318606, -0.3995309409653489, -0.1365529328549276,
				0.1365529328549276, 0.3995309409653489, 0.6328761530318606,
				0.8192793216440067, 0.9448992722228822, 1.0]
		wgt = [0.0151515151515152, 0.0916845174131962, 0.1579747055643701,
				0.2125084177610211, 0.2512756031992013, 0.2714052409106962,
				0.2714052409106962, 0.2512756031992013, 0.2125084177610211,
				0.1579747055643701, 0.0916845174131962, 0.0151515151515152]
	else: raise Warning('Not defined')
	pos = np.array(pos); wgt = np.array(wgt)
	return pos, wgt

# =========================
# OTHERS
# =========================

class solver():
	
	def __init__(self):
		self._nbIter       = 100
		self._thresholdPCG = 1e-12
		return
	
	def CG(self, Afun, b):
		x = np.zeros(np.shape(b))
		r = b; normb = np.linalg.norm(r)
		resPCG = [1.0]
		if normb <= self._thresholdPCG: return
		rsold = np.dot(r, r); p = r

		for i in range(self._nbIter):
			Ap = Afun(p)
			alpha = rsold/np.dot(p, Ap)
			x += alpha*p
			r -= alpha*Ap

			resPCG.append(np.linalg.norm(r)/normb)
			if (resPCG[-1]<=self._thresholdPCG): break

			rsnew = np.dot(r, r)
			p = r + rsnew/rsold*p
			rsold = rsnew

		return x, resPCG
	
	def PCG(self, Afun, Pfun, b):
		x = np.zeros(np.shape(b))
		r = b; normb = np.linalg.norm(r)
		resPCG = [1.0]
		if normb <= self._thresholdPCG: return
		z = Pfun(r)
		rsold = np.dot(r, z); p = z

		for i in range(self._nbIter):
			Ap = Afun(p)
			alpha = rsold/np.dot(p, Ap)
			x += alpha*p
			r -= alpha*Ap

			resPCG.append(np.linalg.norm(r)/normb)
			if (resPCG[-1]<=self._thresholdPCG): break

			z = Pfun(r)
			rsnew = np.dot(r, z)
			p = z + rsnew/rsold*p
			rsold = rsnew

		return x, resPCG
	
	def BiCGSTAB(self, Afun, b):
		x = np.zeros(np.shape(b))
		r = b; normb = np.linalg.norm(r)
		resPCG = [1.0]
		if normb <= self._thresholdPCG: return
		rhat = r; p = r
		rsold = np.dot(r, rhat)

		for i in range(self._nbIter):
			Ap = Afun(p)
			alpha = rsold/np.dot(Ap, rhat)
			s = r - alpha*Ap
			As = Afun(s)
			omega = np.dot(As, s)/np.dot(As, As)
			x += alpha*p + omega*s
			r = s - omega*As

			resPCG.append(np.linalg.norm(r)/normb)
			if (resPCG[-1]<=self._thresholdPCG): break

			rsnew = np.dot(r, rhat)
			beta = (alpha/omega)*(rsnew/rsold)
			p = r + beta(p - omega*Ap)
			rsold = rsnew

		return x, resPCG
	
	def PBiCGSTAB(self, Afun, Pfun, b):
		x = np.zeros(np.shape(b))
		r = b; normb = np.linalg.norm(r)
		resPCG = [1.0]
		if normb <= self._thresholdPCG: return
		rhat = r; p = r
		rsold = np.dot(r, rhat)

		for i in range(self._nbIter):
			ptilde = Pfun(p)
			Aptilde = Afun(ptilde)
			alpha = rsold/np.dot(Aptilde, rhat)
			s = r - alpha*Aptilde
			stilde = Pfun(s)
			Astilde = Afun(stilde)
			omega = np.dot(Astilde, s)/np.dot(Astilde, Astilde)
			x += alpha*ptilde + omega*stilde
			r = s - omega*Astilde

			resPCG.append(np.linalg.norm(r)/normb)
			if (resPCG[-1]<=self._thresholdPCG): break

			rsnew = np.dot(r, rhat)
			beta = (alpha/omega)*(rsnew/rsold)
			p = r + beta*(p - omega*Aptilde)
			rsold = rsnew

		return x, resPCG
	