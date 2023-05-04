"""
.. This module contains functions that will be used in other modules. 
.. Some of them use fortran functions to accelerate computing time
.. Joaquin Cornejo
"""

from .__init__ import *

# ==========================
# GENERAL FUNCTIONS
# ==========================

def cropImage(filename):
	from PIL import Image
	im = Image.open(filename).convert('RGB')
	na = np.array(im)
	colorY, colorX = np.where(np.all(na!=[255, 255, 255], axis=2))

	# Find first and last row containing colored pixels
	top, bottom = colorY[0], colorY[-1]

	# Extract region of Interest
	ROI = na[top:bottom, :]
	Image.fromarray(ROI).save(filename)
	return 

def relativeError(array_interp, array_th):
	error = array_th - array_interp
	try:    relError = np.linalg.norm(error, np.inf)/np.linalg.norm(array_th, np.inf)
	except: relError = sp.linalg.norm(error, np.inf)/sp.linalg.norm(array_th, np.inf)
	return relError

def sigmoid(x, c1=1, c2=0):
	f = 1.0/(1.0 + np.exp(-c1*(x - c2)))
	return f

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
	if lenOfNewRow != len(indj_in): raise Warning('Size problem')
	if (row2in > len(indi)) or (row2in < 1): raise Warning('Not possible')

	data_out = np.zeros((nnz+lenOfNewRow, ncols))
	indj_out = np.zeros(nnz+lenOfNewRow, dtype=int)
	indi_out = np.ones(len(indi)+1, dtype=int)

	left  = indi[row2in-1] - 1
	indj_out[0:left] = indj[0:left];               data_out[0:left, :] = data[0:left, :]
	indj_out[left:left+lenOfNewRow] = indj_in;     data_out[left:left+lenOfNewRow, :] = data_in
	indj_out[left+lenOfNewRow:]   = indj[left:]; data_out[left+lenOfNewRow:, :]   = data[left:, :]

	left  = row2in
	indi_out[0:left] = indi[0:left]
	for i in range(left, len(indi_out)): indi_out[i] = indi[i-1] + lenOfNewRow

	return indi_out, indj_out, data_out

def array2csr_matrix(data, indi, indj, isfortran=True):
	" Computes csr sparse matrix "

	nb_rows = len(indi) - 1
	if isfortran: 
		nb_cols  = max(indj)
		indjcopy = indj - 1
		indicopy = indi - 1 
	else: 
		nb_cols = max(indj) + 1 
		indjcopy = np.copy(indj)
		indicopy = np.copy(indi)
	sparse_matrix = sp.csr_matrix((data, indjcopy, indicopy), shape=(nb_rows, nb_cols))
									
	return sparse_matrix

# ==========================
# B-SPLINE FUNCTIONS
# ==========================

def createKnotVector(p, nbel, multiplicity=1):
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

def findMultiplicity(p, knotvector, knot, threshold=1e-8):
	""" Finds the multiplicity of a given knot.
		Ex: Given the knot-vector {0, 0, 0, 0.5, 0.5, 1, 1, 1} and x = 0.5, the multiplicity is 2.
	"""
	multiplicity = 0
	for i in range(0, len(knotvector)):
		if (np.abs(knot-knotvector[i])<=threshold):
			multiplicity += 1

	return multiplicity

def evalDersBasisPy(degree, knotvector, knots, multiplicity=1): 
	""" Evaluates B-spline functions at given knots. 
		Knot-vector needs to be regular
	"""

	nbknots = len(knots)
	nbel    = len(np.unique(knotvector)) - 1
	nb_ctrlpts = len(knotvector) - degree - 1

	B0 = sp.lil_matrix((nb_ctrlpts, nbknots))
	B1 = sp.lil_matrix((nb_ctrlpts, nbknots))

	# Set table of functions per element 
	table_functions_physpan = np.zeros((nbel, degree + 2), dtype=int); 
	table_functions_physpan[0, 0] = degree; table_functions_physpan[0, 1:] = np.arange(degree + 1) 

	for _ in range(1, nbel): 
		table_functions_physpan[_, :2] = table_functions_physpan[_-1, :2] + multiplicity
		table_functions_physpan[_, 2:] = table_functions_physpan[_, 1] + np.arange(1, degree + 1) 

	for i, knot in enumerate(knots):
		knot_span = helpers.find_span_linear(degree, knotvector, nb_ctrlpts, knot)
		phy_span  = np.where(table_functions_physpan[:, 0] == knot_span)[0].tolist()
		functions_span = table_functions_physpan[phy_span, 1:][0]
		B0t, B1t = helpers.basis_function_ders(degree, knotvector, knot_span, knot, 1)

		# Set procedure if knot is in the knot-vector
		if knot in np.unique(knotvector)[1:-1]:             
			B0t = B0t[:-multiplicity] 
			B1t = B1t[:-multiplicity] 
			functions_span = functions_span[:-multiplicity]

		B0[np.ix_(functions_span, [i])] = np.asarray(B0t).reshape((-1, 1))
		B1[np.ix_(functions_span, [i])] = np.asarray(B1t).reshape((-1, 1))

	B0, B1 = B0.tocsr(), B1.tocsr()

	return B0, B1

def evalDersBasisFortran(degree, knotvector, knots):
	" Evaluates B-spline functions at given knots using fortran libraries "
	B, indi, indj = basisweights.get_genbasis_csr(degree, knotvector, knots)
	return B, indi, indj

def gaussTable(order):
	" Computes Gauss weights and positions in isoparametric space for a given degree "
	if order <= 10: pos, wgt = basisweights.gauss_quadrature_table(order)
	else: raise Warning('Not degree found')
	return pos, wgt

def lobattoTable(order):
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
	else: raise Warning('Not defined')
	return pos, wgt

# =========================
# MF FUNCTIONS
# =========================

def createTableProperties(function, uref=None, prop=None):
	"Create a table of scalar properties from given function "

	# Set default given u
	if uref is None: uref = np.linspace(-1., 1., 201)

	# Compute v
	if prop is None: v = function(uref)
	else:  v = function(uref, prop)

	# Create table 
	table = np.zeros((len(uref), 2))
	table[:, 0] = uref; table[:, 1] = v

	return table

def scipySolver(A, b, nbIterations=100, epsilon=1e-10, PreCond='ilu', isCG=True):
	""" Solves system using an iterative method : conjugate gradient or bi-conjugate gradient. 
		It can use ILU preconditioner to accelerate convergence. 
	"""

	# Select preconditionner (by the moment only ILU)
	if PreCond == 'ilu': 
		B = sclin.spilu(A)
		Mx = lambda x: B.solve(x)
		M = sclin.LinearOperator(A.shape, Mx)

	# Solve with iterative method
	if isCG: x, info = sclin.cg(A, b, tol=epsilon, maxiter=nbIterations, M=M)
	else: x, info = sclin.bicgstab(A, b, tol=epsilon, maxiter=nbIterations, M=M)

	return x

def genEigenDecomposition(indi, indj, data, robin_condition=[0, 0], coefs=None):
	""" Eigen decomposition generalized KU = MUD
		K: stiffness matrix, K = int B1 B1 dx = W11 * B1
		M: mass matrix, M = int B0 B0 dx = W00 * B0
		U: eigenvectors matrix
		D: diagonal of eigenvalues
	"""

	nc = np.max(indj)
	[B0, B1, W0, W1] = data

	# Create mcoef and kcoef is are none type
	if coefs is None: mcoefs = np.ones(nc); kcoefs = np.ones(nc)
	else: [mcoefs, kcoefs] = coefs

	eigenvalues, eigenvectors = geophy.eigen_decomposition_py(indi, indj, B0, W0, 
										B1, W1, mcoefs, kcoefs, robin_condition)

	return eigenvalues, eigenvectors

def separationOfVariables3d(n_list, CC_matrix, dim=3):
	""" Tensor decomposition of CC matrix to improve Fast diagonalization precontionner
		Based on "Preconditioners for Isogemetric Analysis" by M. Montardini
	"""
	if dim != 3: raise Warning('Method not implemented')

	# Unpack shape
	[n1, n2, n3] = n_list
	try: n1 = n1[0]; n2 = n2[0]; n3 = n3[0]
	except: pass

	# Get diagonal of CC matrix
	coefs = np.zeros((dim, n1*n2*n3))
	for i in range(dim): coefs[i, :] = CC_matrix[i, i, :]

	# Initialize
	u1, u2, u3 = np.ones(n1), np.ones(n2), np.ones(n3)
	w1, w2, w3 = np.ones(n1), np.ones(n2), np.ones(n3)
	Vscript = np.zeros((n1, n2, n3))
	Wscript = np.zeros((2, n1, n2, n3))
	Nscript = np.zeros((n1, n2, n3))
	Mscript = np.zeros((n1, n2, n3))

	# Transform coefs to tensor of coefs
	coefstens = np.zeros((dim, n1, n2, n3))
	for k in range(dim):
		for i3 in range(n3):
			for i2 in range(n2):
				for i1 in range(n1):
					pos = i1 + i2*n1 + i3*n1*n2
					coefstens[k, i1, i2, i3] = coefs[k, pos]

	for iter in range(2):
		for k in range(dim):
			# Set Dscript
			for i3 in range(n3):
				for i2 in range(n2):
					for i1 in range(n1):
						U = [u1[i1], u2[i2], u3[i3]]
						Vscript[i1, i2, i3] = coefstens[k, i1, i2, i3]*U[k]\
												/(U[0]*U[1]*U[2])

			# Update w
			if k == 0:
				for j in range(n1):
					m = Vscript[j, :, :].min()
					M = Vscript[j, :, :].max()
					w1[j] = np.sqrt(m*M)
			elif k == 1: 
				for j in range(n2):
					m = Vscript[:, j, :].min()
					M = Vscript[:, j, :].max()
					w2[j] = np.sqrt(m*M)        
			elif k == 2: 
				for j in range(n3):
					m = Vscript[:, :, j].min()
					M = Vscript[:, :, j].max()
					w3[j] = np.sqrt(m*M)

		for k in range(dim):
			c = -1
			# Compute Wscript temporary
			for l in range(dim):
				if k != l:
					c += 1
					for i3 in range(n3):
						for i2 in range(n2):
							for i1 in range(n1):
								U = [u1[i1], u2[i2], u3[i3]]
								W = [w1[i1], w2[i2], w3[i3]]
								Wscript[c, i1, i2, i3] = coefstens[k, i1, i2, i3]*U[l]*U[k]\
															/(U[0]*U[1]*U[2]*W[k])

			# Compute Nscript and Mscript
			for i3 in range(n3):
				for i2 in range(n2):
					for i1 in range(n1): 
						WWlk = [Wscript[_, i1, i2, i3] for _ in range(2)]
						Nscript[i1, i2, i3] = min(WWlk)
						Mscript[i1, i2, i3] = max(WWlk)

			# Update u
			if k == 0:
				for j in range(n1):
					m = Nscript[j, :, :].min()
					M = Mscript[j, :, :].max()
					u1[j] = np.sqrt(m*M)
			elif k == 1: 
				for j in range(n2):
					m = Nscript[:, j, :].min()
					M = Mscript[:, j, :].max()
					u2[j] = np.sqrt(m*M)        
			elif k == 2: 
				for j in range(n3):
					m = Nscript[:, :, j].min()
					M = Mscript[:, :, j].max()
					u3[j] = np.sqrt(m*M) 

	return u1, u2, u3, w1, w2, w3

def computeEigenDiag(eig_u, eig_v, eig_w, coefs=[1.0, 1.0, 1.0]):
	" Computes diagonal of eigen values "

	def kron3vec(A, B, C):
		" Computes kron product of 3 vectors: A x B x C"
		result = np.kron(np.kron(A, B), C)
		return result

	one_u = np.ones(len(eig_u))
	one_v = np.ones(len(eig_v))
	one_w = np.ones(len(eig_w))

	eig_diag  = coefs[0]*kron3vec(one_w, one_v, eig_u)
	eig_diag += coefs[1]*kron3vec(one_w, eig_v, one_u)
	eig_diag += coefs[2]*kron3vec(eig_w, one_v, one_u)
	
	return eig_diag

def fastDiagonalization(U, V, W, D, array_in, fdtype='steady'):
	" Compute fast diagonalization using Fortran"
	if fdtype == 'interp':
		array_out = heatsolver.fd_interpolation_3d(U, V, W, array_in)
	elif fdtype == 'steady':
		array_out = heatsolver.fd_steady_heat_3d(U, V, W, D, array_in)
	elif fdtype == 'elastic':
		array_out = plasticitysolver.fd_elasticity_3d(U, V, W, D, array_in)
	
	return array_out

def computeMean3d(ncu, ncv, ncw, C):

	def trapezoidal_rule_3d(coefs):
		nru, nrv, nrw = np.shape(coefs)
		indu = [0, nru-1]; indv = [0, nrv-1]; indw = [0, nrw-1]
		integral = 0.0
		for iw in range(1, nrw-1):
			for iv in range(1, nrv-1):
				for iu in range(1, nru-1):
					integral += coefs[iu, iv, iw]

		for iw in range(1, nrw-1):
			for iv in range(1, nrv-1):
				for iu in [0, 1]:
					integral += coefs[indu[iu], iv, iw]/2.0

		for iw in range(1, nrw-1):
			for iv in [0, 1]:
				for iu in range(1, nru-1):
					integral += coefs[iu, indv[iv], iw]/2.0

		for iw in [0, 1]:
			for iv in range(1, nrv-1):
				for iu in range(1, nru-1):
					integral += coefs[iu, iv, indw[iw]]/2.0

		for iw in [0, 1]:
			for iv in [0, 1]:
				for iu in range(1, nru-1):
					integral += coefs[iu, indv[iv], indw[iw]]/4.0

		for iw in [0, 1]:
			for iv in range(1, nrv-1):
				for iu in [0, 1]:
					integral += coefs[indu[iu], iv, indw[iw]]/4.0

		for iw in range(1, nrw-1):
			for iv in [0, 1]:
				for iu in [0, 1]:
					integral += coefs[indu[iu], indv[iv], iw]/4.0
	
		for iw in [0, 1]:
			for iv in [0, 1]:
				for iu in [0, 1]:
					integral += coefs[indu[iu], indv[iv], indw[iw]]/8.0

		integral = integral/((nru-1)*(nrv-1)*(nrw-1))

		return integral

	Cset = np.zeros((3, 3, 3))
	if ncu*ncv*ncw != len(C): raise Warning('Mismatch dimension')
	pos = int((ncu-1)/2); indu = [0, pos, ncu-1]
	pos = int((ncv-1)/2); indv = [0, pos, ncv-1]
	pos = int((ncw-1)/2); indw = [0, pos, ncw-1]
	
	for k in range(3):
		for j in range(3):
			for i in range(3):
				genPos = indu[i] + indv[j]*ncu + indw[k]*ncu*ncv
				Cset[i, j, k] = C[genPos]

	integral = trapezoidal_rule_3d(Cset)

	return integral

def CG(funcAu, b, nbIterPCG=100, threshold=1e-8):   

	x = np.zeros(np.shape(b)); r = b
	resPCG = np.zeros(nbIterPCG+1); resPCG[0] = 1.0
	normb  = np.amax(np.absolute(r))

	if normb <= threshold: return
	rsold = np.dot(r, r); p = r

	for iter in range(nbIterPCG):
		Ap = funcAu(p)
		alpha = rsold/np.dot(p, Ap)
		x = x + alpha*p
		r = r - alpha*Ap

		resPCG[iter+1] = np.amax(np.absolute(r))/normb
		if resPCG[iter+1] <= threshold: break
		
		rsnew = np.dot(r, r)
		p = r + rsnew/rsold * p
		rsold = rsnew

	return x, resPCG

def PCG(funcAu, funcPre, b, nbIterPCG=100, threshold=1e-8):   

	x = np.zeros(np.shape(b)); r = b
	resPCG = np.zeros(nbIterPCG+1); resPCG[0] = 1.0
	normb  = np.amax(np.absolute(r))

	if normb <= threshold: return
	z = funcPre(r)
	rsold = np.dot(r, r); p = z

	for iter in range(nbIterPCG):
		Ap = funcAu(p)
		alpha = rsold/np.dot(p, Ap)
		x = x + alpha*p
		r = r - alpha*Ap

		resPCG[iter+1] = np.amax(np.absolute(r))/normb
		if resPCG[iter+1] <= threshold: break

		z = funcPre(r)
		rsnew = np.dot(r, z)
		p = z + rsnew/rsold * p
		rsold = rsnew

	return x, resPCG