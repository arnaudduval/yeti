"""
.. This module contains functions that will be used in other modules. 
.. Some of them use fortran functions to accelerate computing time
.. Joaquin Cornejo
"""

from .__init__ import *

# ==========================
# GENERAL FUNCTIONS
# ==========================

def relativeError(array_app, array_th):
	error = array_th - array_app
	try:    relError = np.linalg.norm(error, np.inf)/np.linalg.norm(array_th, np.inf)
	except: relError = sp.linalg.norm(error, np.inf)/sp.linalg.norm(array_th, np.inf)
	return relError

def sigmoid(x, c1=1, c2=0):
	f = 1.0/(1.0 + np.exp(-c1*(x-c2)))
	return f

def erase_rows_csr(rows2er, indi_in, indj_in, data_in, isfortran=True):
	" Returns new data after erasing rows in CSR format "
	
	indi_int, indj_int = np.copy(indi_in), np.copy(indj_in) 
	if isfortran: indi_int -= 1; indj_int -= 1 
	indi_outt = np.delete(indi_int, rows2er)
	indj_out, data_out = [], []
	
	# Copy column indices
	for _ in range(len(indi_outt)-1): 
		indj_out.extend(indj_int[indi_outt[_]:indi_outt[_+1]])
	indj_out = np.array(indj_out, dtype=int)

	# Copy data 
	for a_in in data_in: 
		a_in = np.atleast_2d(a_in); a_out = []
		for _ in range(len(indi_outt)-1): 
			a_out.extend(a_in[indi_outt[_]:indi_outt[_+1], :])
		a_out = np.array(a_out)
		data_out.append(a_out)

	# Update row indices
	indi_out = [0]
	for _ in range(len(indi_outt)-1): 
		nnz = indi_outt[_+1] - indi_outt[_]
		newvalue = indi_out[-1] + nnz
		indi_out.append(newvalue)
	indi_out = np.array(indi_out, dtype=int)

	if isfortran: indi_out += 1; indj_out += 1 

	return indi_out, indj_out, data_out

# ==========================
# B-SPLINE FUNCTIONS
# ==========================

def create_knotvector(p, nbel, multiplicity=1):
	" Creates an uniform, open, with maximum regularity, knot-vector "

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

def eval_basis_python(degree, knotvector, knots, multiplicity=1): 
	""" Evaluates B-spline functions at given knots. 
		Knot-vector needs to be regular
	"""

	nbknots = len(knots)
	nbel = len(np.unique(knotvector)) - 1
	nb_ctrlpts = degree + multiplicity*(nbel - 1) + 1

	# Set table of functions per element 
	table_functions_element = np.zeros((nbel, degree + 2), dtype=int); 
	table_functions_element[0, 0] = degree; table_functions_element[0, 1:] = np.arange(degree + 1) 

	for _ in range(1, nbel): 
		table_functions_element[_, :2] = table_functions_element[_-1, :2] + multiplicity
		table_functions_element[_, 2:] = table_functions_element[_, 1] + np.arange(1, degree + 1) 

	# Set B0 and B1
	B0 = sp.lil_matrix((nb_ctrlpts, nbknots))
	B1 = sp.lil_matrix((nb_ctrlpts, nbknots))

	for i, knot in enumerate(knots):
	
		knot_span = helpers.find_span_linear(degree, knotvector, nb_ctrlpts, knot)
		element = np.where(table_functions_element[:, 0] == knot_span)[0].tolist()
		functions_element = table_functions_element[element, 1:][0]
		B0t, B1t = helpers.basis_function_ders(degree, knotvector, knot_span, knot, 1)

		# Set procedure if knot is in the knot-vector
		if knot in np.unique(knotvector)[1:-1]:             
			B0t = B0t[:-multiplicity] 
			B1t = B1t[:-multiplicity] 
			functions_element = functions_element[:-multiplicity]

		B0[np.ix_(functions_element, [i])] = np.asarray(B0t).reshape((-1,1))
		B1[np.ix_(functions_element, [i])] = np.asarray(B1t).reshape((-1,1))

	B0, B1 = B0.tocsr(), B1.tocsr()

	return B0, B1

def eval_basis_fortran(degree, knotvector, knots):
	" Evaluates B-spline functions at given knots using fortran libraries "

	B, indi, indj = basis_weights.get_basis_generalized_csr(degree, knotvector, knots)
	
	return B, indi, indj

# ==========================
# IGA FUNCTIONS
# ==========================

def gaussTable(order):
	" Computes Gauss weights and positions in isoparametric space for a given degree "

	if order == 1:
		pos = [0.0]
		wgt = [2.0]
	elif order == 2:
		pos = [ -0.57735026918962576,
				0.57735026918962576]

		wgt = [ 1.0,
				1.0]
	elif order == 3:
		pos = [ -0.77459666924148337,
				0.0,
				0.77459666924148337]

		wgt = [ 5.0 / 9.0,
				8.0 / 9.0,
				5.0 / 9.0]
	elif order == 4:
		pos = [ -0.8611363115940526,
				-0.3399810435848563,
				0.3399810435848563, 
				0.8611363115940526]

		wgt = [ 0.3478548451374539,
				0.6521451548625461,
				0.6521451548625461,
				0.3478548451374539]
	elif order == 5:
		pos = [ -0.9061798459386640,
				-0.5384693101056831,
				0.0,
				0.5384693101056831,
				0.9061798459386640]

		wgt = [ 0.2369268850561891,
				0.4786286704993665,
				0.5688888888888889,
				0.4786286704993665,
				0.2369268850561891]
	elif order == 6:
		pos = [ -0.9324695142031520,
				-0.6612093864662645,
				-0.2386191860831969,
				0.2386191860831969,
				0.6612093864662645,
				0.9324695142031520]

		wgt = [ 0.1713244923791703,
				0.3607615730481386, 
				0.4679139345726910,
				0.4679139345726910,
				0.3607615730481386,
				0.1713244923791703]
	elif order == 7:
		pos = [ -0.9491079123427585,
				-0.7415311855993944,
				-0.4058451513773972,
				0.0,
				0.4058451513773972,
				0.7415311855993944,
				0.9491079123427585]

		wgt = [ 0.1294849661688697,
				0.2797053914892767,
				0.3818300505051189,
				0.4179591836734694,
				0.3818300505051189,
				0.2797053914892767,
				0.1294849661688697]
	elif order == 8:
		pos = [ -0.9602898564975362,
				-0.7966664774136267,
				-0.5255324099163290,
				-0.1834346424956498,
				0.1834346424956498,
				0.5255324099163290,
				0.7966664774136267,
				0.9602898564975362]

		wgt = [ 0.1012285362903763,
				0.2223810344533745,
				0.3137066458778873,
				0.3626837833783620,
				0.3626837833783620,
				0.3137066458778873,
				0.2223810344533745,
				0.1012285362903763]
	elif order == 9:
		pos = [ -0.9681602395076261,
				-0.8360311073266358,
				-0.6133714327005904,
				-0.3242534234038089,
				0.0,
				0.3242534234038089,
				0.6133714327005904,
				0.8360311073266358,
				0.9681602395076261]

		wgt = [ 0.0812743883615744,
				0.1806481606948574,
				0.2606106964029354,
				0.3123470770400028,
				0.3302393550012597,
				0.3123470770400028,
				0.2606106964029354,
				0.1806481606948574,
				0.0812743883615744]
	elif order == 10:
		pos = [ -0.9739065285171717,
				-0.8650633666889845,
				-0.6794095682990244,
				-0.4333953941292472,
				-0.1488743389816312,
				0.1488743389816312,
				0.4333953941292472,
				0.6794095682990244,
				0.8650633666889845,
				0.9739065285171717]

		wgt = [ 0.0666713443086881,
				0.1494513491505806,
				0.2190863625159820,
				0.2692667193099963,
				0.2955242247147529,
				0.2955242247147529,
				0.2692667193099963,
				0.2190863625159820,
				0.1494513491505806,
				0.0666713443086881]                
	elif order == 11:
		pos = [ -0.9782286581460570,
				-0.8870625997680953,
				-0.7301520055740494,
				-0.5190961292068118,
				-0.2695431559523450,
				0.0,
				0.2695431559523450,
				0.5190961292068118,
				0.7301520055740494,
				0.8870625997680953,
				0.9782286581460570]

		wgt = [ 0.0556685671161737,
				0.1255803694649046,
				0.1862902109277343,
				0.2331937645919905,
				0.2628045445102467,
				0.2729250867779006,
				0.2628045445102467,
				0.2331937645919905,
				0.1862902109277343,
				0.1255803694649046,
				0.0556685671161737]

	elif order == 12:
		pos = [ -0.9815606342467192,
				-0.9041172563704749,
				-0.7699026741943047,
				-0.5873179542866175,
				-0.3678314989981802,
				-0.1252334085114689,
				0.1252334085114689,
				0.3678314989981802,
				0.5873179542866175,
				0.7699026741943047,
				0.9041172563704749,
				0.9815606342467192]

		wgt = [ 0.0471753363865118,
				0.1069393259953184,
				0.1600783285433462,
				0.2031674267230659,
				0.2334925365383548,
				0.2491470458134028,
				0.2491470458134028,
				0.2334925365383548,
				0.2031674267230659,
				0.1600783285433462,
				0.1069393259953184,
				0.0471753363865118]

	elif order == 15:
		pos = [ -0.9879925180204854,
				-0.9372733924007060,
				-0.8482065834104272,
				-0.7244177313601701,
				-0.5709721726085388,
				-0.3941513470775634,
				-0.2011940939974345,
				0.0,
				0.2011940939974345,
				0.3941513470775634,
				0.5709721726085388,
				0.7244177313601701,
				0.8482065834104272,
				0.9372733924007060,
				0.9879925180204854]

		wgt = [ 0.0307532419961173,
				0.0703660474881081,
				0.1071592204671719,
				0.1395706779261543,
				0.1662692058169939,
				0.1861610000155622,
				0.1984314853271116,
				0.2025782419255613,
				0.1984314853271116,	
				0.1861610000155622,
				0.1662692058169939,
				0.1395706779261543,
				0.1071592204671719,
				0.0703660474881081,
				0.0307532419961173]

	else: raise Warning('Not degree found')
	
	# Change type of arrays
	pos, wgt = np.array(pos), np.array(wgt)

	return pos, wgt

def iga_find_positions_weights(degree, knotvector):
	" Computes Gauss positions and weights in parametric space for a given knotvector "

	def iga_find_positions_weights_element(degree, knotvector, el):
		" Computes Gauss weights and positions in parametric space for a given element "

		knots = np.unique(knotvector)
		xg, wg = gaussTable(degree + 1)
		xg = 0.5*((knots[el+1] - knots[el])*xg + knots[el] + knots[el+1])
		wg = 0.5*(knots[el+1] - knots[el])*wg

		return xg, wg

	nbel = len(np.unique(knotvector)) - 1

	# Find quadrature points and weights
	xg = []; wg = []
	for _ in range(nbel): 
		xg_el, wg_el = iga_find_positions_weights_element(degree, knotvector, _)
		xg.extend(xg_el); wg.extend(wg_el)

	xg, wg = np.array(xg), np.array(wg)

	return xg, wg

def iga_find_basis_weights_opt(degree, knotvector):
	" Computes basis and weights in IGA approach "

	qp_position, qp_weight = iga_find_positions_weights(degree, knotvector)
	B0, B1 = eval_basis_python(degree, knotvector, qp_position)

	return qp_position, B0, B1, qp_weight

def iga_find_basis_weights_fortran(degree, knotvector): 
	" Computes basis and weights in IGA approach using fortran "

	# Set number of elements, quadrature points and size of data
	nbel = len(np.unique(knotvector)) - 1
	nb_qp = (degree + 1) * nbel
	nnz_B = (degree + 1) * nb_qp

	# Get basis and weights from fortran
	qp_position, qp_weight, basis, \
	indi, indj, nnz_I = basis_weights.iga_get_data_csr(degree, knotvector, nb_qp, nnz_B)

	return nnz_I, qp_position, qp_weight, basis, indi, indj

# =========================
# WQ FUNCTIONS
# =========================

def wq_get_shape_B(degree, nbel, r, multiplicity=1, maxrule=1): 
	""" Return the shape of basis in WQ approach 
		It considers an open and uniform knot-vector
	"""

	nb_ctrlpts = degree + multiplicity*(nbel - 1) + 1
	nb_qp = 2*(degree + r) + nbel*(maxrule + 1) - 2*maxrule - 3

	# Set table of positions of quadrature points 
	tableOfPointsOnSpan = np.zeros((nbel, 2), dtype=int)
	tableOfPointsOnSpan[0, 0] = 0; tableOfPointsOnSpan[0, 1] = degree + r -1 
	tableOfPointsOnSpan[-1, 0] = nb_qp - (degree + r) ; tableOfPointsOnSpan[-1,1] = nb_qp - 1

	for i in range(1, nbel-1):
		tableOfPointsOnSpan[i, 0] = tableOfPointsOnSpan[i - 1, 1]
		tableOfPointsOnSpan[i, 1] = tableOfPointsOnSpan[i, 0] + 1 + maxrule
	
	# Set table of functions on every element 
	tableOfFunctionsOnElement = np.zeros((nbel, degree + 1), dtype=int); 
	tableOfFunctionsOnElement[0, :] = np.arange(degree + 1) 

	for _ in range(1, nbel): 
		tableOfFunctionsOnElement[_, 0] = tableOfFunctionsOnElement[_-1, 0] + multiplicity
		tableOfFunctionsOnElement[_, 1:] = tableOfFunctionsOnElement[_, 0] + np.arange(1, degree + 1) 

	# Set table of knot-span for every function
	tableOfSpansForFunction = []
	for i in range(nb_ctrlpts):
		indi, _ = np.asarray(tableOfFunctionsOnElement == i).nonzero()
		indi = np.sort(indi)
		tableOfSpansForFunction.append(indi)

	# Find basis shape
	indi_B0, indj_B0 = [], []
	indi_B1, indj_B1 = [], []
	for i in range(nb_ctrlpts) : 
		spanOfFunction = tableOfSpansForFunction[i]
		spanLeft = np.min(spanOfFunction)
		spanRight = np.max(spanOfFunction)

		support = np.arange(tableOfPointsOnSpan[spanLeft,0], \
							tableOfPointsOnSpan[spanRight,1] + 1, dtype=int)
		
		support_B0 = support
		if i == 0: 
			support_B0 = np.delete(support, [-1])
		elif i == nb_ctrlpts - 1:
			support_B0 = np.delete(support, [0])
		else : 
			support_B0 = np.delete(support, [0 , -1])
		
		indi_B0.extend(i*np.ones(len(support_B0), dtype= int))
		indj_B0.extend(support_B0)

		support_B1 = support
		if i == 0 or i == 1: 
			support_B1 = np.delete(support, [-1])
		elif i == nb_ctrlpts - 1 or i == nb_ctrlpts - 2:
			support_B1 = np.delete(support, [0])
		else : 
			support_B1 = np.delete(support, [0 , -1])
		
		indi_B1.extend(i*np.ones(len(support_B1), dtype= int))
		indj_B1.extend(support_B1)
	
	data_B0 = np.ones(len(indi_B0))
	B0shape = sp.csr_matrix((data_B0, (indi_B0, indj_B0)), shape=(nb_ctrlpts, nb_qp))
	
	data_B1 = np.ones(len(indi_B1))
	B1shape = sp.csr_matrix((data_B1, (indi_B1, indj_B1)), shape=(nb_ctrlpts, nb_qp))
	
	return B0shape, B1shape

def wq_solve_equation_system(B, I):
	""" Solves linear system B w = I where B and I are written in CSR format. 
		System can be under, well, over determined
	"""

	# Convert type of B and I to array
	B = B.toarray(); I = I.toarray()

	# Solve system 
	sol1 = np.linalg.lstsq(B, I, rcond=None)[0]
	w = sol1.reshape((1, -1)).tolist()[0]
	
	return w

def wq_find_positions(degree, knotvector, r, maxrule= 1):
	" Return position of quadrature points in WQ approach "
	
	kv_unique = np.unique(knotvector)
	nbel = len(kv_unique) - 1

	# Compute quadrature positions
	qp_position = []
	for i in range(nbel):

		if i == 0 or i == nbel - 1:
			xt = np.linspace(kv_unique[i], kv_unique[i+1], degree+r)[1:-1]
		else:
			xt = np.linspace(kv_unique[i], kv_unique[i+1], 2 + maxrule)[1:-1]
		
		qp_position.extend(xt)
	
	# Include knots of knot-vector
	qp_position.extend(list(kv_unique))

	# Change type and sort
	qp_position = np.array(qp_position)
	qp_position.sort()
	
	return qp_position

def wq_find_weights(degree, knotvector, r):
	" Returns weights at quadrature points in WQ approach using space S^[p-1] method "
	
	nbel = len(np.unique(knotvector))-1
	nb_ctrlpts = degree + nbel

	# ------------
	# Space S^p_r
	# ------------
	# Find positions and weights in IGA approach 
	xcgg, wcgg = iga_find_positions_weights(degree, knotvector)

	# Find basis function values at Gauss points 
	B0cgg_p0, B1cgg_p0 = eval_basis_python(degree, knotvector, xcgg)

	# Find positions in WQ approach
	xwq = wq_find_positions(degree, knotvector, r)

	# Find basis function values at WQ points
	B0wq_p0 = eval_basis_python(degree, knotvector, xwq)[0]

	# ---------------------
	# Sapce S^[p-1]_[r-1]
	# ---------------------
	# Create knot-vector
	knotvector_p1 = knotvector[1:-1]

	# Find basis function values at Gauss points
	B0cgg_p1 = eval_basis_python(degree-1, knotvector_p1, xcgg)[0]

	# Find basis function values at WQ points
	B0wq_p1 = eval_basis_python(degree-1, knotvector_p1, xwq)[0]  

	# ------------------
	# Compute Integrals
	# ------------------
	I00 = B0cgg_p0 @ sp.diags(wcgg) @ B0cgg_p0.T
	I01 = B0cgg_p1 @ sp.diags(wcgg) @ B0cgg_p0.T 
	I10 = B0cgg_p0 @ sp.diags(wcgg) @ B1cgg_p0.T
	I11 = B0cgg_p1 @ sp.diags(wcgg) @ B1cgg_p0.T    

	# Set W00, W01, W10 and W11
	data_W00, data_W01, data_W10, data_W11 = [], [], [], []
	
	# Find shape of B and I
	nb_qp_wq = len(xwq)
	shape_B0_p0, shape_B1_p0 = wq_get_shape_B(degree, nbel, r)
	shape_B_p1 = wq_get_shape_B(degree-1, nbel, r+1)[0]

	# -------------------
	# Computation of W00
	# -------------------
	Ishape = shape_B0_p0 @ shape_B0_p0.T 
	indi = []; indj = []
	for i in range(nb_ctrlpts):
		F = np.nonzero(Ishape[:, i])[0]
		P = np.nonzero(shape_B0_p0[i, :])[1]

		Bt = B0wq_p0[np.ix_(F,P)]
		It = I00[np.ix_(F, [i])]

		data_W00.extend(wq_solve_equation_system(Bt, It))
		indi.extend(i*np.ones(len(P), dtype= int))
		indj.extend(P)
	W00 = sp.csr_matrix((data_W00, (indi, indj)), shape=(nb_ctrlpts, nb_qp_wq))

	# --------------------
	# Computation of W01
	# --------------------
	Ishape = shape_B_p1 @ shape_B0_p0.T 
	indi = []; indj = []
	for i in range(nb_ctrlpts):
		F = np.nonzero(Ishape[:,i])[0]
		P = np.nonzero(shape_B0_p0[i, :])[1]

		Bt = B0wq_p1[np.ix_(F,P)]
		It = I01[np.ix_(F, [i])]
		
		data_W01.extend(wq_solve_equation_system(Bt, It))
		indi.extend(i*np.ones(len(P), dtype= int))
		indj.extend(P)
	W01 = sp.csr_matrix((data_W01, (indi, indj)), shape=(nb_ctrlpts, nb_qp_wq))
	
	# --------------------
	# Computation of W10
	# --------------------
	Ishape = shape_B1_p0 @ shape_B0_p0.T 
	indi = []; indj = []
	for i in range(nb_ctrlpts):
		F = np.nonzero(Ishape[:, i])[0]
		P = np.nonzero(shape_B1_p0[i, :])[1]

		Bt = B0wq_p0[np.ix_(F, P)]
		It = I10[np.ix_(F, [i])]

		data_W10.extend(wq_solve_equation_system(Bt, It))
		indi.extend(i*np.ones(len(P), dtype= int))
		indj.extend(P)

	W10 = sp.csr_matrix((data_W10, (indi, indj)), shape=(nb_ctrlpts, nb_qp_wq))

	# -------------------
	# Computation of W11
	# -------------------
	Ishape = shape_B_p1 @ shape_B1_p0.T
	indi = []; indj = []
	for i in range(nb_ctrlpts):
		F = np.nonzero(Ishape[:, i])[0]
		P = np.nonzero(shape_B1_p0[i, :])[1]

		Bt = B0wq_p1[np.ix_(F,P)]
		It = I11[np.ix_(F, [i])]
			
		data_W11.extend(wq_solve_equation_system(Bt, It))
		indi.extend(i*np.ones(len(P), dtype= int))
		indj.extend(P)
	W11 = sp.csr_matrix((data_W11, (indi, indj)), shape=(nb_ctrlpts, nb_qp_wq))
	
	return W00, W01, W10, W11 

def wq_find_basis_weights_opt(degree, knotvector, r): 
	" Return basis and weights in WQ approach computed in an efficient way "
	
	nbel = len(np.unique(knotvector)) - 1
	xwq = wq_find_positions(degree, knotvector, r)

	if nbel <= degree + 3 : 
		B0, B1 = eval_basis_python(degree, knotvector, xwq)
		W00, W01, W10,  W11 = wq_find_weights(degree, knotvector, r)
	else : 
		# Create model
		nbel_model = degree + 3
		knotvector_model = create_knotvector(degree, nbel_model)
		xwq_model = wq_find_positions(degree, knotvector_model, r)
		B0_model, B1_model = eval_basis_python(degree, knotvector_model, xwq_model)
		W00_model, W01_model, W10_model,  W11_model = wq_find_weights(degree, knotvector_model, r)
		shape_B0_model, shape_B1_model = wq_get_shape_B(degree, nbel, r)

		# Modify the results
		W00_model *= nbel_model/nbel
		W01_model *= nbel_model/nbel
		B1_model *= nbel/nbel_model

		#------------------------------
		# Algorithm to get the results
		#------------------------------
		nbx = 2*(nbel + degree + r) - 5
		nbfunct = degree + nbel

		data_B0, data_B1, data_W00, data_W01, data_W10, data_W11 = [], [], [], [], [], []
		function_B0, support_B0, function_B1, support_B1 = [], [], [], []
		shape_B0, shape_B1 = wq_get_shape_B(degree, nbel, r)

		for i in range(degree+1):
			# For B0 and normal
			indj = shape_B0.getrow(i).nonzero()[1]
			function_B0.extend(i * np.ones(len(indj), dtype=int))
			support_B0.extend(indj)
			data_B0.extend(B0_model[np.ix_([i],indj)].toarray()[0].tolist())
			data_W00.extend(W00_model[np.ix_([i],indj)].toarray()[0].tolist())
			data_W01.extend(W01_model[np.ix_([i],indj)].toarray()[0].tolist())

			# For B0 and flip
			indj_flip = shape_B0.getrow(nbfunct-i-1).nonzero()[1]
			function_B0.extend((nbfunct-i-1) * np.ones(len(indj_flip), dtype=int))
			support_B0.extend(indj_flip)
			data_B0.extend(np.flip(B0_model[np.ix_([i],indj)]).toarray()[0].tolist())         
			data_W00.extend(np.flip(W00_model[np.ix_([i],indj)]).toarray()[0].tolist())   
			data_W01.extend(np.flip(W01_model[np.ix_([i],indj)]).toarray()[0].tolist())   
			
			# For B1 and normal
			indj = shape_B1.getrow(i).nonzero()[1]
			function_B1.extend(i * np.ones(len(indj), dtype=int))
			support_B1.extend(indj)
			data_B1.extend(B1_model[np.ix_([i],indj)].toarray()[0].tolist())
			data_W10.extend(W10_model[np.ix_([i],indj)].toarray()[0].tolist())
			data_W11.extend(W11_model[np.ix_([i],indj)].toarray()[0].tolist())

			# For B1 and flip
			indj_flip = shape_B1.getrow(nbfunct-i-1).nonzero()[1]
			function_B1.extend((nbfunct-i-1) * np.ones(len(indj_flip), dtype=int))
			support_B1.extend(indj_flip)
			data_B1.extend(np.flip(-B1_model[np.ix_([i],indj)]).toarray()[0].tolist())         
			data_W10.extend(np.flip(-W10_model[np.ix_([i],indj)]).toarray()[0].tolist())   
			data_W11.extend(np.flip(-W11_model[np.ix_([i],indj)]).toarray()[0].tolist())   

		for i in range(degree+1,nbfunct-degree-1): 
			# For B0
			indj_model = shape_B0_model.getrow(degree+1).nonzero()[1]
			indj = shape_B0.getrow(i).nonzero()[1]
			function_B0.extend(i * np.ones(len(indj), dtype=int))
			support_B0.extend(indj)
			data_B0.extend(B0_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
			data_W00.extend(W00_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
			data_W01.extend(W01_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())

			# For B1
			indj_model = shape_B1_model.getrow(degree+1).nonzero()[1]
			indj = shape_B1.getrow(i).nonzero()[1]
			function_B1.extend(i * np.ones(len(indj), dtype=int))
			support_B1.extend(indj)
			data_B1.extend(B1_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
			data_W10.extend(W10_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
			data_W11.extend(W11_model[np.ix_([degree+1],indj_model)].toarray()[0].tolist())
	
		B0 = sp.csr_matrix((data_B0, (function_B0, support_B0)), shape=(nbfunct, nbx))
		W00 = sp.csr_matrix((data_W00, (function_B0, support_B0)), shape=(nbfunct, nbx))    
		W01 = sp.csr_matrix((data_W01, (function_B0, support_B0)), shape=(nbfunct, nbx))    

		B1 = sp.csr_matrix((data_B1, (function_B1, support_B1)), shape=(nbfunct, nbx))    
		W10 = sp.csr_matrix((data_W10, (function_B1, support_B1)), shape=(nbfunct, nbx))    
		W11 = sp.csr_matrix((data_W11, (function_B1, support_B1)), shape=(nbfunct, nbx))  

	return xwq, B0, B1, W00, W01, W10, W11

def wq_find_basis_weights_fortran(degree, knotvector): 
	" Computes basis and weights in WQ approach using Fortran "

	# Set guessed size of data 
	nnz_B, nb_qp = basis_weights.wq_get_size_data(degree, knotvector)

	# Get basis and weights from fortran
	qp_position, basis, weights, indi, indj, nnz_I = basis_weights.wq_get_data_csr(
													degree, knotvector, nnz_B, nb_qp)

	return nnz_I, qp_position, basis, weights, indi, indj

# =========================
# MF FUNCTIONS
# =========================
def create_table_properties(function, uref=None, prop=None):
	"Create a table of scalar properties from given function "

	# Set default given x
	if uref is None: uref = np.linspace(-1., 1., 201)

	# Compute y
	if prop is None: y = function(uref)
	else:  y = function(uref, prop)

	# Create table 
	table = np.zeros((len(uref), 2))
	table[:, 0] = uref; table[:, 1] = y

	return table

def solver_scipy(A, b, nbIterations=100, epsilon=1e-10, PreCond='ilu', isCG=True):
	""" Solve system using iterative method : conjugate gradient or bi-conjugate gradient. 
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

def eigen_decomposition(indi, indj, data, robin_condition=[0, 0], coefs=None):
	""" 
		Eigen decomposition generalized KU = MUD
		K: stiffness matrix, K = int B1 B1 dx = W11 * B1
		M: mass matrix, M = int B0 B0 dx = W00 * B0
		U: eigenvectors matrix
		D: diagonal of eigenvalues
	"""

	nb_qp = np.max(indj)
	[B0, B1, W0, W1] = data

	# Create mcoef and kcoef is are none type
	if coefs is None: mcoefs = np.ones(nb_qp); kcoefs = np.ones(nb_qp)
	else: [mcoefs, kcoefs] = coefs

	eigenvalues, eigenvectors = solver.eigen_decomposition_py(indi, indj, B0, W0, 
										B1, W1, mcoefs, kcoefs, robin_condition)

	return eigenvalues, eigenvectors

def tensor_decomposition_3d(n_list, CC_matrix, dim=3):
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

def compute_eig_diag(eig_u, eig_v, eig_w, coefs=[1.0, 1.0, 1.0]):
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

def fast_diagonalization(U, V, W, D, array_in, fdtype='steady'):
	" Compute fast diagonalization using Fortran"
	if fdtype == 'interp':
		array_out = solver.fd_interpolation_3d_py(U, V, W, array_in)
	elif fdtype == 'steady':
		array_out = solver.fd_steady_heat_3d_py(U, V, W, D, array_in)
	elif fdtype == 'elastic':
		array_out = elastoplasticity.fd_elasticity_3d_py(U, V, W, D, array_in)
	
	return array_out

def compute_mean_3d(ncu, ncv, ncw, C):

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