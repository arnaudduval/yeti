"""
.. This module contains functions to treat elastoplasticity problems in 3D. 
.. Joaquin Cornejo
.. 
.. Remarks :: Voigt notation is used. That is 
.. Strain e = [e11, e22, e33, 2e12, 2e13, 2e23]^T
.. Stress s = [s11, s22, s33, s12, s13, s23]^T
.. 1 = [1, 1, 1, 0, 0, 0]^T
.. I = diag(1, 1, 1, 0.5, 0.5, 0.5)
.. And then, some products are:
.. Double contraction st-st s:e = dot(s, e), this only works with stress and strains
.. In Voigt notation if one wants to compute e:e or s:s some scale factor need to be considered.
.. Double contraction C:e = matmul(C, e), this is true only with strains
.. Otherwise some scale factor need to be considered.
"""

import numpy as np

def clean_dirichlet_3d(A, dod):
	""" Set to 0 (Dirichlet condition) the values of an array using the indices in each dimension
		A is actually a vector arranged following each dimension [Au, Av, Aw]
	"""
	for i in range(3): A[i, dod[i]] = 0.0
	return

def block_dot_product(d, A, B):
	""" Computes dot product of A and B. 
		Both are actually vectors arranged following each dimension
		A = [Au, Av, Aw] and B = [Bu, Bv, Bw]. Then A.B = Au.Bu + Av.Bv + Aw.Bw
	"""
	result = 0.0
	for i in range(d): result += A[i, :] @ B[i, :]
	return result

def compute_trace(d, tensor):
	" Computes trace of a second-order stress-like tensor "

	trace = 0.0
	for i in range(d):
		trace += tensor[i]
		
	return trace

def compute_deviatoric(d, tensor):
	" Computes deviatoric of a second-order stress-like tensor "

	ddl = int(d*(d+1)/2)
	one = np.zeros(ddl)
	for i in range(d): one[i] = 1.0

	trace = compute_trace(d, tensor)		
	dev = tensor - 1.0/3.0*trace*one

	return dev

def compute_stress_norm(d, tensor): 
	" Returns frobenius norm of a second-order stress-like tensor "

	ddl = int(d*(d+1)/2)
	norm = 0.0

	for i in range(d): norm += tensor[i]**2
	for i in range(d, ddl): norm += 2.0*tensor[i]**2
	norm = np.sqrt(norm)

	return norm

def create_fourth_order_identity(d=3):
	" Creates a fourth-order identity (Voigt notation) "
	
	if d == 2:
		I = np.diag(np.array([1.0, 1.0, 0.5]))
	elif d == 3:    
		I = np.diag(np.array([1.0, 1.0, 1.0, 0.5, 0.5, 0.5]))
	else: raise Warning('Only 2d or 3d')

	return I

def create_second_order_identity(d=3):
	" Creates a fourth-order identity (Voigt notation) "
	
	if d == 2:
		one = np.array([1.0, 1.0, 0.0])
	elif d == 3:    
		one = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
	else: raise Warning('Only 2d or 3d')

	return one

def create_one_kron_one(d=3): 
	" Creates a one kron one tensor (Voigt notation) "

	ddl = int(d*(d+1)/2)
	onekronone = np.zeros((ddl, ddl))

	for i in range(d):
		for j in range(d):
			onekronone[i,j] = 1.0

	return onekronone

def stkronst(A, B):
	""" Computes kron product of tensors A and B
		A and B are second order stress-like tensor in Voigt notation
	"""

	At = np.atleast_2d(A)
	Bt = np.atleast_2d(B)
	result = np.kron(At.T, Bt)

	return result

def create_incidence_matrix(d=3):
	""" Creates M matrix. E is the passage matrix from derivative to actual symetric values. 
		If we multiply a vector of values u_(i, j) with E matrix, one obtains the vector: 
		us_ij = 0.5*(u_(i,j) + u_(j,i))  
	"""

	ddl = int(d*(d+1)/2)
	EE = np.zeros((ddl, d, d))

	if d==3: 
		EE[0, 0, 0] = 1.0; EE[4, 2, 0] = 1.0; EE[5, 1, 0] = 1.0
		EE[1, 1, 1] = 1.0; EE[3, 2, 1] = 1.0; EE[5, 0, 1] = 1.0
		EE[2, 2, 2] = 1.0; EE[3, 1, 2] = 1.0; EE[4, 0, 2] = 1.0
	elif d == 2: 
		EE[0, 0, 0] = 1.0; EE[2, 1, 0] = 1.0
		EE[1, 1, 1] = 1.0; EE[2, 0, 1] = 1.0

	return EE

def cpp_combined_hardening(inputs, eps, ep_n, a_n, b_n, d=3):
	" Returns closest point proyection (cpp) in combined hardening plasticity criteria"

	CC, sigma_Y, bulk, mu, beta, H, I, Idev, one = inputs

	trace_eps = compute_trace(d, eps)
	e_n1 = compute_deviatoric(d, eps)
	s_trial = 2*mu*(e_n1 - ep_n) @ I
	eta_trial = s_trial - b_n

	# Check status
	norm_eta = compute_stress_norm(d, eta_trial)
	f = norm_eta - np.sqrt(2.0/3.0)*(sigma_Y + beta*H*a_n)

	sigma = s_trial + bulk*trace_eps*one
	Cep = CC
	if f <= 0:
		# Copy old variables
		ep_n1 = ep_n
		a_n1 = a_n
		b_n1 = b_n
		
	else:

		# Consistency parameter
		dgamma = f/(2.0*mu+2.0/3.0*H)
		n_n1 = eta_trial/norm_eta
		a_n1 = a_n + np.sqrt(2.0/3.0)*dgamma

		# Update stress
		b_n1 = b_n + 2.0/3.0*(1-beta)*H*dgamma*n_n1
		ep_n1 = ep_n + dgamma*n_n1
		sigma -= 2*mu*dgamma*n_n1

		# Compute consistent tangent matrix
		c1 = 1.0 - 2.0*mu*dgamma/norm_eta
		c2 = 1.0/(1.0 + H/(3.0*mu)) - (1.0 - c1)
		NNT = stkronst(n_n1, n_n1)
		Cep -= 2*mu*((1 - c1)*Idev + c2*NNT)

	return sigma, ep_n1, a_n1, b_n1, Cep

def compute_plasticity_coef(sigma, Cep, invJ, detJ, d=3):
	" Computes the coefficients to use in internal force vector and stiffness matrix"

	EE = create_incidence_matrix(d)
	coef_Fint = np.zeros((d*d, len(detJ)))
	coef_Stiff = np.zeros((d*d, d*d, len(detJ)))

	for k, det in enumerate(detJ):

		for i in range(d):
			for j in range(d):
				Dij = invJ[:,:,k] @ EE[:,:,i].T @ Cep[:,:,k] @ EE[:,:,j] @ invJ[:,:,k].T
				coef_Stiff[i*d:(i+1)*d, j*d:(j+1)*d, k] = Dij*det

			Si = invJ[:,:,k] @ EE[:,:,i].T @ sigma[:,k]
			coef_Fint[i*d:(i+1)*d, k] = Si*det

	return coef_Fint, coef_Stiff

def compute_stress_vonmises(d, tensor):
	" Computes equivalent stress with Von Mises formula of a second-order stress-like tensor "

	dev = compute_deviatoric(d, tensor)
	vm = compute_stress_norm(d, dev)
	vm = np.sqrt(3.0/2.0)*vm

	return vm