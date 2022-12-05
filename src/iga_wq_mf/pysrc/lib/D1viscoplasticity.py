"""
.. This module contains functions to treat elastoplasticity problems in 1D. 
.. It is adapted to IGA problems (using Gauss quadrature)
.. Joaquin Cornejo
"""

import numpy as np

def cpp_combined_hardening_1D_s1(properties, eps, ep_n0, a_n0, b_n0):
	""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
		It uses combined isotropic/kinematic hardening theory.  
	"""
	E, H, beta, sigma_Y0 = properties

	# Elastic predictor
	sigma_trial = E*(eps - ep_n0)
	eta_trial = sigma_trial - b_n0

	# Check yield status
	f_trial = abs(eta_trial) - (sigma_Y0 + beta*H*a_n0)

	if f_trial <= sigma_Y0*1e-6: # Elastic
		sigma = sigma_trial
		ep_n1 = ep_n0
		a_n1 = a_n0
		b_n1 = b_n0
		
	else: # Plastic
		N = np.sign(eta_trial)
		dgamma = f_trial/(E + H)
		sigma = sigma_trial - dgamma*E*N
		ep_n1 = ep_n0 + dgamma*N
		a_n1 = a_n0 + dgamma
		b_n1 = b_n0 + (1-beta)*H*dgamma*N

	return [sigma, ep_n1, a_n1, b_n1]

def cpp_combined_hardening_1D_s2(properties, eps, ep_n0, a_n0, b_n0):
	""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
		It uses combined isotropic/kinematic hardening theory.  
	"""
	E, H, beta, sigma_Y0 = properties

	# Elastic predictor
	sigma_trial = E*(eps - ep_n0)
	eta_trial = sigma_trial - b_n0

	# Check yield status
	f_trial = abs(eta_trial) - (sigma_Y0 + beta*H*a_n0)

	if f_trial <= sigma_Y0*1e-6: # Elastic
		Cep = E
	else: # Plastic
		Cep = E*H/(E + H)

	return Cep

def compute_static_Fint_1D(DB, W, sigma):
	""" Computes internal force Fint. 
		Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
	"""
	Fint = DB[1] @ np.diag(W) @ sigma.T
	return Fint

def compute_volForce_1D(DB, W, b): 
	" Computes volumetric force in 1D. Usualy b = rho*g*L (Weight) "
	F = DB[0] @ np.diag(W) @ b.T
	return F

def compute_tangent_static_matrix_1D(JJ, DB, W, Cep):
	""" Computes stiffness matrix in elastoplasticity
		S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
	"""
	Scoefs = W*Cep*1.0/JJ
	S = DB[1] @ np.diag(Scoefs) @ DB[1].T 
	return S

def interpolate_strain_1D(JJ, DB, disp):
	" Computes strain field from a given displacement field "
	eps = DB[1].T @ disp / JJ
	return eps

def solve_plasticity_1D(properties, DB=None, W=None, Fext=None, dof=None, tol=1e-8, nbIterNL=10, update=1):
	" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

	JJ, nb_qp = properties[0], properties[-1]
	ep_n0, a_n0, b_n0 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
	ep_n1, a_n1, b_n1 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
	Cep, sigma = np.zeros(nb_qp), np.zeros(nb_qp)
	disp = np.zeros(np.shape(Fext))

	strain = np.zeros((nb_qp, np.shape(Fext)[1]))
	stress = np.zeros((nb_qp, np.shape(Fext)[1]))

	for i in range(1, np.shape(Fext)[1]):

		d_n1 = np.copy(disp[:, i-1]); ddisp = np.zeros(len(dof)) 
		Fstep = np.copy(Fext[:, i])

		print('Step %d' %i)
		for j in range(nbIterNL): # Newton-Raphson

			# Compute strain as function of displacement
			d_n1[dof] = disp[dof, i-1] + ddisp
			eps = interpolate_strain_1D(JJ, DB, d_n1)

			# Find closest point projection 
			for k in range(nb_qp):
				result1 = cpp_combined_hardening_1D_s1(properties[1:-1], 
							eps[k], ep_n0[k], a_n0[k], b_n0[k])

				result2 = cpp_combined_hardening_1D_s2(properties[1:-1], 
							eps[k], ep_n0[k], a_n0[k], b_n0[k])

				sigma[k], ep_n1[k], a_n1[k], b_n1[k] = result1
				Cep[k] = result2

			# Compute Fint
			Fint = compute_static_Fint_1D(DB, W, sigma)
			dF = Fstep[dof] - Fint[dof]
			relerror = np.sqrt(np.dot(dF, dF))
			print('Rhapson with error %.5e' %relerror)
			if relerror <= tol: break

			# Compute stiffness
			S = compute_tangent_static_matrix_1D(JJ, DB, W, Cep)[np.ix_(dof, dof)]
			ddisp += np.linalg.solve(S, dF)

		# Update values in output
		disp[:, i]   = d_n1
		strain[:, i] = eps
		stress[:, i] = sigma

		ep_n0, a_n0, b_n0 = np.copy(ep_n1), np.copy(a_n1), np.copy(b_n1)

	return disp, strain, stress

def interpolate_controlPoints_1D(DB, W, u_ref):
	" Interpolate control point field (from parametric to physical space) "
	masse = DB[0] @ np.diag(W) @ DB[0].T
	force = DB[0] @ np.diag(W) @ u_ref
	u_ctrlpts = np.linalg.solve(masse, force)
	return u_ctrlpts