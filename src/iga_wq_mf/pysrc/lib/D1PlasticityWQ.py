"""
.. This module contains functions to treat elastoplasticity problems in 1D. 
.. It is adapted to IGA problems (using Gauss quadrature)
.. Joaquin Cornejo
"""

import numpy as np

def cpp_combined_hardening_1D_s1(properties, strain, pls, a, b):
	""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
		It uses combined isotropic/kinematic hardening theory.  
	"""
	E, H, beta, sigma_Y0 = properties

	# Elastic predictor
	sigma_trial = E*(strain - pls)
	eta_trial = sigma_trial - b

	# Check yield status
	f_trial = abs(eta_trial) - (sigma_Y0 + beta*H*a)

	if f_trial <= sigma_Y0*1e-6: # Elastic
		stress = sigma_trial
		pls_new = pls
		a_new = a
		b_new = b
		
	else: # Plastic
		N = np.sign(eta_trial)
		dgamma = f_trial/(E + H)
		stress = sigma_trial - dgamma*E*N
		pls_new = pls + dgamma*N
		a_new = a + dgamma
		b_new = b + (1-beta)*H*dgamma*N

	return [stress, pls_new, a_new, b_new]

def cpp_combined_hardening_1D_s2(properties, strain, pls, a, b):
	""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
		It uses combined isotropic/kinematic hardening theory.  
	"""
	E, H, beta, sigma_Y0 = properties

	# Elastic predictor
	sigma_trial = E*(strain - pls)
	eta_trial = sigma_trial - b

	# Check yield status
	f_trial = abs(eta_trial) - (sigma_Y0 + beta*H*a)

	if f_trial <= sigma_Y0*1e-6: # Elastic
		Cep = E
	else: # Plastic
		Cep = E*H/(E + H)

	return Cep

def compute_static_Fint_1D(DW, sigma):
	""" Computes internal force Fint. 
		Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
	"""
	Fint = DW[0] @ sigma.T
	return Fint

def compute_volForce_1D(DW, b): 
	" Computes volumetric force in 1D. Usualy b = rho*g*L (Weight) "
	F = DW[0] @ b.T
	return F

def compute_tangent_static_matrix_1D(JJ, DB, DW, Cep):
	""" Computes stiffness matrix in elastoplasticity
		S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
	"""
	Scoefs = Cep*1.0/JJ
	S = DW[-1] @ np.diag(Scoefs) @ DB[1].T 
	return S

def interpolate_strain_1D(JJ, DW, disp):
	" Computes strain field from a given displacement field "
	eps = DW[-1].T @ disp / JJ
	return eps

def solve_plasticity_1D(properties, DB=None, DW=None, Fext=None, dof=None, tol=1e-8, nbIterNL=10):
	" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

	JJ, nb_qp = properties[0], properties[-1]
	ep_n0, a_n0, b_n0 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
	ep_n1, a_n1, b_n1 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
	Cep, sigma = np.zeros(nb_qp), np.zeros(nb_qp)
	disp = np.zeros(np.shape(Fext))

	strain = np.zeros((nb_qp, np.shape(Fext)[1]))
	stress = np.zeros((nb_qp, np.shape(Fext)[1]))
	energy = np.zeros(np.shape(Fext)[1]-1)
	internal =  np.zeros(np.shape(Fext)[1]-1)

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
			Fint = compute_static_Fint_1D(DW, sigma)
			dF 	 = Fstep[dof] - Fint[dof]
			relerror = np.sqrt(np.dot(dF, dF))
			print('Rhapson with error %.5e' %relerror)
			if relerror <= tol: break

			# Compute stiffness
			Sdof = compute_tangent_static_matrix_1D(JJ, DB, DW, Cep)[np.ix_(dof, dof)]
			ddisp += np.linalg.solve(Sdof, dF)

		# # Verify energy conservation
		# S = compute_tangent_static_matrix_1D(JJ, DB, DW, Cep)
		# Fd = np.dot(Fstep, d_n1)
		# dKd = np.dot(d_n1, np.dot(S, d_n1))
		# energy[i-1] = 0.5*dKd - Fd
		# internal[i-1] = np.dot(W, sigma*(eps-ep_n1)) - dKd

		# Update values in output
		disp[:, i]   = d_n1
		strain[:, i] = eps
		stress[:, i] = sigma

		ep_n0, a_n0, b_n0 = np.copy(ep_n1), np.copy(a_n1), np.copy(b_n1)

	return disp, strain, stress, energy, internal

def interpolate_controlPoints_1D(DB, DW, u_ref):
	" Interpolate control point field (from parametric to physical space) "
	masse = DW[0] @ DB[0].T
	force = DW[0] @ u_ref
	u_ctrlpts = np.linalg.solve(masse, force)
	return u_ctrlpts