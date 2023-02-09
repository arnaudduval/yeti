"""
.. This module contains functions to treat elastoplasticity problems in 1D. 
.. It is adapted to IGA problems (using Gauss quadrature)
.. Joaquin Cornejo
"""

import numpy as np
from lib.__init__ import *
from lib.base_functions import eval_basis_python
from mpl_toolkits.axes_grid1 import make_axes_locatable

class MechaBehavior():
	
	def __init__(self, name:str, kwargs:dict):

		# General init
		self._young = kwargs.get('young', 1.0)
		self.Hfun, self.Hderfun = None, None
		self.Kfun, self.Kderfun = None, None

		# Specific init
		if name.lower() == 'linear':
			self._init_linear(kwargs)

		if name.lower() == 'swift':
			self._init_swift(kwargs)

		if name.lower() == 'voce':
			self._init_voce(kwargs)

		funlist = [self.Hfun, self.Hderfun, self.Kfun, self.Kderfun]
		if any([fun is None for fun in funlist]): raise Warning('Something went wrong')
		return
	
	def _init_linear(self, kwargs:dict):
		sigma_y0  = kwargs.get('sigma_Y0', 1.0)
		theta	  = kwargs.get('theta', 1.0)
		H_bar     = kwargs.get('H_bar', 1.0)
		self.Kfun = lambda a: sigma_y0 + theta*H_bar*a 
		self.Hfun = lambda a: (1-theta)*H_bar*a
		self.Kderfun = lambda a: theta*H_bar
		self.Hderfun = lambda a: (1-theta)*H_bar
		return
	
	def _init_swift(self, kwargs:dict):
		sigma_y0 = kwargs.get('sigma_Y0', 1.0)
		K = kwargs.get('K', 1.0)
		n = kwargs.get('exp', 1.0)
		self.Kfun = lambda a: sigma_y0 + self._young*(a/K)**n
		self.Hfun = lambda a: 0.0
		self.Kderfun = lambda a: (self._young/K)*n*(a/K)**(n-1) if a!=0 else 1e10*self._young
		self.Hderfun = lambda a: 0.0
		return
	
	def _init_voce(self, kwargs:dict):
		sigma_y0  = kwargs.get('sigma_Y0', 1.0)
		theta	  = kwargs.get('theta', 1.0)
		H_bar     = kwargs.get('H_bar', 1.0)
		K_inf     = kwargs.get('K_inf', 1.0)
		delta     = kwargs.get('delta', 1.0)
		self.Kfun = lambda a: sigma_y0 + theta*H_bar*a + K_inf*(1.0 - np.exp(-delta*a))
		self.Hfun = lambda a: (1-theta)*H_bar*a
		self.Kderfun = lambda a: theta*H_bar + K_inf*delta*np.exp(-delta*a) 
		self.Hderfun = lambda a: (1-theta)*H_bar
		return
	
	def compute_deltaGamma(self, eta_trial, a_n0, nbIter=20, threshold=1e-8):
		dgamma = 0.0
		a_n1   = a_n0 
		for i in range(nbIter):
			dH = self.Hfun(a_n1) - self.Hfun(a_n0) 
			G  = -self.Kfun(a_n1) + np.abs(eta_trial) - (self._young*dgamma + dH)
			if G <=threshold: break
			dG = - (self._young + self.Hderfun(a_n1) + self.Kderfun(a_n1))
			dgamma = dgamma - G/dG
			a_n1   = a_n0 + dgamma 
		return dgamma
	
	def return_mapping(self, strain, pls, a, b, threshold=1e-8):
		""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		# Elastic predictor
		sigma_trial = self._young*(strain - pls)
		eta_trial   = sigma_trial - b

		# Check yield status
		f_trial = np.abs(eta_trial) - self.Kfun(a)

		if f_trial <= threshold: # Elastic
			stress = sigma_trial
			pls_new = pls
			a_new = a
			b_new = b
			Cep = self._young

		else: # Plastic
			N = np.sign(eta_trial)
			dgamma = self.compute_deltaGamma(eta_trial, a)
			stress = sigma_trial - dgamma*self._young*N
			pls_new = pls + dgamma*N
			a_new = a + dgamma
			b_new = b + (self.Hfun(a_new)-self.Hfun(a))*N
			somme = self.Kderfun(a_new) + self.Hderfun(a_new)
			Cep = self._young*somme/(self._young + somme)

		return [stress, pls_new, a_new, b_new, Cep]

def interpolate_strain_1D(JJ, DB, disp):
	" Computes strain field from a given displacement field "
	eps = DB[1].T @ disp / JJ
	return eps

# ----

def compute_IGA_Fint_1D(DB, W, sigma):
	""" Computes internal force Fint. 
		Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
	"""
	Fint = DB[1] @ np.diag(W) @ sigma.T
	return Fint

def compute_WQ_Fint_1D(DW, sigma):
	""" Computes internal force Fint. 
		Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
	"""
	Fint = DW[2] @ sigma.T
	# Fint = DW[-1] @ sigma.T
	return Fint

# ----

def compute_IGA_Fvol_1D(DB, W, b): 
	" Computes volumetric force in 1D. Usualy b = rho*g*L (Weight) "
	F = DB[0] @ np.diag(W) @ b.T
	return F

def compute_WQ_Fvol_1D(DW, b): 
	" Computes volumetric force in 1D. Usualy b = rho*g*L (Weight) "
	F = DW[0] @ b.T
	return F

# ----

def compute_IGA_tangent_matrix_1D(JJ, DB, W, Cep):
	""" Computes stiffness matrix in elastoplasticity
		S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
	"""
	Scoefs = W*Cep*1.0/JJ
	S = DB[1] @ np.diag(Scoefs) @ DB[1].T 
	return S

def compute_WQ_tangent_matrix_1D(JJ, DB, DW, Cep):
	""" Computes stiffness matrix in elastoplasticity
		S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
	"""
	Scoefs = Cep*1.0/JJ
	S = DW[-1] @ np.diag(Scoefs) @ DB[1].T 
	return S

# ----

def solve_IGA_plasticity_1D(properties, matlaw:MechaBehavior, DB=None, W=None, Fext=None, dof=None, threshold=1e-8, nbIterNL=10):
	" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

	JJ, nb_qp = properties[0], properties[-1]
	ep_n0, a_n0, b_n0 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
	ep_n1, a_n1, b_n1 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
	Cep, sigma = np.zeros(nb_qp), np.zeros(nb_qp)
	disp = np.zeros(np.shape(Fext))

	strain  = np.zeros((nb_qp, np.shape(Fext)[1]))
	plastic = np.zeros((nb_qp, np.shape(Fext)[1]))
	stress  = np.zeros((nb_qp, np.shape(Fext)[1]))
	moduleE = np.zeros((nb_qp, np.shape(Fext)[1]))

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
				result1 = matlaw.return_mapping(eps[k], ep_n0[k], a_n0[k], b_n0[k])
				sigma[k], ep_n1[k], a_n1[k], b_n1[k], Cep[k] = result1

			# Compute Fint
			Fint = compute_IGA_Fint_1D(DB, W, sigma)
			dF 	 = Fstep[dof] - Fint[dof]
			relerror = np.sqrt(np.dot(dF, dF))
			print('Rhapson with error %.5e' %relerror)
			if relerror <= threshold: break

			# Compute stiffness
			Sdof = compute_IGA_tangent_matrix_1D(JJ, DB, W, Cep)[np.ix_(dof, dof)]
			ddisp += np.linalg.solve(Sdof, dF)
			
		# Update values in output
		disp[:, i]   = d_n1
		strain[:, i] = eps
		stress[:, i] = sigma
		plastic[:, i] = ep_n1
		moduleE[:, i] = Cep

		ep_n0, a_n0, b_n0 = np.copy(ep_n1), np.copy(a_n1), np.copy(b_n1)

	return disp, strain, stress, plastic, moduleE

def solve_WQ_plasticity_1D(properties, matlaw:MechaBehavior, DB=None, DW=None, Fext=None, dof=None, threshold=1e-8, nbIterNL=10):
	" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

	JJ, nb_qp = properties[0], properties[-1]
	ep_n0, a_n0, b_n0 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
	ep_n1, a_n1, b_n1 = np.zeros(nb_qp), np.zeros(nb_qp), np.zeros(nb_qp)
	Cep, sigma = np.zeros(nb_qp), np.zeros(nb_qp)
	disp = np.zeros(np.shape(Fext))

	strain  = np.zeros((nb_qp, np.shape(Fext)[1]))
	plastic = np.zeros((nb_qp, np.shape(Fext)[1]))
	stress  = np.zeros((nb_qp, np.shape(Fext)[1]))
	moduleE = np.zeros((nb_qp, np.shape(Fext)[1]))

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
				result1 = matlaw.return_mapping(eps[k], ep_n0[k], a_n0[k], b_n0[k])
				sigma[k], ep_n1[k], a_n1[k], b_n1[k], Cep[k] = result1

			# Compute Fint
			Fint = compute_WQ_Fint_1D(DW, sigma)
			dF 	 = Fstep[dof] - Fint[dof]
			relerror = np.sqrt(np.dot(dF, dF))
			print('Rhapson with error %.5e' %relerror)
			if relerror <= threshold: break

			# Compute stiffness
			Sdof = compute_WQ_tangent_matrix_1D(JJ, DB, DW, Cep)[np.ix_(dof, dof)]
			ddisp += np.linalg.solve(Sdof, dF)

		# Update values in output
		disp[:, i]   = d_n1
		strain[:, i] = eps
		stress[:, i] = sigma
		plastic[:, i] = ep_n1
		moduleE[:, i] = Cep

		ep_n0, a_n0, b_n0 = np.copy(ep_n1), np.copy(a_n1), np.copy(b_n1)

	return disp, strain, stress, plastic, moduleE

# ----

def interpolate_IGA_CP_1D(DB, W, u_ref):
	" Interpolate control point field (from parametric to physical space) "
	masse = DB[0] @ np.diag(W) @ DB[0].T
	force = DB[0] @ np.diag(W) @ u_ref
	u_ctrlpts = np.linalg.solve(masse, force)
	return u_ctrlpts

def interpolate_WQ_CP_1D(DB, DW, u_ref):
	" Interpolate control point field (from parametric to physical space) "
	masse = DW[0] @ DB[0].T
	force = DW[0] @ u_ref
	u_ctrlpts = np.linalg.solve(masse.toarray(), force)
	return u_ctrlpts

# =============

def plot_results(degree, knotvector, JJ, disp_cp, plastic_cp, stress_cp, folder=None, method='IGA', extension='.png'):
	knots  = np.linspace(0, 1, 101)
	DB     = eval_basis_python(degree, knotvector, knots)
	displacement   = DB[0].T @ disp_cp
	strain_interp  = DB[1].T @ disp_cp
	plastic_interp = DB[0].T @ plastic_cp
	stress_interp  = DB[0].T @ stress_cp

	# Plot fields
	N = np.shape(disp_cp)[1]
	XX, STEPS = np.meshgrid(knots*JJ, np.arange(N))
	names = ['Displacement field', 'Plastic strain field', 'Stress field']
	units = ['m', '\%', 'MPa']
	fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
	for ax, variable, name, unit in zip([ax1, ax2, ax3], [displacement, plastic_interp, stress_interp], names, units):
		im = ax.pcolormesh(XX, STEPS, variable.T, cmap='PuBu_r', shading='linear')
		ax.set_title(name)
		ax.set_ylabel('Step')
		ax.set_xlabel('Position (m)')
		ax.grid(False)
		divider = make_axes_locatable(ax)
		cax = divider.append_axes('right', size='5%', pad=0.05)
		cbar = fig.colorbar(im, cax=cax)
		cbar.ax.set_title(unit)

	fig.tight_layout()
	fig.savefig(folder + 'ElastoPlasticity' + method + extension)

	# Plot stress-strain of single point
	fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
	for ax, pos in zip([ax1, ax2, ax3], [25, 50, 75]):
		ax.plot(strain_interp[pos, :]*100, stress_interp[pos, :])
		ax.set_ylabel('Stress (MPa)')
		ax.set_xlabel('Strain (\%)')
		ax.set_ylim(bottom=0.0, top=1500)
		ax.set_xlim(left=0.0, right=strain_interp.max()*100)

	fig.tight_layout()
	fig.savefig(folder + 'TractionCurve' + method + extension)
	return
