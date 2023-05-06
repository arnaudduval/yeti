"""
.. This module contains functions to treat transient heat problems in 1D. 
.. It also contains functions to treat elastoplasticity problems in 1D. 
.. It is adapted to IGA problems (using Gauss quadrature)
.. Joaquin Cornejo
"""

import numpy as np
from lib.__init__ import *
from lib.lib_base import array2csr_matrix, evalDersBasisFortran
from lib.lib_quadrules import *
from lib.lib_material import plasticLaw

class model1D:
	def __init__(self, **kwargs):
		self._sampleSize = kwargs.get('sample_size', 101)
		self._kwargs = kwargs
		self.set_geometry()
		self.set_igaparameterization()
		self.set_quadratureRule()
		return
	
	def set_geometry(self):
		self._J    = self._kwargs.get('length', 1.0)
		self._detJ = self._J
		self._invJ = 1.0/self._J
		return
	
	def set_igaparameterization(self):
		kwargs           = self._kwargs
		self._degree     = kwargs.get('degree', 2)
		self._knotvector = kwargs.get('knotvector', np.array([0, 0, 0, 0.5, 1, 1, 1]))
		self._nbctrlpts  = len(self._knotvector) - self._degree - 1
		return

	def set_quadratureRule(self):
		kwargs        = self._kwargs
		quadRuleName  = kwargs.get('quadrule', 'wq').lower()
		quadRule      = None
		
		if quadRuleName == 'iga':
			quadRule = GaussQuadrature(self._degree, self._knotvector, kwargs=kwargs)
		elif quadRuleName == 'wq':
			quadRule = WeightedQuadrature(self._degree, self._knotvector, kwargs=kwargs)
		else:
			raise Warning('Not found')

		# nbqp, 
		info = quadRule.getQuadratureRulesInfo()
		quadPtsPos, dersIndices, dersBasis, dersWeights = info
		self._nbqp = len(quadPtsPos); self._qpPar = quadPtsPos; 
		self._basis, self._weights = [], []

		for i in range(2): self._basis.append(array2csr_matrix(dersBasis[:, i], *dersIndices))
		for i in range(4): self._weights.append(array2csr_matrix(dersWeights[:, i], *dersIndices))
		
		return

	def set_DirichletCondition(self, table=[0, 0]):
		dod = []
		dof = [i for i in range(0, self._nbctrlpts)]
		if table[0] == 1:  
			dof.pop(0); dod.append(0)
		if table[-1] == 1: 
			dof.pop(-1); dod.append(self._nbctrlpts-1)
		self._dof = dof
		self._dod = dod
		return 

	def activate_thermal(self):
		prop = self._kwargs.get('property', {})
		self._conductivity = prop.get('conductivity', None)
		self._capacity     = prop.get('capacity', None)
		self._density          = prop.get('density', None)
		return
	
	def activate_mechanical(self):
		prop = self._kwargs.get('property', {})
		self._elasticmodulus = prop.get('elastic_modulus', None)
		self._poissonratio   = prop.get('poisson_ratio', None)
		self._elasticlimit   = prop.get('elastic_limit', None)
		self._density        = prop.get('density', None)
		return
	
	def interpolate_sample(self, u_ctrlpts):
		if np.size(u_ctrlpts, axis=0) != self._nbctrlpts: raise Warning('Not possible')
		nbknots = self._sampleSize
		knots   = np.linspace(0, 1, nbknots)
		basis, indi, indj = evalDersBasisFortran(self._degree, self._knotvector, knots)
		B0       = sp.csr_matrix((basis[:, 0], indj-1, indi-1))
		u_interp = B0.T @ u_ctrlpts
		x_interp = self._detJ * knots
		return u_interp, x_interp

class thermo1D(model1D):
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.activate_thermal()
		self._temporaltheta = kwargs.get('heattheta', 1.0)
		return

	def interpolate_temperature(self, T_ctrlpts):
		" Interpolate temperature in 1D "
		T_interp = self._basis[0].T @ T_ctrlpts
		return T_interp

	def compute_volForce(self, Fprop):
		" Computes 'volumetric' source vector in 1D "
		Fcoefs = Fprop * self._detJ
		F = self._weights[0] @ Fcoefs 
		return F

	def compute_intForce(self, Kprop, Cprop, T, dT):
		"Returns the internal heat force in transient heat"

		# Compute conductivity matrix
		Kcoefs = Kprop * self._invJ
		K = self._weights[-1] @ np.diag(Kcoefs) @ self._basis[1].T 

		# Compute capacity matrix 
		Ccoefs = Cprop * self._detJ
		C = self._weights[0] @ np.diag(Ccoefs) @ self._basis[0].T
		# C = np.diag(C.sum(axis=1))

		# Compute internal heat force 
		Fint = C @ dT + K @ T

		return Fint

	def compute_tangentMatrix(self, Kprop, Cprop, dt):
		""" Computes tangent matrix in transient heat
			S = C + theta dt K
			K = int_Omega dB/dx Kprop dB/dx dx = int_[0, 1] J^-1 dB/dxi Kprop J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
			C = int_Omega B Cprop B dx = int [0, 1] B Cprop det J B dxi
		"""

		# Compute conductivity matrix
		Kcoefs = Kprop * self._invJ
		K = self._weights[-1] @ np.diag(Kcoefs) @ self._basis[1].T 

		# Compute capacitiy matrix 
		Ccoefs = Cprop * self._detJ
		C =self._weights[0] @ np.diag(Ccoefs) @ self._basis[0].T
		# C = np.diag(C.sum(axis=1))

		# Compute tangent matrix 
		M = C + self._temporaltheta*dt*K

		return M

	def solve(self, Fext=None, time_list=None, Tinout=None, threshold=1e-12, nbIterNL=20):
		" Solves transient heat problem in 1D. "

		theta = self._temporaltheta
		dod   = self._dod; dof = self._dof
		VVn0  = np.zeros(len(dof)+len(dod))

		# Compute initial velocity from boundry conditions (for i = 0)
		if np.shape(Tinout)[1] == 2:
			dt = time_list[1] - time_list[0]
			VVn0[dod] = (Tinout[dod, 1] - Tinout[dod, 0])/dt
		elif np.shape(Tinout)[1] >= 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			VVn0[dod] = 1.0/(dt1*(factor-factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
		else: raise Warning('We need more than 2 steps')

		for i in range(1, np.shape(Tinout)[1]):
			# Get delta time
			dt = time_list[i] - time_list[i-1]

			# Get values of last step
			TTn0 = np.copy(Tinout[:, i-1])

			# Get values of new step
			TTn1 = TTn0 + dt*(1-theta)*VVn0; TTn1[dod] = np.copy(Tinout[dod, i])
			TTn1i0 = np.copy(TTn1); VVn1 = np.zeros(len(TTn1))
			VVn1[dod] = 1.0/theta*(1.0/dt*(Tinout[dod, i] - Tinout[dod, i-1]) - (1-theta)*VVn0[dod])
			Fstep = Fext[:, i]

			for j in range(nbIterNL): # Newton-Raphson

				# Interpolate temperature
				T_interp = self.interpolate_temperature(TTn1)

				# Get capacity and conductivity properties
				Kprop = self._conductivity(T_interp)
				Cprop = self._capacity(T_interp)

				# Compute residue
				Fint = self.compute_intForce(Kprop, Cprop, TTn1, VVn1)
				dF = Fstep[dof] - Fint[dof]
				resNL = np.sqrt(np.dot(dF, dF))
				if resNL <= threshold: break

				# Compute tangent matrix
				MM = self.compute_tangentMatrix(Kprop, Cprop, dt=dt)[np.ix_(dof, dof)]

				# Compute delta dT 
				ddVV = np.linalg.solve(MM, dF)

				# Update values
				VVn1[dof] += ddVV
				TTn1[dof] = TTn1i0[dof] + theta*dt*VVn1[dof]

			# print(j+1, relerror)

			# Update values in output
			Tinout[:, i] = np.copy(TTn1)
			VVn0 = np.copy(VVn1)

		return 

class mechamat1D(model1D):
	def __init__(self, **kwargs):
		super().__init__(**kwargs)
		self.activate_mechanical()
		prop = kwargs.get('property', {})
		plasticKwargs  = prop.get('law', None)
		self._isPlasticityPossible = False
		if plasticKwargs is not None: 
			self._isPlasticityPossible = True
			self._mechaBehavLaw        = plasticLaw(self._elasticmodulus, self._elasticlimit, plasticKwargs)
		return

	def returnMappingAlgorithm(self, law:plasticLaw, strain, pls, a, b, threshold=1e-8):
		""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		def computeDeltaGamma(law:plasticLaw, E, eta_trial, a_n0, nbIter=20, threshold=1e-8):
			dgamma = 0.0
			a_n1   = a_n0 
			for i in range(nbIter):
				dH = law._Hfun(a_n1) - law._Hfun(a_n0) 
				G  = -law._Kfun(a_n1) + np.abs(eta_trial) - (E*dgamma + dH)
				if np.abs(G) <=threshold: break
				dG = - (E + law._Hderfun(a_n1) + law._Kderfun(a_n1))
				dgamma -= G/dG
				a_n1   = a_n0 + dgamma 
			return dgamma

		# Elastic predictor
		sigma_trial = self._elasticmodulus*(strain - pls)
		eta_trial   = sigma_trial - b

		# Check yield status
		f_trial = np.abs(eta_trial) - law._Kfun(a)

		if f_trial <= threshold: # Elastic
			stress = sigma_trial
			pls_new = pls
			a_new = a
			b_new = b
			Cep = self._elasticmodulus

		else: # Plastic
			dgamma = computeDeltaGamma(law, self._elasticmodulus, eta_trial, a)
			a_new = a + dgamma
			Normal = np.sign(eta_trial)
			stress = sigma_trial - dgamma*self._elasticmodulus*Normal
			pls_new = pls + dgamma*Normal
			b_new = b + (law._Hfun(a_new) - law._Hfun(a))*Normal
			somme = law._Kderfun(a_new) + law._Hderfun(a_new)
			Cep = self._elasticmodulus*somme/(self._elasticmodulus + somme)

		return [stress, pls_new, a_new, b_new, Cep]

	def interpolate_strain(self, disp):
		" Computes strain field from a given displacement field "
		eps = self._basis[1].T @ disp / self._invJ
		return eps
	
	def interpolate_CPfield(self, u_ref):
		" Interpolate control point field (from parametric to physical space) "
		masse = self._weights[0] @ self._basis[0].T
		force = self._weights[0] @ u_ref
		u_ctrlpts = np.linalg.solve(masse.toarray(), force)
		return u_ctrlpts

	def compute_volForce(self, b): 
		" Computes volumetric force in 1D. Usualy b = rho*g*L (Weight) "
		F = self._weights[0] @ b.T
		return F

	def compute_intForce(self, sigma):
		""" Computes internal force Fint. 
			Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		Fint = self._weights[2] @ sigma.T
		return Fint

	def compute_tangentMatrix(self, Cep):
		""" Computes stiffness matrix in elastoplasticity
			S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		Scoefs = Cep*1.0/self._detJ
		S = self._weights[-1] @ np.diag(Scoefs) @ self._basis[1].T 
		return S

	def solve(self, Fext=None, threshold=1e-8, nbIterNL=10):
		" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

		if not self._isPlasticityPossible: raise Warning('Insert a plastic law')
		law  = self._mechaBehavLaw
		nbqp = self._nbqp
		dof  = self._dof
		pls_n0, a_n0, b_n0 = np.zeros(nbqp), np.zeros(nbqp), np.zeros(nbqp)
		ep_n1, a_n1, b_n1 = np.zeros(nbqp), np.zeros(nbqp), np.zeros(nbqp)
		Cep, sigma = np.zeros(nbqp), np.zeros(nbqp)
		disp = np.zeros(np.shape(Fext))

		strain  = np.zeros((nbqp, np.shape(Fext)[1]))
		plastic = np.zeros((nbqp, np.shape(Fext)[1]))
		stress  = np.zeros((nbqp, np.shape(Fext)[1]))
		moduleE = np.zeros((nbqp, np.shape(Fext)[1]))

		for i in range(1, np.shape(Fext)[1]):

			d_n1 = np.copy(disp[:, i-1]); ddisp = np.zeros(len(dof)) 
			Fstep = np.copy(Fext[:, i])

			print('Step %d' %i)
			for j in range(nbIterNL): # Newton-Raphson

				# Compute strain as function of displacement
				d_n1[dof] = disp[dof, i-1] + ddisp
				eps = self.interpolate_strain(d_n1)

				# Find closest point projection 
				for k in range(nbqp):
					result1 = self.returnMappingAlgorithm(law, eps[k], pls_n0[k], a_n0[k], b_n0[k])
					sigma[k], ep_n1[k], a_n1[k], b_n1[k], Cep[k] = result1

				# Compute Fint
				Fint = self.compute_intForce(sigma)
				dF 	 = Fstep[dof] - Fint[dof]
				relerror = np.sqrt(np.dot(dF, dF))
				print('Rhapson with error %.5e' %relerror)
				if relerror <= threshold: break

				# Compute stiffness
				Sdof = self.compute_tangentMatrix(Cep)[np.ix_(dof, dof)]
				ddisp += np.linalg.solve(Sdof, dF)

			# Update values in output
			disp[:, i]   = d_n1
			strain[:, i] = eps
			stress[:, i] = sigma
			plastic[:, i] = ep_n1
			moduleE[:, i] = Cep

			pls_n0, a_n0, b_n0 = np.copy(ep_n1), np.copy(a_n1), np.copy(b_n1)

		return disp, strain, stress, plastic, moduleE

def plot_results(degree, knotvector, JJ, disp_cp, plastic_cp, stress_cp, folder=None, method='IGA', extension='.png'):
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	knots  = np.linspace(0, 1, 101)
	basis, indi, indj = evalDersBasisFortran(degree, knotvector, knots)
	B0 = sp.csr_matrix((basis[:, 0], indj-1, indi-1)); B1 = sp.csr_matrix((basis[:, -1], indj-1, indi-1))
	displacement   = B0.T @ disp_cp
	strain_interp  = B1.T @ disp_cp
	plastic_interp = B0.T @ plastic_cp
	stress_interp  = B0.T @ stress_cp

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
