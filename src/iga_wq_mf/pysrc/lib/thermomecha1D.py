"""
.. This module contains functions to treat transient heat problems in 1D. 
.. It also contains functions to treat elastoplasticity problems in 1D. 
.. It is adapted to IGA problems (using Gauss quadrature)
.. Joaquin Cornejo
"""

import numpy as np
from lib.__init__ import *
from lib.lib_base import evalDersBasisFortran
from lib.lib_quadrules import GaussQuadrature, WeightedQuadrature, QuadratureRules
from lib.lib_material import plasticLaw

class part1D:
	def __init__(self, kwargs:dict):
		geoArgs  = kwargs.get('geoArgs', {})
		self.__setGeometry(geoArgs)

		quadArgs = kwargs.get('quadArgs', {})
		self.degree     = quadArgs.get('degree', 2)
		self.knotvector = quadArgs.get('knotvector', np.array([0, 0, 0, 0.5, 1, 1, 1]))
		self.nbctrlpts  = len(self.knotvector) - self.degree - 1
		self.__setQuadratureRules(quadArgs)
		return
	
	def __setGeometry(self, geoArgs:dict):
		self.Jqp  = geoArgs.get('length', 1.0)
		self.detJ = self.Jqp
		self.invJ = 1.0/self.Jqp
		return
	
	def __setQuadratureRules(self, quadArgs:dict):
		quadRuleName = quadArgs.get('quadrule', 'wq').lower()
		quadRule     = None
		if quadRuleName == 'iga':
			quadRule = GaussQuadrature(self.degree, self.knotvector, quadArgs=quadArgs)
		elif quadRuleName == 'wq':
			quadRule = WeightedQuadrature(self.degree, self.knotvector, quadArgs=quadArgs)
		else:
			raise Warning('Not found')
		quadRule.getQuadratureRulesInfo()
		self.quadRule = quadRule
		self.basis, self.weights = quadRule.getDenseQuadRules(isFortran=True)	
		return

	def add_DirichletCondition(self, table=[0, 0]):
		dod = []
		dof = [i for i in range(0, self.nbctrlpts)]
		if table[0] == 1:  
			dof.pop(0); dod.append(0)
		if table[-1] == 1: 
			dof.pop(-1); dod.append(self.nbctrlpts-1)
		self.dof = dof
		self.dod = dod
		return 

	def activate_thermal(self, matArgs:dict):
		self.conductivity = matArgs.get('conductivity', None)
		self.capacity     = matArgs.get('capacity', None)
		self.density      = matArgs.get('density', None)
		return
	
	def activate_mechanical(self, matArgs:dict):
		self.elasticmodulus = matArgs.get('elastic_modulus', None)
		self.poissonratio   = matArgs.get('poisson_ratio', None)
		self.elasticlimit   = matArgs.get('elastic_limit', None)
		self.density        = matArgs.get('density', None)
		return
	
	def interpolate_sampleField(self, u_ctrlpts, sampleSize=101):
		if np.size(u_ctrlpts, axis=0) != self.nbctrlpts: raise Warning('Not possible')
		basis, knots = self.quadRule.getSampleBasis(sampleSize=sampleSize)
		u_interp = basis[0].T @ u_ctrlpts
		x_interp = self.detJ * knots
		return u_interp, x_interp
	
	def interpolate_CntrlPtsField(self, u_ref):
		" Interpolate control point field (from parametric to physical space) "
		masse = self.weights[0] @ self.basis[0].T
		force = self.weights[0] @ u_ref
		u_ctrlpts = np.linalg.solve(masse.toarray(), force)
		return u_ctrlpts

class thermo1D(part1D):
	def __init__(self, kwargs:dict):
		super().__init__(kwargs)
		self.activate_thermal(kwargs)
		self._temporaltheta = kwargs.get('heattheta', 1.0)
		return

	def interpolate_temperature(self, T_ctrlpts):
		" Interpolate temperature in 1D "
		T_interp = self.basis[0].T @ T_ctrlpts
		return T_interp

	def compute_volForce(self, Fprop):
		" Computes 'volumetric' source vector in 1D "
		Fcoefs = Fprop * self.detJ
		F = self.weights[0] @ Fcoefs 
		return F

	def compute_intForce(self, Kprop, Cprop, T, dT):
		"Returns the internal heat force in transient heat"

		# Compute conductivity matrix
		Kcoefs = Kprop * self.invJ
		K = self.weights[-1] @ np.diag(Kcoefs) @ self.basis[1].T 

		# Compute capacity matrix 
		Ccoefs = Cprop * self.detJ
		C = self.weights[0] @ np.diag(Ccoefs) @ self.basis[0].T
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
		Kcoefs = Kprop * self.invJ
		K = self.weights[-1] @ np.diag(Kcoefs) @ self.basis[1].T 

		# Compute capacitiy matrix 
		Ccoefs = Cprop * self.detJ
		C =self.weights[0] @ np.diag(Ccoefs) @ self.basis[0].T
		# C = np.diag(C.sum(axis=1))

		# Compute tangent matrix 
		M = C + self._temporaltheta*dt*K

		return M

	def solve(self, Fext=None, time_list=None, Tinout=None, threshold=1e-9, nbIterNL=20):
		" Solves transient heat problem in 1D. "

		theta = self._temporaltheta
		dod   = self.dod; dof = self.dof
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
				Kprop = self.conductivity(T_interp)
				Cprop = self.capacity(T_interp)

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

class mechamat1D(part1D):
	def __init__(self, kwargs:dict):
		super().__init__(kwargs)
		self.activate_mechanical(kwargs)
		plasticVars = kwargs.get('law', None)
		self._isPlasticityPossible = False
		if isinstance(plasticVars, dict):  
			self._isPlasticityPossible = True
			self.mechaBehavLaw        = plasticLaw(self.elasticmodulus, self.elasticlimit, plasticVars)
		return

	def returnMappingAlgorithm(self, law:plasticLaw, strain, pls, a, b, threshold=1e-9):
		""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		def computeDeltaGamma(law:plasticLaw, E, eta_trial, a_n0, nbIter=20, threshold=1e-9):
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
		sigma_trial = self.elasticmodulus*(strain - pls)
		eta_trial   = sigma_trial - b

		# Check yield status
		f_trial = np.abs(eta_trial) - law._Kfun(a)
		stress = sigma_trial
		pls_new = pls
		a_new = a
		b_new = b
		Cep = self.elasticmodulus

		if f_trial > threshold: # Plastic
			dgamma = computeDeltaGamma(law, self.elasticmodulus, eta_trial, a)
			a_new = a + dgamma
			Normal = np.sign(eta_trial)
			stress = sigma_trial - dgamma*self.elasticmodulus*Normal
			pls_new = pls + dgamma*Normal
			b_new = b + (law._Hfun(a_new) - law._Hfun(a))*Normal
			somme = law._Kderfun(a_new) + law._Hderfun(a_new)
			Cep = self.elasticmodulus*somme/(self.elasticmodulus + somme)

		return [stress, pls_new, a_new, b_new, Cep]

	def interpolate_strain(self, disp):
		" Computes strain field from a given displacement field "
		eps = self.basis[1].T @ disp / self.invJ
		return eps

	def compute_volForce(self, b): 
		" Computes volumetric force in 1D. Usualy b = rho*g*L (Weight) "
		F = self.weights[0] @ b.T
		return F

	def compute_intForce(self, sigma):
		""" Computes internal force Fint. 
			Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		Fint = self.weights[2] @ sigma.T
		return Fint

	def compute_tangentMatrix(self, Cep):
		""" Computes stiffness matrix in elastoplasticity
			S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		Scoefs = Cep*1.0/self.detJ
		S = self.weights[-1] @ np.diag(Scoefs) @ self.basis[1].T 
		return S

	def solve(self, Fext=None, threshold=1e-9, nbIterNL=30):
		" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

		if not self._isPlasticityPossible: raise Warning('Insert a plastic law')
		law  = self.mechaBehavLaw
		nbqp = self.quadRule.nbqp
		dof  = self.dof
		pls_n0, a_n0, b_n0 = np.zeros(nbqp), np.zeros(nbqp), np.zeros(nbqp)
		ep_n1, a_n1, b_n1 = np.zeros(nbqp), np.zeros(nbqp), np.zeros(nbqp)
		Cep, sigma = np.zeros(nbqp), np.zeros(nbqp)
		disp = np.zeros(np.shape(Fext))

		eps     = np.zeros((nbqp, np.shape(Fext)[1]))
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
				strain = self.interpolate_strain(d_n1)

				# Find closest point projection 
				for k in range(nbqp):
					tmp = self.returnMappingAlgorithm(law, strain[k], pls_n0[k], a_n0[k], b_n0[k])
					sigma[k], ep_n1[k], a_n1[k], b_n1[k], Cep[k] = tmp

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
			eps[:, i]    = strain
			stress[:, i] = sigma
			plastic[:, i] = ep_n1
			moduleE[:, i] = Cep

			pls_n0, a_n0, b_n0 = np.copy(ep_n1), np.copy(a_n1), np.copy(b_n1)

		return disp, eps, stress, plastic, moduleE

def plot_results(quadRule:QuadratureRules, JJ, disp_cp, plastic_cp, stress_cp, folder=None, method='IGA', extension='.png'):
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	basis, knots = quadRule.getSampleBasis(sampleSize=101)
	displacement   = basis[0].T @ disp_cp
	strain_interp  = basis[1].T @ disp_cp
	plastic_interp = basis[0].T @ plastic_cp
	stress_interp  = basis[0].T @ stress_cp

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

	# # Plot stress-strain of single point
	# fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
	# for ax, pos in zip([ax1, ax2, ax3], [25, 50, 75]):
	# 	ax.plot(strain_interp[pos, :]*100, stress_interp[pos, :])
	# 	ax.set_ylabel('Stress (MPa)')
	# 	ax.set_xlabel('Strain (\%)')
	# 	ax.set_ylim(bottom=0.0, top=1500)
	# 	ax.set_xlim(left=0.0, right=strain_interp.max()*100)

	# fig.tight_layout()
	# fig.savefig(folder + 'TractionCurve' + method + extension)
	return
