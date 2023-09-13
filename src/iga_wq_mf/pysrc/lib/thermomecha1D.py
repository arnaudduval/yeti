"""
.. This module contains functions to treat transient heat problems in 1D. 
.. It also contains functions to treat elastoplasticity problems in 1D. 
.. Joaquin Cornejo
"""

from .__init__ import *
from .lib_quadrules import GaussQuadrature, WeightedQuadrature, QuadratureRules
from .lib_material import plasticLaw
from .lib_base import evalDersBasisPy

class part1D:
	def __init__(self, kwargs:dict):
		geoArgs  = kwargs.get('geoArgs', {})
		self.__setGeometry(geoArgs)

		quadArgs = kwargs.get('quadArgs', {})
		self.degree     = quadArgs.get('degree', 2)
		self.knotvector = quadArgs.get('knotvector', np.array([0, 0, 0, 0.5, 1, 1, 1]))
		self.nbctrlpts  = len(self.knotvector) - self.degree - 1
		self.__setQuadratureRules(quadArgs)
		self.addSolverConstraints(kwargs.get('solverArgs', {}))
		return
	
	def addSolverConstraints(self, solverArgs:dict):
		self._thresholdNR = solverArgs.get('NRThreshold', 1.e-10)
		self._nbIterNR    = solverArgs.get('nbIterationsNR', 50)
		return
	
	def __setGeometry(self, geoArgs:dict):
		self.Jqp  = geoArgs.get('length', 1.0)
		self.detJ = self.Jqp
		self.invJ = 1.0/self.Jqp
		return
	
	def __setQuadratureRules(self, quadArgs:dict):
		quadRuleName = quadArgs.get('quadrule', '').lower()
		if quadRuleName == 'iga':
			quadRule = GaussQuadrature(self.degree, self.knotvector, quadArgs=quadArgs)
		elif quadRuleName == 'wq':
			quadRule = WeightedQuadrature(self.degree, self.knotvector, quadArgs=quadArgs)
		else: raise Warning('Not found')
		quadRule.getQuadratureRulesInfo()
		self.quadRule = quadRule
		self.basis, self.weights = quadRule.getDenseQuadRules(isFortran=True)	
		self.qpPhy = quadRule.quadPtsPos*self.Jqp
		self.nbqp  = len(self.qpPhy)
		return

	def add_DirichletCondition(self, table=[0, 0]):
		dod = np.array([], dtype=int)
		dof = np.arange(self.nbctrlpts, dtype=int)
		if table[0] == 1:  
			dof = np.delete(dof, 0); dod = np.append(dod, 0)
		if table[-1] == 1: 
			dof = np.delete(dof, -1); dod = np.append(dod, self.nbctrlpts-1)
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
	
	def interpolateMeshgridField(self, u_ctrlpts, sampleSize=101, isSample=True):
		if isSample: basis, knots = self.quadRule.getSampleBasis(sampleSize=sampleSize)
		else: 		 basis 		  = self.quadRule.getDenseQuadRules()[0]
		u_interp = basis[0].T @ u_ctrlpts
		if isSample: x_interp = self.detJ * knots
		else: 		 x_interp = self.qpPhy
		return u_interp, x_interp
	
	def L2projectionCtrlpts(self, u_ref):
		" Interpolate control point field (from parametric to physical space) "
		masse = self.weights[0] @ np.diag(self.detJ*np.ones(self.nbqp)) @ self.basis[0].T
		force = self.weights[0] @ np.diag(self.detJ*np.ones(self.nbqp)) @ u_ref
		masseSparse = sp.csr_matrix(masse)
		u_ctrlpts = sp.linalg.spsolve(masseSparse, force)
		return u_ctrlpts
	
	def L2NormOfError(self, u_ctrlpts, L2NormArgs:dict):
		""" Computes the norm L2 of the error. The exactfun is the function of the exact solution. 
			and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
			whether the default quadrature is weighted quadrature. 
		"""	
		
		# Compute u interpolation
		quadRule = GaussQuadrature(self.degree, self.knotvector, quadArgs={'type':'leg'})
		quadPts  = quadRule.getQuadratureRulesInfo()[0]
		denseBasis = quadRule.getDenseQuadRules()[0]
		parametricWeights = quadRule._parametricWeights
		
		qpPhy = quadPts*self.detJ
		u_interp = denseBasis[0].T @ u_ctrlpts

		# Compute u exact
		u_exact  = None
		exactfun = L2NormArgs.get('exactFunction', None)
		if callable(exactfun): u_exact = exactfun(qpPhy)
		refPart  = L2NormArgs.get('referencePart', None); u_ref = L2NormArgs.get('u_ref', None)
		if isinstance(refPart, part1D) and isinstance(u_ref, np.ndarray):
			denseBasisExact = evalDersBasisPy(refPart.degree, refPart.knotvector, quadPts)
			u_exact = denseBasisExact[0].T @ u_ref
		if u_exact is None: raise Warning('Not possible')

		# Compute error
		ue_diff_uh2 = (u_exact - u_interp)**2 * self.detJ
		ue2         = u_exact**2 * self.detJ

		tmp1 = np.einsum('i,i->', parametricWeights, ue2)
		tmp2 = np.einsum('i,i->', parametricWeights, ue_diff_uh2)
		error = np.sqrt(tmp2/tmp1)
		return error

class thermo1D(part1D):
	def __init__(self, kwargs:dict):
		super().__init__(kwargs)	
		return
	
	def activate_thermal(self, matArgs:dict):
		super().activate_thermal(matArgs)
		self._temporaltheta = matArgs.get('heattheta', 1.0)
		return 

	def interpolate_temperature(self, T_ctrlpts):
		" Interpolate temperature in 1D "
		T_interp = self.basis[0].T @ T_ctrlpts
		return T_interp

	def compute_volForce(self, Fprop):
		" Computes 'volumetric' source vector in 1D "
		Fcoefs = Fprop*self.detJ
		Fvol = self.weights[0] @ Fcoefs 
		return Fvol

	def compute_intForce(self, Kprop, Cprop, T, dT, isLumped=False):
		"Returns the internal heat force in transient heat"

		Kcoefs = Kprop*self.invJ
		K = self.weights[-1] @ np.diag(Kcoefs) @ self.basis[1].T 
		Ccoefs = Cprop * self.detJ
		C = self.weights[0] @ np.diag(Ccoefs) @ self.basis[0].T
		if isLumped: C = np.diag(C.sum(axis=1))
		Fint = C @ dT + K @ T

		return Fint

	def compute_tangentMatrix(self, Kprop, Cprop, dt, isLumped=False):
		""" Computes tangent matrix in transient heat
			S = C + theta dt K
			K = int_Omega dB/dx Kprop dB/dx dx = int_[0, 1] J^-1 dB/dxi Kprop J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
			C = int_Omega B Cprop B dx = int [0, 1] B Cprop det J B dxi
		"""

		Kcoefs = Kprop * self.invJ
		K = self.weights[-1] @ np.diag(Kcoefs) @ self.basis[1].T 
		Ccoefs = Cprop * self.detJ
		C =self.weights[0] @ np.diag(Ccoefs) @ self.basis[0].T
		if isLumped: C = np.diag(C.sum(axis=1))
		M = C + self._temporaltheta*dt*K

		return M

	def solve(self, Fext=None, time_list=None, Tinout=None, isLumped=False):
		" Solves transient heat problem in 1D. "

		theta = self._temporaltheta
		dod   = self.dod; dof = self.dof
		V_n0  = np.zeros(self.nbctrlpts)

		# Compute initial velocity using interpolation
		if np.shape(Tinout)[1] == 2:
			dt = time_list[1] - time_list[0]
			V_n0[dod] = (Tinout[dod, 1] - Tinout[dod, 0])/dt
		elif np.shape(Tinout)[1] >= 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			V_n0[dod] = 1.0/(dt1*(factor-factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
		else: raise Warning('We need more than 2 steps')

		for i in range(1, np.shape(Tinout)[1]):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]

			# Get values of last step
			d_n0 = np.copy(Tinout[:, i-1])

			# Get values of new step
			d0_n1 = d_n0 + dt*(1 - theta)*V_n0
			V_n1 = np.zeros(self.nbctrlpts); V_n1[dod] = 1.0/theta*(1.0/dt*(Tinout[dod, i] - Tinout[dod, i-1]) - (1 - theta)*V_n0[dod])
			Fext_n1 = Fext[:, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR): 
				
				# Interpolate temperature
				dj_n1 = d0_n1 + theta*dt*V_n1
				temperature = self.interpolate_temperature(dj_n1)

				# Compute internal force
				Kprop = self.conductivity(temperature)
				Cprop = self.capacity(temperature)
				Fint_dj = self.compute_intForce(Kprop, Cprop, dj_n1, V_n1, isLumped=isLumped)

				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0
				
				# Direct solver
				tangentM = sp.csr_matrix(self.compute_tangentMatrix(Kprop, Cprop, dt=dt, isLumped=isLumped)[np.ix_(dof, dof)])
				deltaV = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update values
				V_n1[dof] += deltaV

				resNRj = abs(theta*dt*np.dot(V_n1, r_dj))
				if j == 0: resNR0 = resNRj
				print('Rhapson with error %.5e' %resNRj)
				if j > 0 and resNRj <= self._thresholdNR*resNR0: break

			# Update values in output
			Tinout[:, i] = np.copy(dj_n1)
			V_n0 = np.copy(V_n1)

		return 

class mechamat1D(part1D):
	def __init__(self, kwargs:dict):
		super().__init__(kwargs)		
		return
	
	def activate_mechanical(self, matArgs:dict):
		super().activate_mechanical(matArgs)
		self.plasticLaw = None
		self._isPlasticityPossible = False
		tmp = matArgs.get('plasticLaw', None)
		if isinstance(tmp, dict):  
			self._isPlasticityPossible = True
			self.plasticLaw            = plasticLaw(self.elasticmodulus, self.elasticlimit, tmp)
		return 

	def returnMappingAlgorithm(self, strain, pls, a, b, threshold=1e-9):
		""" Return mapping algorithm for one-dimensional rate-independent plasticity. 
			It uses combined isotropic/kinematic hardening theory.  
		"""

		def computeDeltaGamma(law:plasticLaw, E, a_n0, eta_trial, nbIter=100, threshold=1e-10):
			dgamma = np.zeros(np.size(a_n0))
			a_n1   = a_n0 
			for i in range(nbIter):
				dH = law._Hfun(a_n1) - law._Hfun(a_n0) 
				G  = -law._Kfun(a_n1) + np.abs(eta_trial) - (E*dgamma + dH)
				if np.all(np.abs(G)<=threshold): break
				dG = - (E + law._Hderfun(a_n1) + law._Kderfun(a_n1))
				dgamma -= G/dG
				a_n1 = a_n0 + dgamma
			return dgamma

		nnz = np.size(strain)
		output = np.zeros((5, nnz))

		# Compute trial stress
		s_trial = self.elasticmodulus*(strain - pls)
		
		# Computed shifted stress 
		eta_trial = s_trial - b

		# Check yield status
		f_trial = np.abs(eta_trial) - self.plasticLaw._Kfun(a)
		Cep   = self.elasticmodulus*np.ones(nnz)
		pls_new = pls
		a_new = a
		b_new = b
		stress  = s_trial

		plsInd = np.nonzero(f_trial>threshold)[0]
		if np.size(plsInd) > 0:

			# Compute plastic-strain increment
			dgamma_plsInd = computeDeltaGamma(self.plasticLaw, self.elasticmodulus, a[plsInd], eta_trial[plsInd])
			
			# Update internal hardening variable
			a_new[plsInd] = a[plsInd] + dgamma_plsInd
			
			# Compute df/dsigma
			Normal_plsInd = np.sign(eta_trial[plsInd])

			# Update stress
			stress[plsInd] = s_trial[plsInd] - dgamma_plsInd*self.elasticmodulus*Normal_plsInd
			
			# Update plastic strain
			pls_new[plsInd] = pls[plsInd] + dgamma_plsInd*Normal_plsInd

			# Update backstress
			b_new[plsInd] = b[plsInd] + (self.plasticLaw._Hfun(a_new[plsInd]) - self.plasticLaw._Hfun(a[plsInd]))*Normal_plsInd
			
			# Update tangent coefficients
			somme = self.plasticLaw._Kderfun(a_new[plsInd]) + self.plasticLaw._Hderfun(a_new[plsInd])
			Cep[plsInd] = self.elasticmodulus*somme/(self.elasticmodulus + somme)

		output[0, :] = stress; output[1, :] = pls_new; output[2, :] = a_new
		output[3, :] = b_new; output[4, :] = Cep
		return output

	def interpolate_strain(self, disp):
		" Computes strain field from a given displacement field "
		strain = self.basis[1].T @ disp * self.invJ
		return strain

	def compute_volForce(self, b): 
		" Computes volumetric force in 1D. "
		Fcoefs = self.detJ*np.ones(self.nbqp)
		Fvol = self.weights[0] @ np.diag(Fcoefs) @ b.T
		return Fvol

	def compute_intForce(self, sigma):
		""" Computes internal force Fint. 
			Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		Fint = self.weights[-1] @ sigma.T
		return Fint

	def compute_tangentMatrix(self, Cep):
		""" Computes stiffness matrix in elastoplasticity
			S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		Scoefs = Cep*self.invJ
		S = self.weights[-1] @ np.diag(Scoefs) @ self.basis[1].T 
		return S

	def solve(self, Fext=None):
		" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

		if not self._isPlasticityPossible: raise Warning('Insert a plastic law')
		nbqp = self.nbqp; dof = self.dof; dod = self.dod
		pls_n0, a_n0, b_n0 = np.zeros(nbqp), np.zeros(nbqp), np.zeros(nbqp)
		pls_n1, a_n1, b_n1 = np.zeros(nbqp), np.zeros(nbqp), np.zeros(nbqp)
		stress, Cep = np.zeros(nbqp), np.zeros(nbqp)
		Alldisplacement = np.zeros(np.shape(Fext))

		Allstrain  = np.zeros((nbqp, np.shape(Fext)[1]))
		Allplastic = np.zeros((nbqp, np.shape(Fext)[1]))
		Allstress  = np.zeros((nbqp, np.shape(Fext)[1]))
		AllCep = np.zeros((nbqp, np.shape(Fext)[1]))

		for i in range(1, np.shape(Fext)[1]):
			
			# Get values of last step
			d_n0 = np.copy(Alldisplacement[:, i-1])
			
			# Get values of new step
			V_n1 = np.zeros(self.nbctrlpts) 
			Fext_n1 = np.copy(Fext[:, i])

			print('Step: %d' %i)
			for j in range(self._nbIterNR): # Newton-Raphson 
				
				# Compute strain at each quadrature point
				dj_n1 = d_n0 + V_n1
				strain = self.interpolate_strain(dj_n1)

				# Find closest point projection 
				output = self.returnMappingAlgorithm(strain, pls_n0, a_n0, b_n0, threshold=1e-10)
				stress, pls_n1, a_n1, b_n1, Cep = output[0, :], output[1, :], output[2, :], output[3, :], output[4, :]

				# Compute internal force 
				Fint_dj = self.compute_intForce(stress)
				
				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				# Direct solver
				tangentS = sp.csr_matrix(self.compute_tangentMatrix(Cep)[np.ix_(dof, dof)])
				deltaV = sp.linalg.spsolve(tangentS, r_dj[dof])
				
				# Update values
				V_n1[dof] += deltaV

				# Compute residue of Newton Raphson using an energetic approach
				resNRj = abs(np.dot(V_n1, r_dj))
				if j == 0: resNR0 = resNRj
				print('Rhapson with error %.5e' %resNRj)
				if j > 0 and resNRj <= self._thresholdNR*resNR0: break

			Alldisplacement[:, i] = dj_n1
			Allstrain[:, i] = strain
			Allstress[:, i] = stress
			Allplastic[:, i] = pls_n1
			AllCep[:, i] = Cep

			pls_n0, a_n0, b_n0 = np.copy(pls_n1), np.copy(a_n1), np.copy(b_n1)

		return Alldisplacement, Allstrain, Allstress, Allplastic, AllCep

def plot_results(quadRule:QuadratureRules, JJ, disp_cp, plastic_cp, stress_cp, folder=None, method='iga', extension='.png'):
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	basis, knots = quadRule.getSampleBasis(sampleSize=101)
	displacement   = basis[0].T @ disp_cp
	strain_interp  = basis[1].T @ disp_cp / JJ
	plastic_interp = basis[0].T @ plastic_cp
	stress_interp  = basis[0].T @ stress_cp

	# Plot fields
	N = np.shape(disp_cp)[1]
	XX, STEPS = np.meshgrid(knots*JJ, np.arange(N))
	names = ['Displacement field', 'Plastic strain field', 'Stress field']
	units = ['mm', '', 'MPa']
	fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
	for ax, variable, name, unit in zip([ax1, ax2, ax3], [displacement, plastic_interp, stress_interp], names, units):
		im = ax.pcolormesh(XX, STEPS, variable.T, cmap='PuBu_r', shading='linear')
		ax.set_title(name)
		ax.set_ylabel('Step')
		ax.set_xlabel('Position (mm)')
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
		ax.plot(strain_interp[pos, :], stress_interp[pos, :])
		ax.set_ylabel('Stress (MPa)')
		ax.set_xlabel('Total strain (-)')
		ax.set_ylim(bottom=0.0, top=200)
		ax.set_xlim(left=0.0, right=strain_interp.max())

	fig.tight_layout()
	fig.savefig(folder + 'TractionCurve' + method + extension)
	return
