"""
.. This module contains functions to treat transient heat problems in 1D. 
.. It also contains functions to treat elastoplasticity problems in 1D. 
.. Joaquin Cornejo
"""

from . import *
from .lib_part import part1D
from .lib_quadrules import GaussQuadrature
from .lib_material import mechamat
from .lib_base import evalDersBasisFortran, array2csr_matrix

class problem1D(part1D):
	def __init__(self, part:BSpline.Curve, kwargs:dict):
		
		part1D.__init__(self, part, kwargs)
		self.addSolverConstraints(kwargs.get('solverArgs', {}))
		return
	
	def addSolverConstraints(self, solverArgs:dict):
		self._thresholdNewton = solverArgs.get('thresholdNewton', 1e-6)
		self._nIterNewton     = solverArgs.get('nIterationsNewton', 20)
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
		self.relaxation   = matArgs.get('relaxation', None)
		return
	
	def compute_volForce(self, volfun):
		" Computes 'volumetric' source vector in 1D "
		prop = volfun(self.qpPhy)*self.detJ
		volForce = self._denseweights[0] @ prop 
		return volForce
	
	def interpolateMeshgridField(self, u_ctrlpts, sampleSize=101, isSample=True):
		if isSample: basis = self.quadRule.getSampleBasis(sampleSize=sampleSize)[0]
		else: 		 basis 		  = self.quadRule.getDenseQuadRules()[0]
		u_interp = basis[0].T @ u_ctrlpts
		if isSample: x_interp = basis[0].T @ self.ctrlpts
		else: 		 x_interp = self.qpPhy
		return u_interp, x_interp
	
	def L2projectionCtrlpts(self, u_atqp):
		" Interpolate control point field (from parametric to physical space) "
		masse = self._denseweights[0] @ np.diag(self.detJ) @ self._densebasis[0].T
		force = self._denseweights[0] @ np.diag(self.detJ) @ u_atqp
		massesp   = sp.csr_matrix(masse)
		u_ctrlpts = sp.linalg.spsolve(massesp, force)
		return u_ctrlpts
	
	def normOfError(self, u_ctrlpts, normArgs:dict, isRelative=True):
		""" Computes the norm L2 or H1 of the error. The exactfun is the function of the exact solution. 
			and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
			whether the default quadrature is weighted quadrature. 
		"""	
		typeNorm = normArgs.get('type', 'l2').lower()
		if all(norm != typeNorm  for norm in ['l2', 'h1', 'semih1']): raise Warning('Unknown norm')

		# Compute u interpolation
		quadRule = GaussQuadrature(self.degree, self.knotvector, quadArgs={'type':'leg'})
		quadPts  = quadRule.getQuadratureRulesInfo()[0]
		denseBasis = quadRule.getDenseQuadRules()[0]
		parametricWeights = quadRule._parametricWeights
		
		qpPhy = denseBasis[0].T @ self.ctrlpts 
		detJ  = denseBasis[1].T @ self.ctrlpts
		u_interp = denseBasis[0].T @ u_ctrlpts
		uders_interp = denseBasis[1].T @ u_ctrlpts / detJ

		# Compute u exact
		u_exact, uders_exact = None, None
		exactfun = normArgs.get('exactFunction', None)
		exactfunders = normArgs.get('exactFunctionDers', None)
		if callable(exactfun): u_exact = exactfun(qpPhy)
		if callable(exactfunders): uders_exact = exactfunders(qpPhy)
		
		part_ref = normArgs.get('part_ref', None); u_ref = normArgs.get('u_ref', None)
		if isinstance(part_ref, problem1D) and isinstance(u_ref, np.ndarray):
			denseBasisExact = []
			basis_csr, indi_csr, indj_csr = evalDersBasisFortran(part_ref.degree, part_ref.knotvector, quadPts)
			for i in range(2): denseBasisExact.append(array2csr_matrix(basis_csr[:, i], indi_csr, indj_csr))
			u_exact = denseBasisExact[0].T @ u_ref
			invJExact = denseBasisExact[1].T @ part_ref.ctrlpts
			uders_exact = denseBasisExact[1].T @ u_ref / invJExact
			
		# Compute error
		uedfuf2_l2, uedfuf2_sh1 = 0., 0.
		ue2_l2, ue2_sh1 = 0., 0.

		if typeNorm == 'l2' or typeNorm == 'h1':
			uedfuf2_l2 += (u_exact - u_interp)**2
			ue2_l2     += u_exact**2

		if typeNorm == 'h1' or typeNorm == 'semih1':
			uedfuf2_sh1 += (uders_exact - uders_interp)**2
			ue2_sh1     += uders_exact**2

		norm1 = (uedfuf2_l2 + uedfuf2_sh1)*detJ
		norm2 = (ue2_l2 + ue2_sh1)*detJ

		tmp1 = np.einsum('i,i->', parametricWeights, norm1)
		tmp2 = np.einsum('i,i->', parametricWeights, norm2)

		if isRelative: error = np.sqrt(tmp1/tmp2)
		else:          error = np.sqrt(tmp1)
		return error

class heatproblem1D(problem1D):
	def __init__(self, part, kwargs:dict):
		problem1D.__init__(self, part, kwargs)	
		return
	
	def activate_thermal(self, matArgs:dict):
		problem1D.activate_thermal(self, matArgs)
		if self.density is None: self.density = lambda x: np.ones(len(x))
		return 
	
	def compute_mfCapacity(self, Cprop, array_in, isLumped=False):
		Ccoefs = Cprop*self.detJ
		C = self._denseweights[0] @ np.diag(Ccoefs) @ self._densebasis[0].T
		if isLumped: C = np.diag(C.sum(axis=1))
		array_out = C @ array_in
		return array_out
	
	def compute_mfConductivity(self, Kprop, array_in):
		Kcoefs = Kprop*self.invJ
		K = self._denseweights[-1] @ np.diag(Kcoefs) @ self._densebasis[1].T 
		array_out = K @ array_in
		return array_out

	def interpolate_temperature(self, T_ctrlpts):
		" Interpolate temperature in 1D "
		T_interp = self._densebasis[0].T @ T_ctrlpts
		return T_interp

	def compute_FourierMatrix(self, Kprop, Cprop, dt, alpha=1.0, isLumped=False):
		""" Computes tangent matrix in transient heat
			S = C + theta dt K
			K = int_Omega dB/dx Kprop dB/dx dx = int_[0, 1] J^-1 dB/dxi Kprop J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
			C = int_Omega B Cprop B dx = int [0, 1] B Cprop det J B dxi
		"""

		Kcoefs = Kprop*self.invJ
		K = self._denseweights[-1] @ np.diag(Kcoefs) @ self._densebasis[1].T 
		Ccoefs = Cprop*self.detJ
		C = self._denseweights[0] @ np.diag(Ccoefs) @ self._densebasis[0].T
		if isLumped: C = np.diag(C.sum(axis=1))
		tangentM = C + alpha*dt*K

		return tangentM

	def compute_CattaneoMatrix(self, Kprop, Cprop, Mprop, dt, beta=0.25, gamma=0.5, isLumped=False):
		Kcoefs = Kprop*self.invJ
		K = self._denseweights[-1] @ np.diag(Kcoefs) @ self._densebasis[1].T 
		Ccoefs = Cprop*self.detJ
		C = self._denseweights[0] @ np.diag(Ccoefs) @ self._densebasis[0].T
		Mcoefs = Mprop*self.detJ
		M = self._denseweights[0] @ np.diag(Mcoefs) @ self._densebasis[0].T
		if isLumped: M = np.diag(M.sum(axis=1))
		tangentM = M + gamma*dt*C + beta*dt**2*K

		return tangentM

	def solveFourierTransientProblem(self, Tinout, Fext_list, time_list, alpha=1.0, isLumped=False):
		" Solves transient heat problem in 1D. "

		dod  = self.dod; dof = self.dof
		V_n0 = np.zeros(self.nbctrlpts)

		# Compute initial velocity using interpolation
		if np.shape(Tinout)[1] >= 2:
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

			# Predict values of new step
			dj_n1 = d_n0 + (1 - alpha)*dt*V_n0
			Vj_n1 = np.zeros(self.nbctrlpts)
			
			# Overwrite inactive control points
			Vj_n1[dod] = 1.0/(alpha*dt)*(Tinout[dod, i] - dj_n1[dod])
			dj_n1[dod] = Tinout[dod, i]

			Fext_n1 = np.copy(Fext_list[:, i])

			print('Step: %d' %i)
			for j in range(self._nIterNewton): 
				
				# Interpolate temperature
				temperature = self.interpolate_temperature(dj_n1)

				# Compute internal force
				Kprop = self.conductivity(temperature)
				Cprop = self.capacity(temperature)*self.density(temperature)
				Fint_dj = self.compute_mfCapacity(Cprop, Vj_n1, isLumped=isLumped) + self.compute_mfConductivity(Kprop, dj_n1)

				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				resNRj = np.sqrt(np.dot(Vj_n1, r_dj))
				if j == 0: resNR0 = resNRj
				print('NR error %.5e' %resNRj)
				if resNRj <= max([self._thresholdNewton*resNR0, 1e-14]): break
				
				# Solver for active control points
				tangentM = sp.csr_matrix(self.compute_FourierMatrix(Kprop, Cprop, dt=dt, alpha=alpha, isLumped=isLumped)[np.ix_(dof, dof)])
				deltaV = np.zeros(self.nbctrlpts); deltaV[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])
				
				# Update active control points
				dj_n1 += alpha*dt*deltaV
				Vj_n1 += deltaV

			Tinout[:, i] = np.copy(dj_n1)
			V_n0 = np.copy(Vj_n1)

		return 

	def solveCattaneoTransientProblem(self, Tinout, Fext_list, time_list, beta=0.25, gamma=0.5, isLumped=False):
		" Solves transient heat problem in 1D using Cattaneo approach. "

		dod  = self.dod; dof = self.dof
		V_n0 = np.zeros(self.nbctrlpts)
		A_n0 = np.zeros(self.nbctrlpts)

		# Compute initial velocity using interpolation
		if np.shape(Tinout)[1] >= 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			V_n0[dod] = 1.0/(dt1*(factor-factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
			A_n0[dod] = 2.0/(dt1*dt2)*((Tinout[dod, 2] - factor*Tinout[dod, 1])/(factor - 1) + Tinout[dod, 0])
		else: raise Warning('We need more than 2 steps')

		for i in range(1, np.shape(Tinout)[1]):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]

			# Get values of last step
			d_n0 = np.copy(Tinout[:, i-1])

			# Predict values of new step
			dj_n1 = d_n0 + dt*V_n0 + 0.5*dt**2*(1 - 2*beta)*A_n0
			Vj_n1 = V_n0 + (1 - gamma)*dt*A_n0
			Aj_n1 = np.zeros(self.nbctrlpts) 

			# Overwrite inactive control points
			Aj_n1[dod] = (Tinout[dod, i] - dj_n1[dod])/(beta*dt**2)
			Vj_n1[dod] = Vj_n1[dod] + gamma*dt*Aj_n1[dod]
			dj_n1[dod] = Tinout[dod, i]

			Fext_n1 = np.copy(Fext_list[:, i])

			print('Step: %d' %i)
			for j in range(self._nIterNewton): 
				
				# Interpolate temperature
				temperature = self.interpolate_temperature(dj_n1)

				# Compute internal force
				Kprop = self.conductivity(temperature)
				Cprop = self.capacity(temperature)*self.density(temperature)
				Mprop = self.relaxation(temperature)*Cprop
				Fint_dj = (self.compute_mfCapacity(Mprop, Aj_n1, isLumped=isLumped) 
						+ self.compute_mfCapacity(Cprop, Vj_n1, isLumped=False) 
						+ self.compute_mfConductivity(Kprop, dj_n1)
				)

				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				resNRj = np.sqrt(np.dot(r_dj, r_dj))
				if j == 0: resNR0 = resNRj
				print('NR error %.5e' %resNRj)
				if resNRj <= max([self._thresholdNewton*resNR0, 1e-14]): break

				# Solve for active control points
				tangentM = sp.csr_matrix(self.compute_CattaneoMatrix(Kprop, Cprop, Mprop, dt=dt, 
										beta=beta, gamma=gamma, isLumped=isLumped)[np.ix_(dof, dof)])
				deltaA = np.zeros(self.nbctrlpts); deltaA[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update active control points
				dj_n1 += beta*dt**2*deltaA
				Vj_n1 += gamma*dt*deltaA
				Aj_n1 += deltaA

			Tinout[:, i] = np.copy(dj_n1)
			V_n0 = np.copy(Vj_n1)
			A_n0 = np.copy(Aj_n1)

		return

class mechaproblem1D(problem1D):
	def __init__(self, part, kwargs:dict):
		problem1D.__init__(self, part, kwargs)		
		return
	
	def activate_mechanical(self, mechamaterial:mechamat):
		assert isinstance(mechamaterial, mechamat)
		self.mechamat = mechamaterial
		return 

	def interpolate_strain(self, disp):
		" Computes strain field from a given displacement field "
		strain = self._densebasis[1].T @ disp * self.invJ
		strain = np.reshape(strain, (1, len(strain)))
		return strain

	def compute_MechStaticIntForce(self, stress):
		""" Computes internal force Fint. 
			Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		Fint = self._denseweights[-1] @ np.ravel(stress).T
		return Fint

	def compute_tangentMatrix(self, Cep):
		""" Computes stiffness matrix in elastoplasticity
			S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		coefs = np.ravel(Cep)*self.invJ
		tangentM = self._denseweights[-1] @ np.diag(coefs) @ self._densebasis[1].T 
		return tangentM

	def solvePlasticityProblem(self, dispinout, Fext_list):
		" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

		nbChaboche = self.mechamat._chabocheNBparameters
		nbqp = self.nbqp; dof = self.dof; dod = self.dod
		pls_n0, a_n0, b_n0 = np.zeros((1, nbqp)), np.zeros((1, nbqp)), np.zeros((nbChaboche, 1, nbqp))
		pls_n1, a_n1, b_n1 = np.zeros((1, nbqp)), np.zeros((1, nbqp)), np.zeros((nbChaboche, 1, nbqp))

		Allstrain  = np.zeros((nbqp, np.shape(Fext_list)[1]))
		Allplseq = np.zeros((nbqp, np.shape(Fext_list)[1]))
		Allstress  = np.zeros((nbqp, np.shape(Fext_list)[1]))
		AllCep = np.zeros((nbqp, np.shape(Fext_list)[1]))

		for i in range(1, np.shape(Fext_list)[1]):
			
			# Get values of last step
			d_n0 = np.copy(dispinout[:, i-1])
			
			# Predict values of new step
			dj_n1 = np.copy(d_n0)

			# Overwrite inactive control points
			dj_n1[dod] = np.copy(dispinout[dod, i])

			Fext_n1 = np.copy(Fext_list[:, i])

			print('Step: %d' %i)
			for j in range(100): # Newton-Raphson 
				
				# Compute strain at each quadrature point
				strain = self.interpolate_strain(dj_n1)

				# Find closest point projection 
				output, isElasticLoad = self.mechamat.J2returnMappingAlgorithm1D(strain, pls_n0, a_n0, b_n0)
				stress = output['stress']; pls_n1 = output['pls']; a_n1 = output['alpha']
				b_n1 = output['beta']; Cep = output['mechArgs']

				# Compute internal force 
				Fint_dj = self.compute_MechStaticIntForce(stress)
				
				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				resNRj = np.sqrt(np.dot(r_dj, r_dj))
				if j == 0: resNR0 = resNRj
				print('NR error %.5e' %resNRj)
				if resNRj <= max([self._thresholdNewton*resNR0, 1e-14]): break
				if j > 0 and isElasticLoad: break

				# Solver for active control points
				tangentM = sp.csr_matrix(self.compute_tangentMatrix(Cep)[np.ix_(dof, dof)])
				deltaD = np.zeros(self.nbctrlpts); deltaD[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update active control points
				dj_n1 += deltaD
				if np.sqrt(np.dot(deltaD, deltaD)) <= 1e-12: break

			dispinout[:, i] = dj_n1
			Allstrain[:, i] = np.ravel(strain)
			Allstress[:, i] = np.ravel(stress)
			Allplseq[:, i]  = a_n1
			AllCep[:, i]    = np.ravel(Cep)

			pls_n0, a_n0, b_n0 = np.copy(pls_n1), np.copy(a_n1), np.copy(b_n1)

		return Allstrain, Allstress, Allplseq, AllCep

