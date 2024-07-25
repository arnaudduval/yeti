"""
.. This module contains functions to treat transient heat problems in 1D. 
.. It also contains functions to treat elastoplasticity problems in 1D. 
.. Joaquin Cornejo
"""

from . import *
from .lib_part import part1D
from .lib_quadrules import GaussQuadrature
from .lib_material import heatmat, mechamat
from .lib_boundary import boundaryCondition
from .lib_base import evalDersBasisCSRFortran, array2csr_matrix, macaulayfunc

class problem1D():
	def __init__(self, part:part1D, boundary:boundaryCondition, solverArgs={}):
		self.part = part
		self.boundary = boundary
		self.addSolverConstraints(solverArgs=solverArgs)
		return
	
	def addSolverConstraints(self, solverArgs:dict):
		self._thresNL = solverArgs.get('thres_nonlinear', 1e-8)
		self._itersNL = solverArgs.get('iters_nonlinear', 20)
		self._safeguard = 1e-14
		return

	def compute_volForce(self, volfun, args=None):
		" Computes 'volumetric' source vector in 1D "
		if args is None: args={'position':self.part.qpPhy}
		prop  = volfun(args)*self.part.detJ
		force = self.part._denseweights[0] @ prop
		return force
		
	def interpolateMeshgridField(self, u_ctrlpts, sampleSize=101, isSample=True):
		if isSample: basis = self.part.quadRule.getSampleBasis(sampleSize=sampleSize)[0]
		else: 		 basis = self.part.quadRule.getDenseQuadRules()[0]
		u_interp = basis[0].T @ u_ctrlpts
		if isSample: x_interp = basis[0].T @ self.part.ctrlpts
		else: 		 x_interp = self.part.qpPhy
		return u_interp, x_interp
	
	def L2projectionCtrlpts(self, u_atqp):
		" Interpolate control point field (from parametric to physical space) "
		force = self.part._denseweights[0] @ np.diag(self.part.detJ) @ u_atqp
		masse = self.part._denseweights[0] @ np.diag(self.part.detJ) @ self.part._densebasis[0].T
		u_ctrlpts = sp.linalg.spsolve(sp.csr_matrix(masse), force)
		return u_ctrlpts
	
	def normOfError(self, u_ctrlpts, normArgs:dict):
		""" Computes the norm L2 or H1 of the error. The exactfun is the function of the exact solution. 
			and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
			whether the default quadrature is weighted quadrature. 
		"""	
		typeNorm = normArgs.get('type', 'l2').lower()
		if all(norm != typeNorm  for norm in ['l2', 'h1', 'semih1']): raise Warning('Unknown norm')

		# Compute u interpolation
		quadRule = GaussQuadrature(self.part.degree, self.part.knotvector, quadArgs={'type':'leg'})
		quadPts  = quadRule.getQuadratureRulesInfo()[0]
		denseBasis, parametricWeights = quadRule._denseBasis, quadRule._parametricWeights
		
		qpPhy = denseBasis[0].T @ self.part.ctrlpts 
		detJ  = denseBasis[1].T @ self.part.ctrlpts
		u_interp = denseBasis[0].T @ u_ctrlpts
		uders_interp = denseBasis[1].T @ u_ctrlpts / detJ

		# Compute u exact
		u_exact, uders_exact = None, None
		exactfun = normArgs.get('exactFunction', None)
		exactfunders = normArgs.get('exactFunctionDers', None)
		exactextraArgs = normArgs.get('exactExtraArgs', None)
		if exactextraArgs is not None:
			assert isinstance(exactextraArgs, dict), 'Error type of extra args'
			if not 'position' in exactextraArgs.keys(): exactextraArgs['position'] = qpPhy
			if callable(exactfun): u_exact = exactfun(exactextraArgs)
			if callable(exactfunders): uders_exact = exactfunders(exactextraArgs)
		else:
			if callable(exactfun): u_exact = exactfun(qpPhy)
			if callable(exactfunders): uders_exact = exactfunders(qpPhy)

		part_ref = normArgs.get('part_ref', None); u_ref = normArgs.get('u_ref', None)
		if isinstance(part_ref, problem1D) and isinstance(u_ref, np.ndarray):
			denseBasisExact = []
			basis_csr, indi_csr, indj_csr = evalDersBasisCSRFortran(part_ref.part.degree, part_ref.part.knotvector, quadPts)
			for i in range(2): denseBasisExact.append(array2csr_matrix(basis_csr[:, i], indi_csr, indj_csr))
			u_exact = denseBasisExact[0].T @ u_ref
			invJExact = denseBasisExact[1].T @ part_ref.part.ctrlpts
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
		
		abserror = np.sqrt(tmp1)
		relerror = np.sqrt(tmp1/tmp2) if tmp2!=0 else np.sqrt(tmp1)
		if tmp2 == 0: print('Warning: Dividing by zero')

		return abserror, relerror

class heatproblem1D(problem1D):
	def __init__(self, heat_material:heatmat, part:part1D, boundary:boundaryCondition, solverArgs={}):
		problem1D.__init__(self, part, boundary, solverArgs=solverArgs)	
		self.heatmaterial = heat_material
		if self.heatmaterial.density is None: self.heatmaterial.density = lambda x: np.ones(self.part.nbqp)
		return
	
	def compute_mfCapacity(self, array_in, args=None, isLumped=False):
		if args is None: args = {'position': self.part.qpPhy}
		prop = self.heatmaterial.capacity(args)*self.heatmaterial.density(args)*self.part.detJ
		matrix = self.part._denseweights[0] @ np.diag(prop) @ self.part._densebasis[0].T
		if isLumped: matrix = np.diag(matrix.sum(axis=1))
		array_out = matrix @ array_in
		return array_out
	
	def compute_mfConductivity(self, array_in, args=None):
		if args is None: args = {'position': self.part.qpPhy}
		prop = self.heatmaterial.conductivity(args)*self.part.invJ
		matrix = self.part._denseweights[-1] @ np.diag(prop) @ self.part._densebasis[1].T 
		array_out = matrix @ array_in
		return array_out

	def interpolate_temperature(self, u_ctrlpts):
		" Interpolate temperature in 1D "
		u_interp = self.part._densebasis[0].T @ u_ctrlpts
		return u_interp

	def compute_FourierTangentMatrix(self, dt, alpha=0.5, args=None, isLumped=False):
		""" Computes tangent matrix in transient heat
			S = C + theta dt K
			K = int_Omega dB/dx Kprop dB/dx dx = int_[0, 1] J^-1 dB/dxi Kprop J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
			C = int_Omega B Cprop B dx = int [0, 1] B Cprop det J B dxi
		"""
		if args is None: args = {'position': self.part.qpPhy}
		Kprop = self.heatmaterial.conductivity(args)*self.part.invJ
		K = self.part._denseweights[-1] @ np.diag(Kprop) @ self.part._densebasis[1].T 
		Cprop = self.heatmaterial.capacity(args)*self.heatmaterial.density(args)*self.part.detJ
		C = self.part._denseweights[0] @ np.diag(Cprop) @ self.part._densebasis[0].T
		if isLumped: C = np.diag(C.sum(axis=1))
		matrix = C + dt*alpha*K
		return matrix
	
	def compute_CattaneoMatrix(self, dt, beta=0.25, gamma=0.5, args=None, isLumped=False):
		if args is None: args = {'position': self.part.qpPhy}
		Kprop = self.heatmaterial.conductivity(args)*self.part.invJ
		K = self.part._denseweights[-1] @ np.diag(Kprop) @ self.part._densebasis[1].T 
		Cprop = self.heatmaterial.capacity(args)*self.heatmaterial.density(args)*self.part.detJ
		C = self.part._denseweights[0] @ np.diag(Cprop) @ self.part._densebasis[0].T
		Mprop = self.heatmaterial.relaxation(args)*Cprop
		M = self.part._denseweights[0] @ np.diag(Mprop) @ self.part._densebasis[0].T
		if isLumped: M = np.diag(M.sum(axis=1))
		matrix = M + gamma*dt*C + beta*dt**2*K
		return matrix

	def solveFourierTransientProblem(self, Tinout, Fext_list, time_list, alpha=0.5, isLumped=False, extraArgs=None):
		" Solves transient heat problem in 1D. "

		if extraArgs is None: extraArgs = {'stabilized': False, 'forces': None}
		useStabilization = extraArgs.get('stabilized', False); forces = extraArgs.get('forces', None)
		if useStabilization and forces is None: raise Warning('Do not forget to set volumetric force if stabilization is needed')
		def computeVelocity(problem:heatproblem1D, Fext, args=None, isLumped=False):
			if args is None: args = {'position': problem.part.qpPhy}
			dof = problem.boundary.thdof
			prop = problem.heatmaterial.capacity(args)*problem.heatmaterial.density(args)*problem.part.detJ
			matrix = problem.part._denseweights[0] @ np.diag(prop) @ problem.part._densebasis[0].T
			if isLumped: matrix = np.diag(matrix.sum(axis=1))
			tmp = np.linalg.solve(matrix[np.ix_(dof, dof)], Fext[dof])
			velocity = np.zeros(shape=np.shape(Fext)); velocity[dof] = tmp
			return velocity

		nbctrlpts_total = self.part.nbctrlpts_total; nsteps = len(time_list)
		dod = self.boundary.thdod; dof = self.boundary.thdof

		# Compute initial velocity using interpolation
		assert nsteps > 2, 'At least 2 steps'
		V_n0 = np.zeros(nbctrlpts_total)
		dt = time_list[1] - time_list[0]
		V_n0[dod] = 1.0/dt*(Tinout[dod, 1] - Tinout[dod, 0])

		temperature = self.interpolate_temperature(Tinout[:, 0])
		args ={'temperature':temperature, 'position':self.part.qpPhy}
		tmp = (Fext_list[:, 0] - self.compute_mfCapacity(V_n0, args=args, isLumped=isLumped)
				- self.compute_mfConductivity(Tinout[:, 0], args=args))
		V_n0[dof] = computeVelocity(self, tmp, args=args, isLumped=isLumped)[dof]

		for i in range(1, nsteps):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			d_n0 = np.copy(Tinout[:, i-1])

			# Predict values of new step
			dj_n1 = d_n0 + (1 - alpha)*dt*V_n0
			Vj_n1 = np.zeros(nbctrlpts_total)
			
			# Overwrite inactive control points
			Vj_n1[dod] = 1.0/(alpha*dt)*(Tinout[dod, i] - dj_n1[dod])
			dj_n1[dod] = Tinout[dod, i]

			Fext_n1 = np.copy(Fext_list[:, i])

			print('Step: %d' %i)
			for j in range(self._itersNL): 
				
				# Interpolate temperature
				temperature = self.interpolate_temperature(dj_n1)

				# Compute internal force
				args = {'position':self.part.qpPhy, 'temperature':temperature}
				Fint_dj = self.compute_mfCapacity(Vj_n1, args=args, isLumped=isLumped) + self.compute_mfConductivity(dj_n1, args=args)

				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				resNLj = np.sqrt(np.dot(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('NonLinear error: %.5e' %resNLj)
				if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break
				
				# Solver for active control points
				tangentM = sp.csr_matrix(self.compute_FourierTangentMatrix(dt, alpha=alpha, args=args, isLumped=isLumped)[np.ix_(dof, dof)])
				deltaV = np.zeros(nbctrlpts_total); deltaV[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update active control points
				Vj_n1 += deltaV
				dj_n1 += alpha*dt*deltaV

			Tinout[:, i] = np.copy(dj_n1)
			V_n0 = np.copy(Vj_n1)

		return 

	def solveCattaneoTransientProblem(self, Tinout, Fext_list, time_list, beta=0.25, gamma=0.5, isLumped=False):
		" Solves transient heat problem in 1D using Cattaneo approach. "

		nbctrlpts_total = self.part.nbctrlpts_total; nsteps = len(time_list)
		dod  = self.boundary.thdod; dof = self.boundary.thdof

		# Compute initial velocity using interpolation
		assert nsteps > 2, 'At least 2 steps'
		V_n0 = np.zeros(nbctrlpts_total); A_n0 = np.zeros(nbctrlpts_total)
		dt1 = time_list[1] - time_list[0]
		dt2 = time_list[2] - time_list[0]
		factor = dt2/dt1
		V_n0[dod] = 1.0/(dt1*(factor-factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
		A_n0[dod] = 2.0/(dt1*dt2)*((Tinout[dod, 2] - factor*Tinout[dod, 1])/(factor - 1) + Tinout[dod, 0])
		# TO DO: ADD A PROCESS TO COMPUTE THE VELOCITY AND ACCELERATION IN DOF

		for i in range(1, nsteps):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]

			# Get values of last step
			d_n0 = np.copy(Tinout[:, i-1])

			# Predict values of new step
			dj_n1 = d_n0 + dt*V_n0 + 0.5*dt**2*(1 - 2*beta)*A_n0
			Vj_n1 = V_n0 + (1 - gamma)*dt*A_n0
			Aj_n1 = np.zeros(nbctrlpts_total) 

			# Overwrite inactive control points
			Aj_n1[dod] = (Tinout[dod, i] - dj_n1[dod])/(beta*dt**2)
			Vj_n1[dod] = Vj_n1[dod] + gamma*dt*Aj_n1[dod]
			dj_n1[dod] = Tinout[dod, i]

			Fext_n1 = np.copy(Fext_list[:, i])

			print('Step: %d' %i)
			for j in range(self._itersNL): 
				
				# Interpolate temperature
				temperature = self.interpolate_temperature(dj_n1)

				# Compute internal force
				args = {'position':self.part.qpPhy, 'temperature':temperature}
				Fint_dj = (self.compute_mfCapacity(Aj_n1, args=args, isLumped=isLumped) 
						+ self.compute_mfCapacity(Vj_n1, args=args, isLumped=False) 
						+ self.compute_mfConductivity(dj_n1, args=args)
				)

				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				resNLj = np.sqrt(np.dot(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('NonLinear error: %.5e' %resNLj)
				if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break

				# Solve for active control points
				tangentM = sp.csr_matrix(self.compute_CattaneoMatrix(dt=dt, beta=beta, gamma=gamma, 
												args=args, isLumped=isLumped)[np.ix_(dof, dof)])
				deltaA = np.zeros(nbctrlpts_total); deltaA[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update active control points
				dj_n1 += beta*dt**2*deltaA
				Vj_n1 += gamma*dt*deltaA
				Aj_n1 += deltaA

			Tinout[:, i] = np.copy(dj_n1)
			V_n0 = np.copy(Vj_n1)
			A_n0 = np.copy(Aj_n1)

		return

class mechaproblem1D(problem1D):
	def __init__(self, mechanical_material:mechamat, part:part1D, boundary:boundaryCondition, solverArgs={}):
		problem1D.__init__(self, part, boundary, solverArgs=solverArgs)
		self.mechamaterial = mechanical_material		
		if self.mechamaterial.density is None: self.mechamaterial.density = lambda x: np.ones(self.part.nbqp)
		return

	def interpolate_strain(self, disp):
		" Computes strain field from a given displacement field "
		strain = self.part._densebasis[1].T @ disp * self.part.invJ
		strain = np.reshape(strain, (1, len(strain)))
		return strain

	def compute_MechStaticIntForce(self, stress):
		""" Computes internal force Fint. 
			Fint = int_Omega dB/dx sigma dx = int_[0, 1] J^-1 dB/dxi sigma detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		Fint = self.part._denseweights[2] @ np.ravel(stress).T
		return Fint

	def compute_PlasticityTangentMatrix(self, Cep):
		""" Computes stiffness matrix in elastoplasticity
			S = int_Omega dB/dx Dalg dB/dx dx = int_[0, 1] J^-1 dB/dxi Dalg J^-1 dB/dxi detJ dxi.
			But in 1D: detJ times J^-1 get cancelled.
		"""
		prop = np.ravel(Cep)*self.part.invJ
		matrix = self.part._denseweights[-1] @ np.diag(prop) @ self.part._densebasis[1].T 
		return matrix

	def solvePlasticityProblem(self, dispinout, Fext_list):
		" Solves elasto-plasticity problem in 1D. It considers Dirichlet boundaries equal to 0 "

		nbChaboche = self.mechamaterial._chabocheNBparameters
		nbqp = self.part.nbqp; dof = self.boundary.thdof; dod = self.boundary.thdod
		plasticstrain_n0, plseq_n0, back_n0 = np.zeros((1, nbqp)), np.zeros((1, nbqp)), np.zeros((nbChaboche, 1, nbqp))
		plasticstrain_n1, plseq_n1, back_n1 = np.zeros((1, nbqp)), np.zeros((1, nbqp)), np.zeros((nbChaboche, 1, nbqp))

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
			for j in range(self._itersNL): # Newton-Raphson 
				
				# Compute strain at each quadrature point
				strain = self.interpolate_strain(dj_n1)

				# Find closest point projection 
				output, isElasticLoad = self.mechamaterial.J2returnMappingAlgorithm1D(strain, plasticstrain_n0, plseq_n0, back_n0)
				stress = output['stress']; plasticstrain_n1 = output['plastic']; plseq_n1 = output['plseq']
				back_n1 = output['back']; Cep = output['mechArgs']

				# Compute internal force 
				Fint_dj = self.compute_MechStaticIntForce(stress)
				
				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				resNLj = np.sqrt(np.dot(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('NonLinear error: %.5e' %resNLj)
				if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break
				if j > 0 and isElasticLoad: break

				# Solver for active control points
				tangentM = sp.csr_matrix(self.compute_PlasticityTangentMatrix(Cep)[np.ix_(dof, dof)])
				deltaD = np.zeros(self.part.nbctrlpts_total); deltaD[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update active control points
				dj_n1 += deltaD

			dispinout[:, i] = dj_n1
			Allstrain[:, i] = np.ravel(strain)
			Allstress[:, i] = np.ravel(stress)
			Allplseq[:, i]  = plseq_n1
			AllCep[:, i]    = np.ravel(Cep)

			plasticstrain_n0, plseq_n0, back_n0 = np.copy(plasticstrain_n1), np.copy(plseq_n1), np.copy(back_n1)

		return Allstrain, Allstress, Allplseq, AllCep

	def solveViscoPlasticityProblem(self, dispinout, Fext_list, time_list):

		""" Solves elasto-viscoplasticity problem in 1D. 
			For the moment, it only works with isotropic hardening without kinematic hardening. 
		"""

		assert self.mechamaterial._chabocheNBparameters <= 1, 'Try another method'
		assert not np.any(self.mechamaterial._chabocheTable), 'Try another method'
		assert self.mechamaterial._isoHardening._isoname == 'linear', 'Try another method'
		nbChaboche = self.mechamaterial._chabocheNBparameters
		nbqp = self.part.nbqp; dof = self.boundary.thdof; dod = self.boundary.thdod
		plasticstrain_n0, plseq_n0 = np.zeros((1, nbqp)), np.zeros((1, nbqp))
		plasticstrain_n1, plseq_n1 = np.zeros((1, nbqp)), np.zeros((1, nbqp))
		back_tmp = np.zeros((nbChaboche, 1, nbqp))

		Allstrain  = np.zeros((nbqp, np.shape(Fext_list)[1]))
		Allplseq = np.zeros((nbqp, np.shape(Fext_list)[1]))
		Allstress  = np.zeros((nbqp, np.shape(Fext_list)[1]))
		AllCep = np.zeros((nbqp, np.shape(Fext_list)[1]))

		for i in range(1, np.shape(Fext_list)[1]):

			# Get time step
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			d_n0 = np.copy(dispinout[:, i-1])
			
			# Predict values of new step
			dj_n1 = np.copy(d_n0)

			# Overwrite inactive control points
			dj_n1[dod] = np.copy(dispinout[dod, i])

			Fext_n1 = np.copy(Fext_list[:, i])

			print('Step: %d' %i)
			for j in range(self._itersNL): # Newton-Raphson 
				
				# Compute strain at each quadrature point
				strain = self.interpolate_strain(dj_n1)

				# Find closest point projection 
				output, isElasticLoad = self.mechamaterial.J2returnMappingAlgorithm1D(strain, plasticstrain_n0, plseq_n0, back_tmp)
				stress_inf = output['stress']; plasticstrain_inf = output['plastic']; plseq_inf = output['plseq']
				Cep_inf = output['mechArgs']; plsInd = output['plsInd']

				stress = np.copy(stress_inf)
				plseq_n1 = np.copy(plseq_inf)
				plasticstrain_n1 = np.copy(plasticstrain_inf)
				Cep = np.copy(Cep_inf)

				# Perform viscolpastic regularization
				if np.size(plsInd) > 0:
					eta = self.mechamaterial.viscoparameter
					E = self.mechamaterial.elasticModulus
					H = self.mechamaterial._isoHardening._isohardfun(plseq_n1[:, plsInd])
					tau = eta/(E + H)

					stress[:, plsInd] = (E*(strain[:, plsInd] - plasticstrain_n0[:, plsInd]) 
										+ dt/tau*stress_inf[:, plsInd])/(1 + dt/tau)
					plseq_n1[:, plsInd] = (plseq_n0[:, plsInd] + dt/tau*plseq_inf[:, plsInd])/(1 + dt/tau)
					plasticstrain_n1[:, plsInd] = strain[:, plsInd] - 1/E*stress[:, plsInd]
					Cep[0, plsInd] = (E + dt/tau*Cep_inf[0, plsInd])/(1 + dt/tau)

				# Compute internal force 
				Fint_dj = self.compute_MechStaticIntForce(stress)
				
				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				resNLj = np.sqrt(np.dot(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('NonLinear error: %.5e' %resNLj)
				if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break
				if j > 0 and isElasticLoad: break

				# Solver for active control points
				tangentM = sp.csr_matrix(self.compute_PlasticityTangentMatrix(Cep)[np.ix_(dof, dof)])
				deltaD = np.zeros(self.part.nbctrlpts_total); deltaD[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update active control points
				dj_n1 += deltaD

			dispinout[:, i] = dj_n1
			Allstrain[:, i] = np.ravel(strain)
			Allstress[:, i] = np.ravel(stress)
			Allplseq[:, i]  = plseq_n1
			AllCep[:, i]    = np.ravel(Cep)

			plasticstrain_n0, plseq_n0 = np.copy(plasticstrain_n1), np.copy(plseq_n1)
		return Allstrain, Allstress, Allplseq, AllCep

class stproblem1D():

	def __init__(self, part:part1D, tspan:part1D, boundary:boundaryCondition, solverArgs={}):
		self.part = part
		self.time = tspan
		self.boundary = boundary
		self.addSolverConstraints(solverArgs=solverArgs)
		return
	
	def addSolverConstraints(self, solverArgs:dict):
		self._thresNL = solverArgs.get('thres_nonlinear', 1e-8)
		self._itersNL = solverArgs.get('iters_nonlinear', 20)
		self._safeguard = 1e-14
		return

	def compute_volForce(self, volfun, args=None):
		" Computes 'volumetric' source vector in 1D "
		if args is None: args={'position':self.part.qpPhy}
		prop  = np.kron(self.time.detJ, self.part.detJ)*volfun(args)
		force = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ prop
		return force
	
	def normOfError(self, u_ctrlpts, normArgs:dict):
		""" Computes the norm L2 or H1 of the error. The exactfun is the function of the exact solution. 
			and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
			whether the default quadrature is weighted quadrature. 
		"""	
		typeNorm = normArgs.get('type', 'l2').lower()
		if all(norm != typeNorm  for norm in ['l2', 'h1', 'semih1']): raise Warning('Unknown norm')

		# Compute u interpolation
		quadRulePart = GaussQuadrature(self.part.degree, self.part.knotvector, quadArgs={'type':'leg'})
		quadRulePart.getQuadratureRulesInfo()

		quadRuleTime = GaussQuadrature(self.time.degree, self.time.knotvector, quadArgs={'type':'leg'})
		quadRuleTime.getQuadratureRulesInfo()
		parametricWeights = np.kron(quadRuleTime._parametricWeights, quadRulePart._parametricWeights)

		sptimectrlpts = np.zeros((2, self.part.nbctrlpts_total*self.time.nbctrlpts_total))
		for i in range(self.time.nbctrlpts_total):
			iold = i*self.part.nbctrlpts_total; inew = (i + 1)*self.part.nbctrlpts_total
			sptimectrlpts[:, iold:inew] = np.vstack([self.part.ctrlpts, self.time.ctrlpts[i]*np.ones(self.part.nbctrlpts_total)])
		
		qpPhy = np.zeros((2, quadRulePart.nbqp*quadRuleTime.nbqp))
		for i in range(2): qpPhy[i, :] = sp.kron(quadRuleTime._denseBasis[0], quadRulePart._denseBasis[0]).T @ sptimectrlpts[i, :]

		Jqp = np.zeros((2, 2, quadRulePart.nbqp*quadRuleTime.nbqp))
		for j in range(2):
			for i in range(2):
				beta = np.zeros(2, dtype=int); beta[j] = 1
				Jqp[i, j, :] = sp.kron(quadRuleTime._denseBasis[beta[1]], quadRulePart._denseBasis[beta[0]]).T @ sptimectrlpts[i, :]

		detJ = Jqp[0, 0, :]*Jqp[1, 1, :] - Jqp[0, 1, :]*Jqp[1, 0, :]
		invJ = np.zeros((2, 2, quadRulePart.nbqp*quadRuleTime.nbqp))
		invJ[0, 0, :] = Jqp[1, 1, :]/detJ; invJ[1, 1, :] = Jqp[0, 0, :]/detJ
		invJ[0, 1, :] = -Jqp[0, 1, :]/detJ; invJ[1, 0, :] = -Jqp[1, 0, :]/detJ 

		u_interp = sp.kron(quadRuleTime._denseBasis[0], quadRulePart._denseBasis[0]).T @ u_ctrlpts
		derstemp = np.zeros((1, 2, quadRulePart.nbqp*quadRuleTime.nbqp))
		derstemp[0, 0, :] = sp.kron(quadRuleTime._denseBasis[0], quadRulePart._denseBasis[1]).T @ u_ctrlpts
		derstemp[0, 1, :] = sp.kron(quadRuleTime._denseBasis[1], quadRulePart._denseBasis[0]).T @ u_ctrlpts
		uders_interp = np.einsum('ijl,jkl->ikl', derstemp, invJ); uders_interp = uders_interp[0, :, :]

		# Compute u exact
		u_exact, uders_exact = None, None
		exactfun = normArgs.get('exactFunction', None)
		exactfunders = normArgs.get('exactFunctionDers', None)
		exactextraArgs = normArgs.get('exactExtraArgs', None)
		if exactextraArgs is not None:
			assert isinstance(exactextraArgs, dict), 'Error type of extra args'
			if not 'position' in exactextraArgs.keys(): exactextraArgs['position'] = qpPhy
			if callable(exactfun): u_exact = exactfun(exactextraArgs)
			if callable(exactfunders): uders_exact = exactfunders(exactextraArgs)
		else:
			if callable(exactfun): u_exact = exactfun(qpPhy)
			if callable(exactfunders): uders_exact = exactfunders(qpPhy)

		part_ref = normArgs.get('part_ref', None); time_ref = normArgs.get('time_ref', None); u_ref = normArgs.get('u_ref', None)
		if isinstance(part_ref, part1D) and isinstance(time_ref, part1D) and isinstance(u_ref, np.ndarray):
			# Compute u exact
			quadRulePart = GaussQuadrature(part_ref.degree, part_ref.knotvector, quadArgs={'type':'leg'})
			quadRulePart.getQuadratureRulesInfo()

			quadRuleTime = GaussQuadrature(time_ref.degree, time_ref.knotvector, quadArgs={'type':'leg'})
			quadRuleTime.getQuadratureRulesInfo()

			sptimectrlpts = np.zeros((2, part_ref.nbctrlpts_total*time_ref.nbctrlpts_total))
			for i in range(time_ref.nbctrlpts_total):
				iold = i*part_ref.nbctrlpts_total; inew = (i + 1)*part_ref.nbctrlpts_total
				sptimectrlpts[:, iold:inew] = np.vstack([part_ref.ctrlpts, time_ref.ctrlpts[i]*np.ones(part_ref.nbctrlpts_total)])
			
			qpPhyExact = np.zeros((2, quadRulePart.nbqp*quadRuleTime.nbqp))
			for i in range(2): qpPhyExact[i, :] = sp.kron(quadRuleTime._denseBasis[0], quadRulePart._denseBasis[0]).T @ sptimectrlpts

			JqpExact = np.zeros((2, 2, quadRulePart.nbqp*quadRuleTime.nbqp))
			for j in range(2):
				for i in range(2):
					beta = np.zeros(2, dtype=int); beta[j] = 1
					JqpExact[i, j, :] = sp.kron(quadRuleTime._denseBasis[beta[1]], quadRulePart._denseBasis[beta[0]]).T @ sptimectrlpts[i]

			detJExact = JqpExact[0, 0, :]*JqpExact[1, 1, :] - JqpExact[0, 1, :]*JqpExact[1, 0, :]
			invJExact = np.zeros((2, 2, quadRulePart.nbqp*quadRuleTime.nbqp))
			invJExact[0, 0, :] = JqpExact[1, 1, :]/detJExact; invJExact[1, 1, :] = JqpExact[0, 0, :]/detJExact
			invJExact[0, 1, :] = -JqpExact[0, 1, :]/detJExact; invJExact[1, 0, :] = -JqpExact[1, 0, :]/detJExact 

			u_exact = sp.kron(quadRuleTime._denseBasis[0], quadRulePart._denseBasis[0]).T @ u_ctrlpts
			derstemp = np.zeros((1, 2, quadRulePart.nbqp*quadRuleTime.nbqp))
			derstemp[0, 0, :] = sp.kron(quadRuleTime._denseBasis[0], quadRulePart._denseBasis[1]).T @ u_ctrlpts
			derstemp[0, 1, :] = sp.kron(quadRuleTime._denseBasis[1], quadRulePart._denseBasis[0]).T @ u_ctrlpts
			uders_exact = np.einsum('ijl,jkl->ikl', derstemp, invJExact); uders_exact = uders_exact[0, :, :]
			
		# Compute error
		uedfuf2_l2, uedfuf2_sh1 = 0., 0.
		ue2_l2, ue2_sh1 = 0., 0.

		if typeNorm == 'l2' or typeNorm == 'h1':
			uedfuf2_l2 += (u_exact - u_interp)**2
			ue2_l2     += u_exact**2

		if typeNorm == 'h1' or typeNorm == 'semih1':
			uedfuf2_sh1 += np.einsum('il->l', (uders_exact - uders_interp)**2)
			ue2_sh1     += np.einsum('il->l', uders_exact**2)

		norm1 = (uedfuf2_l2 + uedfuf2_sh1)*detJ
		norm2 = (ue2_l2 + ue2_sh1)*detJ

		tmp1 = np.einsum('i,i->', parametricWeights, norm1)
		tmp2 = np.einsum('i,i->', parametricWeights, norm2)
		
		abserror = np.sqrt(tmp1)
		relerror = np.sqrt(tmp1/tmp2) if tmp2!=0 else np.sqrt(tmp1)
		if tmp2 == 0: print('Warning: Dividing by zero')

		return abserror, relerror
	
class stheatproblem1D(stproblem1D):
	def __init__(self, heat_material:heatmat, part:part1D, tspan:part1D, boundary:boundaryCondition, solverArgs={}):
		stproblem1D.__init__(self, part, tspan, boundary, solverArgs)
		self.heatmaterial = heat_material
		if self.heatmaterial.density is None: self.heatmaterial.density = lambda x: np.ones(self.part.nbqp_total)
		return
	
	def compute_mfSTConductivity(self, array_in, args=None):
		if args is None: 
			temperature = np.ones(self.part.nbqp*self.time.nbqp) 
			if self.heatmaterial._isConductivityIsotropic: args = temperature
			else: args = {'temperature': temperature}
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T @ array_in
		tmp2 = (self.heatmaterial.conductivity(args)*np.kron(self.time.detJ, self.part.invJ))*tmp1
		array_out = sp.kron(self.time._denseweights[0], self.part._denseweights[-1]) @ tmp2
		return array_out
		
	def compute_mfSTCapacity(self, array_in, args=None):
		if args is None: 
			temperature = np.ones(self.part.nbqp_total*self.time.nbqp) 
			if self.heatmaterial._isCapacityIsotropic: args = temperature
			else: args = {'temperature': temperature}
		tmp1 = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T @ array_in
		tmp2 = self.heatmaterial.capacity(args)*np.kron(np.ones(self.time.nbqp), self.part.detJ)*tmp1
		array_out = sp.kron(self.time._denseweights[1], self.part._denseweights[0]) @ tmp2
		return array_out

	def interpolate_STtemperature_gradients(self, u_ctrlpts):
		sptimectrlpts = np.zeros((2, self.part.nbctrlpts_total*self.time.nbctrlpts_total))
		for i in range(self.time.nbctrlpts_total):
			iold = i*self.part.nbctrlpts_total; inew = (i + 1)*self.part.nbctrlpts_total
			sptimectrlpts[:, iold:inew] = np.vstack([self.part.ctrlpts, self.time.ctrlpts[i]*np.ones(self.part.nbctrlpts_total)])
		
		Jqp = np.zeros((2, 2, self.part.nbqp*self.time.nbqp))
		for j in range(2):
			for i in range(2):
				beta = np.zeros(2, dtype=int); beta[j] = 1
				Jqp[i, j, :] = sp.kron(self.time._densebasis[beta[1]], self.part._densebasis[beta[0]]).T @ sptimectrlpts[i]

		detJ = Jqp[0, 0, :]*Jqp[1, 1, :] - Jqp[0, 1, :]*Jqp[1, 0, :]
		invJ = np.zeros((2, 2, self.part.nbqp*self.time.nbqp))
		invJ[0, 0, :] = Jqp[1, 1, :]/detJ; invJ[1, 1, :] = Jqp[0, 0, :]/detJ
		invJ[0, 1, :] = -Jqp[0, 1, :]/detJ; invJ[1, 0, :] = -Jqp[1, 0, :]/detJ 

		u_interp = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T @ u_ctrlpts
		derstemp = np.zeros((1, 2, self.part.nbqp*self.time.nbqp))
		derstemp[0, 0, :] = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T @ u_ctrlpts
		derstemp[0, 1, :] = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T @ u_ctrlpts
		uders_interp = np.einsum('ijl,jkl->ikl', derstemp, invJ); uders_interp = uders_interp[0, :, :]
		return u_interp, uders_interp
	
	def compute_TangentMatrix(self, args=None, isfull=False):
		assert args is not None, 'Please enter a valid argument'
		if args is None: temperature = np.ones(self.part.nbqp*self.time.nbqp) 
		
		if args is None: 
			if self.heatmaterial._isConductivityIsotropic: args = temperature
			else: args = {'temperature': temperature}
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T
		tmp2 = sp.diags(self.heatmaterial.conductivity(args)*np.kron(self.time.detJ, self.part.invJ)) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[-1]) @ tmp2

		if args is None: 
			if self.heatmaterial._isCapacityIsotropic: args = temperature
			else: args = {'temperature': temperature}
		tmp1 = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T
		tmp2 = sp.diags(self.heatmaterial.capacity(args)*np.kron(np.ones(self.time.nbqp), self.part.detJ)) @ tmp1
		matrix += sp.kron(self.time._denseweights[1], self.part._denseweights[0]) @ tmp2

		if not isfull: return matrix

		if args is None: 
			if self.heatmaterial._isConductivityIsotropic: args = temperature
			else: args = {'temperature': temperature}
			gradTemperature = np.ones((2, self.part.nbqp_total*self.time.nbqp))
		else: gradTemperature = args['gradients']
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(self.heatmaterial.conductivityDers(args)*gradTemperature[0, :]*np.kron(self.time.detJ, np.ones(self.part.nbqp))) @ tmp1
		matrix += sp.kron(self.time._denseweights[0], self.part._denseweights[2]) @ tmp2

		if args is None: 
			if self.heatmaterial._isCapacityIsotropic: args = temperature
			else: args = {'temperature': temperature}
			gradTemperature = np.ones((2, self.part.nbqp_total*self.time.nbqp))
		else: gradTemperature = args['gradients']
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags((self.heatmaterial.capacityDers(args)*gradTemperature[-1, :]*np.kron(self.time.detJ, self.part.detJ))) @ tmp1
		matrix += sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix
	
	def solveFourierSTHeatProblem(self, Tinout, Fext, isfull=False, isadaptive=False):

		dod = self.boundary.thdod; dof = self.boundary.thdof
		nbctrlpts_total = self.part.nbctrlpts_total*self.time.nbctrlpts_total
		AllresNewton, Allsol, Alldelta = [], [], []
		for j in range(self._itersNL):

			# Compute temperature at each quadrature point
			temperature, gradtemperature = self.interpolate_STtemperature_gradients(Tinout)

			# Compute internal force
			Fint_dj = (self.compute_mfSTCapacity(Tinout, args={'temperature':temperature}) 
						+ self.compute_mfSTConductivity(Tinout, args={'temperature':temperature}))

			# Compute residue
			r_dj = Fext - Fint_dj
			r_dj[dod] = 0.0

			resNLj1 = np.sqrt(np.dot(r_dj, r_dj))
			print('Nonlinear error: %.3e' %resNLj1)

			if j == 0: resNL0 = np.copy(resNLj1)
			if resNLj1 <= max([self._safeguard, self._thresNL*resNL0]): break
			AllresNewton.append(resNLj1); Allsol.append(np.copy(Tinout))

			# Solve for active control points
			tangentM = self.compute_TangentMatrix(
										args={'temperature':temperature, 
											'gradients':gradtemperature}, 
										isfull=isfull).todense()
			tangentM = sp.csr_matrix(tangentM[np.ix_(dof, dof)])
			deltaD = np.zeros(nbctrlpts_total)
			# deltaD[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])
			deltaD[dof] = np.linalg.solve(tangentM.todense(), r_dj[dof])
		
			# Update active control points
			Tinout += deltaD
			Alldelta.append(deltaD)
			
		output = {'NewtonRes':AllresNewton, 'Solution':Allsol, 'Delta':Alldelta}
		return output

class stmechaproblem1D(stproblem1D):

	def __init__(self, mechanical_material:mechamat, part:part1D, tspan:part1D, boundary:boundaryCondition, solverArgs={}):
		stproblem1D.__init__(self, part, tspan, boundary, solverArgs)
		self.mechamaterial = mechanical_material
		if self.mechamaterial.density is None: self.mechamaterial.density = lambda x: np.ones(self.part.nbqp_total)
		return

	def compute_volForce(self, volfun, args=None):
		" Computes 'volumetric' source vector in 1D "
		if args is None: args={'position':self.part.qpPhy, 'time':self.time.qpPhy}
		prop  = np.kron(np.ones(self.time.nbqp), self.part.detJ)*volfun(args)
		force = sp.kron(self.time._denseweights[2], self.part._denseweights[0]) @ prop
		return force

	def compute_internalVariables(self, disp_cp, plastic_cp, plseq_cp, 
								pseudo1_cp, pseudo2_cp, lagrange_cp, threshold=1e-8, factor=1):
		# Compute internal variables
		E = self.mechamaterial.elasticModulus*np.ones(self.part.nbqp*self.time.nbqp)
		strain = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T @ disp_cp
		plastic = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T @ plastic_cp
		plseq = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T @ plseq_cp
		pseudo1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T @ pseudo1_cp
		pseudo2 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T @ pseudo2_cp
		lagrange = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T @ lagrange_cp

		stress = E*(strain - plastic)
		fyield = np.abs(pseudo1) + pseudo2 - self.mechamaterial.elasticLimit
		tmp = lagrange + factor*fyield
		gamma, heaviside = macaulayfunc(tmp), np.heaviside(tmp, threshold*self.mechamaterial.elasticLimit)
		K = self.mechamaterial._isoHardening._isohardfunders(plseq)

		# Compute flux of internal variables
		fluxplastic = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T @ plastic_cp
		fluxplseq = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T @ plseq_cp

		args = {'stress': stress, 'strain': strain,
				'plastic': plastic, 'fluxplastic': fluxplastic,
				'plseq': plseq, 'fluxplseq': fluxplseq,
				'gamma': gamma, 'heaviside': heaviside,
				'pseudo1': pseudo1, 'pseudo2': pseudo2,
				'Emod': E, 'Kmod': K, 'lagrange': lagrange,
				'factor': factor
				}
		return args

	def compute_MechStaticResidual1(self, Fext, args={}):
		stress = args.get('stress')
		res = (Fext - sp.kron(self.time._denseweights[2], self.part._denseweights[2]) @ stress)
		return res
	
	def compute_MechStaticResidual2(self, args={}):
		pseudo1 = args.get('pseudo1'); E = args.get('Emod'); strain = args.get('strain'); plastic = args.get('plastic')
		tmp2 = -pseudo1 + E*(strain - plastic) 
		res = sp.kron(self.time._denseweights[2], self.part._denseweights[0]) @ tmp2
		return res
	
	def compute_MechStaticResidual3(self, args={}):
		K = args.get('Kmod'); pseudo2 = args.get('pseudo2'); plseq = args.get('plseq')
		tmp2 = -pseudo2 - K*plseq
		res = sp.kron(self.time._denseweights[2], self.part._denseweights[0]) @ tmp2
		return res
	
	def compute_MechStaticResidual4(self, args={}):
		gamma = args.get('gamma'); pseudo1 = args.get('pseudo1'); fluxplastic = args.get('fluxplastic')
		tmp2 = gamma*np.sign(pseudo1) - fluxplastic
		res = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return res

	def compute_MechStaticResidual5(self, args={}):
		gamma = args.get('gamma'); fluxplseq = args.get('fluxplseq')
		tmp2 = gamma - fluxplseq
		res = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return res
	
	def compute_MechStaticResidual6(self, args={}):
		gamma = args.get('gamma'); lagrange = args.get('lagrange'); factor=args.get('factor')
		tmp2 = 1/factor*(gamma - lagrange)
		res = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return res
	
	# ==================================================
	def compute_TangentStiffness11(self, args={}):
		E = args.get('Emod'); Elocal = args.get('Elocal')
		prop = -E
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[3]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness12(self, args={}):
		E = args.get('Emod')
		prop = E
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[2]) @ tmp2
		return matrix.todense()	
	# ==================================================
	
	def compute_TangentStiffness21(self, args={}):
		E = args.get('Emod')
		prop = E
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[1]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness22(self, args={}):
		E = args.get('Emod')
		prop = -E
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness24(self, args={}):
		prop = -np.ones(self.time.nbqp*self.part.nbqp)
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	# ==================================================	
	
	def compute_TangentStiffness33(self, args={}):
		K = args.get('Kmod')
		prop = -K
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[0]) @ tmp2
		return matrix.todense()

	def compute_TangentStiffness35(self, args={}):
		prop = -np.ones(self.time.nbqp*self.part.nbqp)
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	# ==================================================	
	
	def compute_TangentStiffness42(self, args={}):
		prop = -np.ones(self.time.nbqp*self.part.nbqp)
		tmp1 = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[1], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness44(self, args={}):
		factor = args.get('factor'); H = args.get('heaviside')
		prop = factor*H
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness45(self, args={}):
		factor = args.get('factor'); H = args.get('heaviside'); pseudo1 = args.get('pseudo1')
		prop = factor*H*np.sign(pseudo1)
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness46(self, args={}):
		H = args.get('heaviside'); pseudo1 = args.get('pseudo1')
		prop = H*np.sign(pseudo1)
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()	
	# ==================================================	
	
	def compute_TangentStiffness53(self, args={}):
		prop = -np.ones(self.time.nbqp*self.part.nbqp)
		tmp1 = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[1], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness54(self, args={}):
		H = args.get('heaviside'); pseudo1 = args.get('pseudo1'); factor = args.get('factor')
		prop = H*np.sign(pseudo1)*factor
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()	
	
	def compute_TangentStiffness55(self, args={}):
		H = args.get('heaviside'); factor = args.get('factor')
		prop = H*factor
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness56(self, args={}):
		H = args.get('heaviside')
		prop = H
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	# ==================================================	
	
	def compute_TangentStiffness64(self, args={}):
		H = args.get('heaviside'); pseudo1 = args.get('pseudo1')
		prop = H*np.sign(pseudo1)
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness65(self, args={}):
		H = args.get('heaviside')
		prop = H
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def compute_TangentStiffness66(self, args={}):
		H = args.get('heaviside'); factor = args.get('factor')
		prop = 1/factor*(H - 1)
		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
		tmp2 = sp.diags(prop) @ tmp1
		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
		return matrix.todense()
	
	def solveViscoPlasticityProblem(self, dispinout, Fext, init=1e0):
		assert self.mechamaterial._chabocheNBparameters <= 1, 'Try another method'
		assert not np.any(self.mechamaterial._chabocheTable), 'Try another method'
		assert self.mechamaterial._isoHardening._isoname == 'linear', 'Try another method'

		boundary = deepcopy(self.boundary)
		boundary.clear_Dirichlet()
		table = np.zeros(shape=np.shape(self.boundary.thDirichletTable)); table[-1, 0] = 1
		boundary.add_DirichletConstTemperature(table)
		dof_intvar, dod_intvar = boundary.thdof, boundary.thdod
		
		nbctrlpts_total = self.boundary._nbctrlpts_total
		dof = self.boundary.thdof; dod = self.boundary.thdod
		disp_cp, plastic_cp, plseq_cp = np.copy(dispinout), np.zeros(nbctrlpts_total), np.zeros(nbctrlpts_total)
		pseudo1_cp, pseudo2_cp, lagrange_cp = np.zeros(nbctrlpts_total), np.zeros(nbctrlpts_total), np.zeros(nbctrlpts_total)
		
		factor = np.copy(init)
		for j in range(self._itersNL): # Newton-Raphson 
			
			# Compute internal variables
			# factor *= 10
			args = self.compute_internalVariables(disp_cp, plastic_cp, plseq_cp, 
												pseudo1_cp, pseudo2_cp, lagrange_cp, factor=factor)

			# Compute internal forces
			res1 = self.compute_MechStaticResidual1(Fext, args); res1[dod] = 0.0
			res2 = self.compute_MechStaticResidual2(args); res2[dod_intvar] = 0.0
			res3 = self.compute_MechStaticResidual3(args); res3[dod_intvar] = 0.0
			res4 = self.compute_MechStaticResidual4(args); res4[dod_intvar] = 0.0
			res5 = self.compute_MechStaticResidual5(args); res5[dod_intvar] = 0.0
			res6 = self.compute_MechStaticResidual6(args); res6[dod_intvar] = 0.0

			# Compute residue
			resNLj = np.sqrt(
						np.dot(res1, res1)+np.dot(res2, res2)+np.dot(res3, res3)
						+np.dot(res4, res4)+np.dot(res5, res5)+np.dot(res6, res6)
					)
			if j == 0: resNL0 = resNLj
			print('NonLinear error: %.5e' %resNLj)
			if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break

			# Solver for active control points
			K11 = self.compute_TangentStiffness11(args=args)[np.ix_(dof, dof)]
			K12 = self.compute_TangentStiffness12(args=args)[np.ix_(dof, dof_intvar)]
			#
			K21 = self.compute_TangentStiffness21(args=args)[np.ix_(dof_intvar, dof)]
			K22 = self.compute_TangentStiffness22(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K24 = self.compute_TangentStiffness24(args=args)[np.ix_(dof_intvar, dof_intvar)]
			#
			K33 = self.compute_TangentStiffness33(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K35 = self.compute_TangentStiffness35(args=args)[np.ix_(dof_intvar, dof_intvar)]
			#
			K42 = self.compute_TangentStiffness42(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K44 = self.compute_TangentStiffness44(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K45 = self.compute_TangentStiffness45(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K46 = self.compute_TangentStiffness46(args=args)[np.ix_(dof_intvar, dof_intvar)]
			#
			K53 = self.compute_TangentStiffness53(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K54 = self.compute_TangentStiffness54(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K55 = self.compute_TangentStiffness55(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K56 = self.compute_TangentStiffness56(args=args)[np.ix_(dof_intvar, dof_intvar)]
			#
			K64 = self.compute_TangentStiffness64(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K65 = self.compute_TangentStiffness65(args=args)[np.ix_(dof_intvar, dof_intvar)]
			K66 = self.compute_TangentStiffness66(args=args)[np.ix_(dof_intvar, dof_intvar)]
			
			ZZ1 = np.zeros((len(dof), len(dof_intvar)))
			ZZ2 = np.zeros((len(dof_intvar), len(dof)))
			ZZ3 = np.zeros((len(dof_intvar), len(dof_intvar)))

			if j == 0:
				K = -np.block([	[K11, K12, ZZ1, ZZ1, ZZ1, ZZ1], 
								[K21, K22, ZZ3, 0*K24, ZZ3, ZZ3], 
								[ZZ2, ZZ3, K33, ZZ3, 0*K35, ZZ3],
								[ZZ2, 0*K42, ZZ3, K44, K45, K46],
								[ZZ2, ZZ3, 0*K53, K54, K55, K56], 
								[ZZ2, ZZ3, ZZ3, K64, K65, K66]])
			else:
				K = -np.block([	[K11, K12, ZZ1, ZZ1, ZZ1, ZZ1], 
								[K21, K22, ZZ3, K24, ZZ3, ZZ3], 
								[ZZ2, ZZ3, K33, ZZ3, K35, ZZ3],
								[ZZ2, K42, ZZ3, K44, K45, K46],
								[ZZ2, ZZ3, K53, K54, K55, K56], 
								[ZZ2, ZZ3, ZZ3, K64, K65, K66]])

			F = np.hstack([res1[dof], res2[dof_intvar], res3[dof_intvar],
							res4[dof_intvar], res5[dof_intvar], res6[dof_intvar]])

			# sol = sclin.lstsq(K, F)[0]
			# sol = sp.linalg.lsqr(sp.csr_matrix(K), F)[0]
			sol = np.linalg.pinv(K) @ F

			deltaD1 = sol[:len(dof)]
			deltaD2 = sol[len(dof):len(dof)+len(dof_intvar)]
			deltaD3 = sol[len(dof)+len(dof_intvar):len(dof)+2*len(dof_intvar)]
			deltaD4 = sol[len(dof)+2*len(dof_intvar):len(dof)+3*len(dof_intvar)]
			deltaD5 = sol[len(dof)+3*len(dof_intvar):len(dof)+4*len(dof_intvar)]
			deltaD6 = sol[-len(dof_intvar):]

			# Update active control points
			disp_cp[dof] += deltaD1
			plastic_cp[dof_intvar] += deltaD2
			plseq_cp[dof_intvar] += deltaD3
			pseudo1_cp[dof_intvar] += deltaD4
			pseudo2_cp[dof_intvar] += deltaD5
			lagrange_cp[dof_intvar] += deltaD6
		return disp_cp

# class stmechaproblem1D(stproblem1D):

# 	def __init__(self, mechanical_material:mechamat, part:part1D, tspan:part1D, boundary:boundaryCondition, solverArgs={}):
# 		stproblem1D.__init__(self, part, tspan, boundary, solverArgs)
# 		self.mechamaterial = mechanical_material
# 		if self.mechamaterial.density is None: self.mechamaterial.density = lambda x: np.ones(self.part.nbqp_total)
# 		return

# 	def compute_volForce(self, volfun, args=None):
# 		" Computes 'volumetric' source vector in 1D "
# 		if args is None: args={'position':self.part.qpPhy, 'time':self.time.qpPhy}
# 		prop  = np.kron(np.ones(self.time.nbqp), self.part.detJ)*volfun(args)
# 		force = sp.kron(self.time._denseweights[2], self.part._denseweights[0]) @ prop
# 		return force
	
# 	def __evalyielfunction(self, stress, plseq):
# 		fyield = np.abs(stress) - self.mechamaterial._isoHardening._isohardfun(plseq)
# 		return fyield

# 	def compute_internalVariables(self, disp_cp, plseq_cp, plastic_cp, threshold=1e-8):
# 		# Compute internal variables
# 		E = self.mechamaterial.elasticModulus
# 		strain = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T @ disp_cp*np.kron(np.ones(self.time.nbqp), self.part.invJ)
# 		plastic = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T @ plastic_cp
# 		plseq = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T @ plseq_cp
# 		stress = E*(strain - plastic)
# 		fyield = self.__evalyielfunction(stress, plseq)

# 		# Compute flux of internal variables
# 		K = self.mechamaterial._isoHardening._isohardfunders(plseq)
# 		eta = self.mechamaterial.viscoparameter

# 		fluxstrain = sp.kron(self.time._densebasis[1], self.part._densebasis[1]).T @ disp_cp*np.kron(self.time.invJ, self.part.invJ)
# 		fluxplastic = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T @ plastic_cp * np.kron(self.time.invJ, np.ones(self.part.nbqp))
# 		fluxplseq = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T @ plseq_cp * np.kron(self.time.invJ, np.ones(self.part.nbqp))
# 		gamma = 1e-5/eta*macaulayfunc(fyield)
# 		Elocal = 1e-5*E/eta*np.heaviside(fyield, threshold*self.mechamaterial.elasticLimit)
# 		Klocal = 1e-5*K/eta*np.heaviside(fyield, threshold*self.mechamaterial.elasticLimit)
# 		args = {'stress': stress, 'strain': strain, 'fluxstrain': fluxstrain,
# 				'plastic': plastic, 'fluxplastic': fluxplastic,
# 				'plseq': plseq, 'fluxplseq': fluxplseq,
# 				'gamma': gamma, 'Elocal': Elocal, 'Klocal': Klocal,
# 				'Emod': E, 'Kmod': K,
# 				}
# 		return args

# 	# ==================================================
# 	def compute_MechStaticResidual1(self, Fext, args={}):
# 		fluxplastic = args.get('fluxplastic'); stress = args.get('stress')
# 		E = args.get('Emod'); gamma = args.get('gamma')
# 		res = (sp.kron(self.time._denseweights[2], self.part._denseweights[2]) @ stress 
# 			+ sp.kron(self.time._denseweights[0], self.part._denseweights[2]) @ (E*(fluxplastic - gamma*np.sign(stress))*np.kron(self.time.detJ, np.ones(self.part.nbqp)))
# 			- Fext
# 			)
# 		return res
	
# 	def compute_TangentStiffness11(self, args={}):
# 		E = args.get('Emod'); Elocal = args.get('Elocal')
# 		prop = E*np.kron(np.ones(self.time.nbqp), self.part.invJ)
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[-1]) @ tmp2

# 		prop = -Elocal*np.kron(self.time.detJ, self.part.invJ)
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix += sp.kron(self.time._denseweights[0], self.part._denseweights[-1]) @ tmp2
	
# 		return matrix.todense()
	
# 	def compute_TangentStiffness12(self, args={}):
# 		E = args.get('Emod'); Elocal = args.get('Elocal')
# 		prop = -E*np.ones(self.part.nbqp*self.time.nbqp)
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[2], self.part._denseweights[2]) @ tmp2

# 		prop = E*np.ones(self.part.nbqp*self.time.nbqp)
# 		tmp1 = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix += sp.kron(self.time._denseweights[1], self.part._denseweights[2]) @ tmp2
		
# 		prop = Elocal*np.kron(self.time.detJ, np.ones(self.part.nbqp))
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix += sp.kron(self.time._denseweights[0], self.part._denseweights[2]) @ tmp2
		
# 		return matrix.todense()
	
# 	def compute_TangentStiffness13(self, args={}):
# 		K = args.get('Kmod'); Elocal = args.get('Elocal'); stress = args.get('stress')
# 		prop = K*Elocal*np.sign(stress)*np.kron(self.time.detJ, np.ones(self.part.nbqp))
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[2]) @ tmp2
# 		return matrix.todense()
	
# 	# ==================================================	
# 	def compute_MechStaticResidual2(self, args={}):
# 		gamma = args.get('gamma'); stress = args.get('stress'); fluxplastic = args.get('fluxplastic'); E = args.get('Emod')
# 		tmp2 = E*(gamma*np.sign(stress) - fluxplastic)*np.kron(self.time.detJ, self.part.detJ)
# 		res = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
# 		return res
	
# 	def compute_TangentStiffness21(self, args={}):
# 		Elocal = args.get('Elocal')
# 		prop = Elocal*np.kron(self.time.detJ, np.ones(self.part.nbqp))
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[1]) @ tmp2
# 		return matrix.todense()
	
# 	def compute_TangentStiffness22(self, args={}):
# 		E = args.get('Emod'); Elocal = args.get('Elocal')
# 		prop = -E*np.kron(np.ones(self.time.nbqp), self.part.detJ)
# 		tmp1 = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[1], self.part._denseweights[0]) @ tmp2

# 		prop = -Elocal*np.kron(self.time.detJ, self.part.detJ)
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix += sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
# 		return matrix.todense()
	
# 	def compute_TangentStiffness23(self, args={}):
# 		stress = args.get('stress'); K = args.get('Kmod'); Elocal = args.get('Elocal')
# 		prop = -K*Elocal*np.sign(stress)*np.kron(self.time.detJ, self.part.detJ)
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
# 		return matrix.todense()
	
# 	# ==================================================	
# 	def compute_MechStaticResidual3(self, args={}):
# 		K = args.get('Kmod'); fluxplseq = args.get('fluxplseq'); gamma = args.get('gamma')
# 		tmp2 = K*(gamma - fluxplseq)*np.kron(self.time.detJ, self.part.detJ)
# 		res = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
# 		return res
	
# 	def compute_TangentStiffness31(self, args={}):
# 		stress = args.get('stress'); Klocal = args.get('Klocal')
# 		prop = Klocal*np.sign(stress)*np.kron(self.time.detJ, np.ones(self.part.nbqp))
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[1]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[1]) @ tmp2
# 		return matrix.todense()
	
# 	def compute_TangentStiffness32(self, args={}):
# 		stress = args.get('stress'); Klocal = args.get('Klocal')
# 		prop = -Klocal*np.sign(stress)*np.kron(self.time.detJ, self.part.detJ)

# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
# 		return matrix.todense()
	
# 	def compute_TangentStiffness33(self, args={}):

# 		Klocal = args.get('Klocal'); K = args.get('Kmod')
# 		prop = -K*np.kron(np.ones(self.time.nbqp), self.part.detJ)
# 		tmp1 = sp.kron(self.time._densebasis[1], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix = sp.kron(self.time._denseweights[1], self.part._denseweights[0]) @ tmp2

# 		prop = -Klocal*K*np.kron(self.time.detJ, self.part.detJ)
# 		tmp1 = sp.kron(self.time._densebasis[0], self.part._densebasis[0]).T
# 		tmp2 = sp.diags(prop) @ tmp1
# 		matrix += sp.kron(self.time._denseweights[0], self.part._denseweights[0]) @ tmp2
# 		return matrix.todense()
	
# 	def solveViscoPlasticityProblem(self, dispinout, Fext):
# 		assert self.mechamaterial._chabocheNBparameters <= 1, 'Try another method'
# 		assert not np.any(self.mechamaterial._chabocheTable), 'Try another method'
# 		assert self.mechamaterial._isoHardening._isoname == 'linear', 'Try another method'

# 		boundary = deepcopy(self.boundary)
# 		boundary.clear_Dirichlet()
# 		table = np.zeros(shape=np.shape(self.boundary.thDirichletTable)); table[-1, 0] = 1
# 		boundary.add_DirichletConstTemperature(table)
# 		dof_intvar, dod_intvar = boundary.thdof, boundary.thdod
		
# 		nbctrlpts_total = self.boundary._nbctrlpts_total
# 		dof = self.boundary.thdof; dod = self.boundary.thdod
# 		disp_cp, plastic_cp, plseq_cp = np.zeros(nbctrlpts_total), np.zeros(nbctrlpts_total), np.copy(dispinout)
		
# 		for j in range(self._itersNL): # Newton-Raphson 
			
# 			# Compute internal variables
# 			args = self.compute_internalVariables(disp_cp, plseq_cp, plastic_cp)

# 			# Compute internal forces
# 			res1 = self.compute_MechStaticResidual1(Fext, args); res1[dod] = 0.0
# 			res2 = self.compute_MechStaticResidual2(args); res2[dod_intvar] = 0.0
# 			res3 = self.compute_MechStaticResidual3(args); res3[dod_intvar] = 0.0

# 			# Compute residue
# 			resNLj = np.sqrt(np.dot(res1, res1)+np.dot(res2, res2)+np.dot(res3, res3))
# 			if j == 0: resNL0 = resNLj
# 			print('NonLinear error: %.5e' %resNLj)
# 			if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break

# 			# Solver for active control points
# 			K11 = self.compute_TangentStiffness11(args=args)[np.ix_(dof, dof)]
# 			K12 = self.compute_TangentStiffness12(args=args)[np.ix_(dof, dof_intvar)]
# 			K13 = self.compute_TangentStiffness13(args=args)[np.ix_(dof, dof_intvar)]
# 			K21 = self.compute_TangentStiffness21(args=args)[np.ix_(dof_intvar, dof)]
# 			K22 = self.compute_TangentStiffness22(args=args)[np.ix_(dof_intvar, dof_intvar)]
# 			K23 = self.compute_TangentStiffness23(args=args)[np.ix_(dof_intvar, dof_intvar)]
# 			K31 = self.compute_TangentStiffness31(args=args)[np.ix_(dof_intvar, dof)]
# 			K32 = self.compute_TangentStiffness32(args=args)[np.ix_(dof_intvar, dof_intvar)]
# 			K33 = self.compute_TangentStiffness33(args=args)[np.ix_(dof_intvar, dof_intvar)]

# 			K = -np.block([	[K11, K12, K13], 
# 							[K21, K22, K23], 
# 							[K31, K32, K33]])
# 			F = np.hstack([res1[dof], res2[dof_intvar], res3[dof_intvar]])
# 			# sol = sp.linalg.spsolve(sp.csr_matrix(K), F)
# 			sol = sp.linalg.gmres(sp.csr_matrix(K), F)[0]

# 			deltaD1 = sol[:len(dof)]
# 			deltaD2 = sol[len(dof):len(dof)+len(dof_intvar)]
# 			deltaD3 = sol[-len(dof_intvar):]

# 			# Update active control points
# 			disp_cp[dof] += deltaD1
# 			plastic_cp[dof_intvar] += deltaD2
# 			plseq_cp[dof_intvar] += deltaD3

# 		return disp_cp
