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
from .lib_base import evalDersBasisFortran, array2csr_matrix

class problem1D():
	def __init__(self, part:part1D, boundary:boundaryCondition, solverArgs={}):
		self.part = part
		self.boundary = boundary
		self.addSolverConstraints(solverArgs=solverArgs)
		return
	
	def addSolverConstraints(self, solverArgs:dict):
		self._thresNL = solverArgs.get('thres_nonlinear', 1e-10)
		self._itersNL = solverArgs.get('iters_nonlinear', 20)
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

		problem_ref = normArgs.get('part_ref', None); u_ref = normArgs.get('u_ref', None)
		if isinstance(problem_ref, problem1D) and isinstance(u_ref, np.ndarray):
			denseBasisExact = []
			basis_csr, indi_csr, indj_csr = evalDersBasisFortran(problem_ref.part.degree, problem_ref.part.knotvector, quadPts)
			for i in range(2): denseBasisExact.append(array2csr_matrix(basis_csr[:, i], indi_csr, indj_csr))
			u_exact = denseBasisExact[0].T @ u_ref
			invJExact = denseBasisExact[1].T @ problem_ref.part.ctrlpts
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

	def solveFourierTransientProblem(self, Tinout, Fext_list, time_list, alpha=0.5, isLumped=False):
		" Solves transient heat problem in 1D. "

		def computeVelocity(problem:heatproblem1D, Fext, args=None, isLumped=False):
			if args is None: args = {'position': problem.part.qpPhy}
			prop = problem.heatmaterial.capacity(args)*problem.heatmaterial.density(args)*problem.part.detJ
			matrix = problem.part._denseweights[0] @ np.diag(prop) @ problem.part._densebasis[0].T
			if isLumped: matrix = np.diag(matrix.sum(axis=1))
			velocity = np.linalg.solve(matrix, Fext)
			return velocity

		nbctrlpts_total = self.part.nbctrlpts_total; nsteps = len(time_list)
		dod = self.boundary.thdod; dof = self.boundary.thdof

		# Compute initial velocity using interpolation
		assert nsteps > 2, 'At least 2 steps'
		temperature = self.interpolate_temperature(Tinout[:, 0])
		args = {'temperature':temperature, 'position':self.part.qpPhy}
		tmp = (Fext_list[:, 0] - self.compute_mfConductivity(Tinout[:, 0], args=args))
		V_n0 = computeVelocity(self, tmp, args=args, isLumped=isLumped)

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
				print('Nonlinear error %.5e' %resNLj)
				if resNLj <= max([self._thresNL*resNL0, 1e-12]): break
				
				# Solver for active control points
				tangentM = sp.csr_matrix(self.compute_FourierTangentMatrix(dt, alpha=alpha, args=args, isLumped=isLumped)[np.ix_(dof, dof)])
				deltaV = np.zeros(nbctrlpts_total); deltaV[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update active control points
				dj_n1 += alpha*dt*deltaV
				Vj_n1 += deltaV

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
				print('Nonlinear error %.5e' %resNLj)
				if resNLj <= max([self._thresNL*resNL0, 1e-12]): break

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
		self.mechamat = mechanical_material		
		if self.mechamat.density is None: self.mechamat.density = lambda x: np.ones(self.part.nbqp)
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
		Fint = self.part._denseweights[-1] @ np.ravel(stress).T
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

		nbChaboche = self.mechamat._chabocheNBparameters
		nbqp = self.part.nbqp; dof = self.boundary.thdof; dod = self.boundary.thdod
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

				resNLj = np.sqrt(np.dot(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('Nonlinear error %.5e' %resNLj)
				if resNLj <= max([self._thresNL*resNL0, 1e-12]): break
				if j > 0 and isElasticLoad: break

				# Solver for active control points
				tangentM = sp.csr_matrix(self.compute_PlasticityTangentMatrix(Cep)[np.ix_(dof, dof)])
				deltaD = np.zeros(self.part.nbctrlpts_total); deltaD[dof] = sp.linalg.spsolve(tangentM, r_dj[dof])

				# Update active control points
				dj_n1 += deltaD

			dispinout[:, i] = dj_n1
			Allstrain[:, i] = np.ravel(strain)
			Allstress[:, i] = np.ravel(stress)
			Allplseq[:, i]  = a_n1
			AllCep[:, i]    = np.ravel(Cep)

			pls_n0, a_n0, b_n0 = np.copy(pls_n1), np.copy(a_n1), np.copy(b_n1)

		return Allstrain, Allstress, Allplseq, AllCep
