from .__init__ import *
from .lib_base import get_faceInfo, get_INCTable, evalDersBasisFortran
from .lib_quadrules import GaussQuadrature
from .lib_material import (thermomat, mechamat,
							clean_dirichlet, block_dot_product)
from .lib_part import part
from .lib_boundary import boundaryCondition

class problem():
	def __init__(self, part:part, boundary:boundaryCondition, solverArgs:dict):
		self.material = None
		self.part     = part
		self.boundary = boundary
		self.addSolverConstraints(solverArgs)
		return
	
	def _getInputs(self):
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, *self.part.weights]
		return inpts

	def addSolverConstraints(self, solverArgs:dict):
		self._nbIterPCG    = solverArgs.get('nbIterationsPCG', 50)
		self._nbIterNR     = solverArgs.get('nbIterationsNR', 30)
		self._thresholdPCG = solverArgs.get('PCGThreshold', 1e-10)
		self._thresholdNR  = solverArgs.get('NRThreshold', 1e-8)
		self._methodPCG    = solverArgs.get('PCGmethod', 'JMC')
		return
	
	def __getInfo4surfForce(self, nbFacePosition):
		INC_ctrlpts = get_INCTable(self.part.nbctrlpts)
		direction, side = get_faceInfo(nbFacePosition)
		if direction>=2*self.part.dim: raise Warning('Not possible')
		valrange = [i for i in range(self.part.dim)]
		valrange.pop(direction)
		if side == 0:   CPList = np.where(INC_ctrlpts[:, direction] == 0)[0]
		elif side == 1: CPList = np.where(INC_ctrlpts[:, direction] == self.part.nbctrlpts[direction]-1)[0]
		CPList = list(np.sort(CPList))

		nnz, indices, basis, weights = [], [], [], []
		for _ in valrange:
			nnz.append(self.part.nbqp[_]); basis.append(self.part.basis[_]); weights.append(self.part.weights[_])
			indices.append(self.part.indices[2*_]); indices.append(self.part.indices[2*_+1]) 
		inpts = [*nnz, *indices, *basis, self.part.ctrlpts[:, CPList]]
		if self.part.dim == 2: 
			Jqp = geophy.eval_jacobien_1d(*inpts)
			qpPhy = geophy.interpolate_meshgrid_1d(*inpts)
		elif self.part.dim == 3:
			Jqp = geophy.eval_jacobien_2d(*inpts)
			qpPhy = geophy.interpolate_meshgrid_2d(*inpts)
		return nnz, indices, basis, weights, Jqp, qpPhy, CPList

	def compute_surfForce(self, surffun, nbFacePosition):
		""" Computes the surface foce over the boundary of a geometry. 
			The surffun is a Neumann like function, ie, in transfer heat q = - (k grad(T)).normal
			and in elasticity t = sigma.normal
		"""
		nnz, indices, _, weights, Jqp, qpPhy, CPList = self.__getInfo4surfForce(nbFacePosition)
		prop = surffun(qpPhy); prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
		inpts = [*nnz, *indices, *weights, Jqp, prop]
		if   self.part.dim == 2: tmp = geophy.get_forcesurf_2d(*inpts)
		elif self.part.dim == 3: tmp = geophy.get_forcesurf_3d(*inpts)
		surfForce = np.zeros((nr, self.part.nbctrlpts_total))
		surfForce[:, CPList] = tmp
		if nr == 1: surfForce = np.ravel(surfForce)
		return surfForce, CPList
	
	def compute_surfForce_woNormal(self, gradfun, nbFacePosition):
		nnz, indices, _, weights, Jqp, qpPhy, CPList = self.__getInfo4surfForce(nbFacePosition)

		normal_qpPhy = geophy.eval_normal(Jqp)
		if np.any(np.array([0, 3, 4], dtype=int) == nbFacePosition): normal_qpPhy = -normal_qpPhy 
		# By the moment I do not why in this surfaces the normal has the wrong sign (the same for 2 and 3 dimensions)
		grad_qpPhy = gradfun(qpPhy)
		if np.array(grad_qpPhy).ndim == 3:
			prop = np.einsum('ijl,jl->il', grad_qpPhy, normal_qpPhy)
		elif np.array(grad_qpPhy).ndim == 2:
			prop = np.einsum('jl,jl->l', grad_qpPhy, normal_qpPhy)
		else: raise Warning('Size problem')

		prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
		inpts = [*nnz, *indices, *weights, Jqp, prop]
		if   self.part.dim == 2: tmp = geophy.get_forcesurf_2d(*inpts)
		elif self.part.dim == 3: tmp = geophy.get_forcesurf_3d(*inpts)
		surfForce = np.zeros((nr, self.part.nbctrlpts_total))
		surfForce[:, CPList] = tmp
		if nr == 1: surfForce = np.ravel(surfForce)
		return surfForce, CPList
	
	def compute_volForce(self, volfun): 
		" Computes the volume force over a geometry "
		prop = volfun(self.part.qpPhy)
		prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.detJ, prop]
		if self.part.dim == 2: volForce = geophy.get_forcevol_2d(*inpts)
		if self.part.dim == 3: volForce = geophy.get_forcevol_3d(*inpts)
		if nr == 1: volForce = np.ravel(volForce)
		return volForce	
	
	def L2NormOfError(self, u_ctrlpts, L2NormArgs:dict):
		""" Computes the norm L2 of the error. The exactfun is the function of the exact solution. 
			and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
			whether the default quadrature is weighted quadrature. 
		"""	
		# Compute u interpolation
		nr = np.size(np.atleast_2d(u_ctrlpts), axis=0)
		nbqp, quadPts, indices, basis, parametricWeights = [], [], [], [], []
		for i in range(self.part.dim):
			quadRule = GaussQuadrature(self.part.degree[i], self.part.knotvector[i], quadArgs={'type':'leg'})
			quadPtsByDir, indicesByDir, basisByDir, _ = quadRule.getQuadratureRulesInfo()
			indi, indj = indicesByDir; parweightsByDir = quadRule._parametricWeights
			
			nbqp.append(quadRule.nbqp); quadPts.append(quadPtsByDir); indices.append(indi); indices.append(indj)
			basis.append(basisByDir); parametricWeights.append(parweightsByDir)

		inpts = [*nbqp, *indices, *basis]
		if self.part.dim == 2:
			Jqp = geophy.eval_jacobien_2d(*inpts, self.part.ctrlpts)
			detJ, _ = geophy.eval_inverse_det(Jqp)
			qpPhy = geophy.interpolate_meshgrid_2d(*inpts, self.part.ctrlpts)
			u_interp = geophy.interpolate_meshgrid_2d(*inpts, np.atleast_2d(u_ctrlpts))
		if self.part.dim == 3:
			Jqp = geophy.eval_jacobien_3d(*inpts, self.part.ctrlpts)
			detJ, _ = geophy.eval_inverse_det(Jqp)
			qpPhy = geophy.interpolate_meshgrid_3d(*inpts, self.part.ctrlpts)
			u_interp = geophy.interpolate_meshgrid_3d(*inpts, np.atleast_2d(u_ctrlpts))
		if nr == 1: u_interp = np.ravel(u_interp)

		# Compute u exact
		u_exact  = None
		exactfun = L2NormArgs.get('exactFunction', None)
		if callable(exactfun): u_exact = exactfun(qpPhy)
		refPart  = L2NormArgs.get('referencePart', None); u_ref = L2NormArgs.get('u_ref', None)
		if isinstance(refPart, part) and isinstance(u_ref, np.ndarray):
			nbqpExact, basisExact, indicesExact = [], [], []
			for i in range(self.part.dim):
				basis, indi, indj = evalDersBasisFortran(refPart.degree[i], refPart.knotvector[i], quadPts[i])
				nbqpExact.append(len(quadPts[i])); basisExact.append(basis); indicesExact.append(indi); indicesExact.append(indj)
			inpts = [*nbqpExact, *indicesExact, *basisExact, np.atleast_2d(u_ref)]
			if self.part.dim == 2:   u_exact = geophy.interpolate_meshgrid_2d(*inpts)    
			elif self.part.dim == 3: u_exact = geophy.interpolate_meshgrid_3d(*inpts)
			if nr == 1: u_exact = np.ravel(u_exact)
		if u_exact is None: raise Warning('Not possible')

		# Compute error
		ue_diff_uh2 = (u_exact - u_interp)**2 
		ue2         = u_exact**2
		if nr > 1: 
			ue_diff_uh2 = np.ravel(np.sum(ue_diff_uh2, axis=0))
			ue2         = np.ravel(np.sum(ue2, axis=0))
		ue_diff_uh2 = ue_diff_uh2 * detJ
		ue2         = ue2*detJ

		ue_diff_uh2 = np.reshape(ue_diff_uh2, tuple(nbqp), order='F')
		ue2         = np.reshape(ue2, tuple(nbqp), order='F')
		if self.part.dim == 2: 
			tmp1 = np.einsum('i,j,ij->', parametricWeights[0], parametricWeights[1], ue2)
			tmp2 = np.einsum('i,j,ij->', parametricWeights[0], parametricWeights[1], ue_diff_uh2)
		if self.part.dim == 3: 
			tmp1 = np.einsum('i,j,k,ijk->', parametricWeights[0], parametricWeights[1], parametricWeights[2], ue2)
			tmp2 = np.einsum('i,j,k,ijk->', parametricWeights[0], parametricWeights[1], parametricWeights[2], ue_diff_uh2)
		error = np.sqrt(tmp2/tmp1)

		return error

	def L2projectionCtrlptsVol(self, volfun):
		" Given the solution field (function) over a physical space, it computes the L2 projection, ie. the value at control points. "
		volForce = self.compute_volForce(volfun); volForce = np.atleast_2d(volForce); nr = np.size(volForce, axis=0)
		inpts = [*self._getInputs(), self.part.detJ, volForce, self._nbIterPCG, self._thresholdPCG]
		if self.part.dim == 2: u_interp, _ = geophy.l2projection_ctrlpts_2d(*inpts)
		if self.part.dim == 3: u_interp, _ = geophy.l2projection_ctrlpts_3d(*inpts)
		if nr == 1: u_interp = np.ravel(u_interp)
		return u_interp
	
class heatproblem(problem):
	def __init__(self, material:thermomat, part:part, boundary:boundaryCondition, solverArgs={}):
		super().__init__(part, boundary, solverArgs)
		self.material = material
		return
	
	def compute_mfConductivity(self, array_in, args=None):
		if args is None: args = self.part.qpPhy
		prop = self.material.conductivity(args)
		inpts = [*super()._getInputs(), self.part.invJ, self.part.detJ, prop]
		if self.part.dim == 2: array_out = heatsolver.mf_get_ku_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = heatsolver.mf_get_ku_3d(*inpts, array_in)
		return array_out
	
	def compute_mfCapacity(self, array_in, args=None, isLumped=False): 
		if args is None: args = self.part.qpPhy
		prop = self.material.capacity(args)
		inpts = [*super()._getInputs(), isLumped, self.part.invJ, self.part.detJ, prop]
		if self.part.dim == 2: array_out = heatsolver.mf_get_cu_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = heatsolver.mf_get_cu_3d(*inpts, array_in)
		return array_out

	def interpolate_temperature(self, T_ctrlpts):
		" Computes the temperature at the quadrature points "
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, np.atleast_2d(T_ctrlpts)]
		if self.part.dim == 2:   T_interp = geophy.interpolate_meshgrid_2d(*inpts)
		elif self.part.dim == 3: T_interp = geophy.interpolate_meshgrid_3d(*inpts)
		T_interp = np.ravel(T_interp)
		return T_interp

	def solveSteadyHeatProblemFT(self, Fext, args=None):
		dod = deepcopy(self.boundary.thdod) + 1
		if args is None: args = self.part.qpPhy
		prop = self.material.conductivity(args)
		inpts = [*super()._getInputs(), dod, self.boundary.thDirichletTable, self.part.invJ, self.part.detJ, 
				prop, Fext,  self._nbIterPCG, self._thresholdPCG, self._methodPCG]
		if self.part.dim == 2: temperature, residue = heatsolver.solver_steady_heat_2d(*inpts)
		if self.part.dim == 3: temperature, residue = heatsolver.solver_steady_heat_3d(*inpts)
		return temperature, residue
	
	def solveLinearTransientHeatProblemFT(self, Fext, thetadt, args=None, isLumped=False):
		dod = deepcopy(self.boundary.thdod) + 1
		if args is None: args = self.part.qpPhy
		Cprop = self.material.capacity(args)
		Kprop = self.material.conductivity(args)
		inpts = [*super()._getInputs(), isLumped, dod, self.boundary.thDirichletTable, self.part.invJ, self.part.detJ,
				Cprop, Kprop, thetadt, Fext, self._nbIterPCG, self._thresholdPCG, self._methodPCG]
		if self.part.dim == 2: temperature, residue = heatsolver.solver_lineartransient_heat_2d(*inpts)
		if self.part.dim == 3: temperature, residue = heatsolver.solver_lineartransient_heat_3d(*inpts)
		return temperature, residue

	def solveNLTransientHeatProblemPy(self, Tinout, time_list, Fext_list, theta=1.0, isLumped=False):
		nbctrlpts_total = self.part.nbctrlpts_total; nsteps = len(time_list)
		dod = self.boundary.getThermalBoundaryConditionInfo()[0]

		# Compute inital velocity using interpolation
		V_n0 = np.zeros(nbctrlpts_total)
		if nsteps == 2:
			dt = time_list[1] - time_list[0]
			V_n0[dod] = 1.0/dt*(Tinout[dod, 1] - Tinout[dod, 0])
		elif nsteps > 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			V_n0[dod] = 1.0/(dt1*(factor - factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
		else: raise Warning('At least 2 steps')
		
		AllresPCG = []
		for i in range(1, nsteps):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			d_n0 = np.copy(Tinout[:, i-1])

			# Get values of new step
			d0_n1 = d_n0 + dt*(1.0 - theta)*V_n0
			V_n1  = np.zeros(nbctrlpts_total); V_n1[dod]  = 1.0/theta*(1.0/dt*(Tinout[dod, i] - Tinout[dod, i-1]) - (1 - theta)*V_n0[dod])
			Fext_n1 = Fext_list[:, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR):

				# Compute temperature at each quadrature point
				dj_n1 = d0_n1 + theta*dt*V_n1
				temperature = self.interpolate_temperature(dj_n1)
			
				# Compute internal force
				args = np.row_stack((self.part.qpPhy, temperature))
				Fint_dj = self.compute_mfCapacity(V_n1, args=args, isLumped=isLumped) + self.compute_mfConductivity(dj_n1, args=args)

				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				# Iterative solver
				resPCGj = np.array([i, j+1])
				deltaV, resPCG = self.solveLinearTransientHeatProblemFT(r_dj, theta*dt, args=args, isLumped=isLumped)
				resPCGj = np.append(resPCGj, resPCG); AllresPCG.append(resPCGj)

				# Update values
				V_n1 += deltaV

				# Compute residue of Newton Raphson using an energetic approach
				resNRj = abs(theta*dt*np.dot(V_n1, r_dj))
				if j == 0: resNR0 = resNRj
				print('NR error: %.5e' %resNRj)
				if j>0 and resNRj<=self._thresholdNR*resNR0: break

			Tinout[:, i] = np.copy(dj_n1)
			V_n0 = np.copy(V_n1)

		return AllresPCG

class mechaproblem(problem):
	def __init__(self, material:mechamat, part:part, boundary:boundaryCondition, solverArgs={}):
		super().__init__(part, boundary, solverArgs)
		self.material = material
		return
	
	def compute_mfStiffness(self, array_in, mechArgs=None):
		if mechArgs is None: 
			dimen = self.part.dim; nvoigt = int(dimen*(dimen+1)/2)
			mechArgs = np.zeros((nvoigt+3, self.part.nbqp_total))
			mechArgs[0, :] = self.material.lame_lambda
			mechArgs[1, :] = self.material.lame_mu
		inpts = [*super()._getInputs(), self.part.invJ, self.part.detJ, mechArgs]
		if   self.part.dim == 2: array_out = plasticitysolver.mf_get_su_2d(*inpts, array_in)
		elif self.part.dim == 3: array_out = plasticitysolver.mf_get_su_3d(*inpts, array_in)
		return array_out
	
	def interpolate_strain(self, displacement):
		" Compute strain field from displacement field "
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, self.part.invJ, displacement]
		if   self.part.dim == 2: strain = plasticitysolver.interpolate_strain_2d(*inpts)
		elif self.part.dim == 3: strain = plasticitysolver.interpolate_strain_3d(*inpts)
		return strain
	
	def compute_intForce(self, stress):
		"Compute internal force using sigma coefficients "
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.invJ, self.part.detJ, stress]
		if   self.part.dim == 2: intForce = plasticitysolver.get_intforce_2d(*inpts)
		elif self.part.dim == 3: intForce = plasticitysolver.get_intforce_3d(*inpts)
		return intForce

	def solveElasticityProblemFTm(self, Fext, mechArgs=None):
		dod = deepcopy(self.boundary.mchdod)
		for i, tmp in enumerate(dod):
			tmp = tmp + 1; dod[i] = tmp

		dimen  = self.part.dim
		nvoigt = int(dimen*(dimen+1)/2)
		prop = [self.material.elasticmodulus, self.material.poissonratio, self.material.elasticlimit]
		if mechArgs is None:
			mechArgs = np.zeros((nvoigt+3, self.part.nbqp_total))
			mechArgs[0, :] = self.material.lame_lambda
			mechArgs[1, :] = self.material.lame_mu
		inpts = [*super()._getInputs(), *dod, self.boundary.mchDirichletTable, 
				self.part.invJ, self.part.detJ, prop, mechArgs, Fext, self._nbIterPCG, self._thresholdPCG, self._methodPCG]
		if   self.part.dim == 2: displacement, resPCG = plasticitysolver.solver_elasticity_2d(*inpts)
		elif self.part.dim == 3: displacement, resPCG = plasticitysolver.solver_elasticity_3d(*inpts)
		
		return displacement, resPCG
	
	def solvePlasticityProblemPy(self, Fext_list): 

		if not self.material._isPlasticityPossible: raise Warning('Plasticity not defined')
		dimen  = self.part.dim
		nvoigt = int(dimen*(dimen+1)/2)
		nbqp_total = self.part.nbqp_total
		nsteps = np.shape(Fext_list)[2]

		# Internal variables
		pls_n0 = np.zeros((nvoigt, nbqp_total))
		a_n0   = np.zeros(nbqp_total)
		b_n0   = np.zeros((nvoigt, nbqp_total))
		pls_n1 = np.zeros((nvoigt, nbqp_total))
		a_n1   = np.zeros(nbqp_total)
		b_n1   = np.zeros((nvoigt, nbqp_total))
		stress = np.zeros((nvoigt, nbqp_total))
		mechArgs = np.zeros((nvoigt+3, nbqp_total))
		
		# Output variables
		Alldisplacement = np.zeros(np.shape(Fext_list))
		Allstress 		= np.zeros((nvoigt, nbqp_total, nsteps))
		Allstrain 		= np.zeros((nvoigt, nbqp_total, nsteps))
		AllresPCG 		= []

		for i in range(1, nsteps):
			
			# Get values of last step
			d_n0  = Alldisplacement[:, :, i-1]

			# Get values of new step
			V_n1    = np.zeros(np.shape(d_n0))
			Fext_n1 = Fext_list[:, :, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR):

				# Compute strain at each quadrature point
				dj_n1  = d_n0 + V_n1
				strain = self.interpolate_strain(dj_n1)
	
				# Closest point projection in perfect plasticity
				output = self.material.returnMappingAlgorithm(strain, pls_n0, a_n0, b_n0)
				stress, pls_n1, a_n1 = output[:nvoigt, :], output[nvoigt:2*nvoigt, :], output[2*nvoigt, :]
				b_n1, mechArgs = output[2*nvoigt+1:3*nvoigt+1, :], output[3*nvoigt+1:, :]

				# Compute internal force 
				Fint_dj = self.compute_intForce(stress)
				
				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				clean_dirichlet(r_dj, self.boundary.mchdod) 
				
				# Iterative solver
				resPCGj = np.array([i, j+1])
				deltaV, resPCG = self.solveElasticityProblemFTm(Fext=r_dj, mechArgs=mechArgs)
				resPCGj = np.append(resPCGj, resPCG); AllresPCG.append(resPCGj)

				# Update values
				V_n1 += deltaV
				
				# Compute residue of Newton Raphson using an energetic approach
				resNRj = abs(block_dot_product(dimen, V_n1, r_dj))
				if j == 0: resNR0 = resNRj
				print('NR error: %.5e' %resNRj)
				if j>0 and resNRj<=self._thresholdNR*resNR0: break

			Alldisplacement[:, :, i] = dj_n1
			Allstress[:, :, i] = stress	
			Allstrain[:, :, i] = strain

			pls_n0 = np.copy(pls_n1)
			a_n0 = np.copy(a_n1)
			b_n0 = np.copy(b_n1)

		return Alldisplacement, AllresPCG, {'stress': Allstress, 'totalstrain': Allstrain}

