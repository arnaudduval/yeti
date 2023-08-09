from .__init__ import *
from .lib_base import eraseRowsCSR
from .lib_quadrules import GaussQuadrature
from .lib_material import (thermomat, mechamat, 
							clean_dirichlet, block_dot_product)
from .lib_part import part
from .lib_boundary import boundaryCondition, get_INCTable

def get_faceInfo(nb):
	direction = int(np.floor(nb/2))
	if nb%2 == 1: side = 1
	else: side = 0
	return direction, side

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

	def eval_volForce(self, volfun): 
		prop = volfun(self.part.qpPhy)
		prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.detJ, prop]
		if self.part.dim == 2: volForce = geophy.get_forcevol_2d(*inpts)
		if self.part.dim == 3: volForce = geophy.get_forcevol_3d(*inpts)
		if nr == 1: volForce = np.ravel(volForce)
		return volForce	
	
	def eval_surfForce(self, surffun, nbFacePosition):

		INC_ctrlpts = get_INCTable(self.part.nbctrlpts)
		direction, side = get_faceInfo(nbFacePosition)

		if direction>=2*self.part.dim: raise Warning('Not possible')
		valrange = [i for i in range(self.part.dim)]
		valrange.pop(direction)

		# Get control points and quadrature points list
		if side == 0: CPList = np.where(INC_ctrlpts[:, direction] == 0)[0]
		elif side == 1:  CPList = np.where(INC_ctrlpts[:, direction] == self.part.nbctrlpts[direction]-1)[0]
		
		CPList = list(np.sort(CPList))
		CtrlPts = self.part.ctrlpts[:, CPList]

		nnz, indices, basis, weights = [], [], [], []
		for _ in valrange:
			nnz.append(self.part.nbqp[_]); basis.append(self.part.basis[_]); weights.append(self.part.weights[_])
			indices.append(self.part.indices[2*_]); indices.append(self.part.indices[2*_+1]) 

		inpts = [*nnz, *indices, *basis, CtrlPts]
		if self.part.dim == 2: 
			JJ = geophy.eval_jacobien_1d(*inpts)
			qpPhy = geophy.interpolate_meshgrid_1d(*inpts)
		elif self.part.dim == 3:
			JJ = geophy.eval_jacobien_2d(*inpts)
			qpPhy = geophy.interpolate_meshgrid_2d(*inpts)

		prop = surffun(qpPhy); prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
		inpts = [*nnz, *indices, *weights, JJ, prop]
		if   self.part.dim == 2: tmp = geophy.get_forcesurf_2d(*inpts)
		elif self.part.dim == 3: tmp = geophy.get_forcesurf_3d(*inpts)
		surfForce = np.zeros((nr, self.part.nbctrlpts_total))
		surfForce[:, CPList] = tmp
		if nr == 1: surfForce = np.ravel(surfForce)
		return surfForce
	
	def L2NormOfError(self, exactfun, u_ctrlpts):
		""" Computes the norm L2 of the error. The fun_exact is the function of the exact solution. 
			and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
			whether the default quadrature is weighted quadrature. 
		"""
		u_tmp = np.atleast_2d(u_ctrlpts)
		nr    = np.size(u_tmp, axis=0)
		model = self.part

		nbqp, indices, basis, parweightsbyDir = [], [], [], []
		for i in range(model.dim):
			quadRule = GaussQuadrature(model.degree[i], model.knotvector[i], quadArgs={'type':'leg'})
			_, dersIndices, dersBasis, _ = quadRule.getQuadratureRulesInfo()
			indi, indj = dersIndices; parweights = quadRule._parametricWeights
			
			nbqp.append(quadRule.nbqp); indices.append(indi); indices.append(indj)
			basis.append(dersBasis); parweightsbyDir.append(parweights)

		inputs = [*nbqp, *indices, *basis]
		if model.dim == 2:
			Jqp = geophy.eval_jacobien_2d(*inputs, model.ctrlpts)
			detJ, _ = geophy.eval_inverse_det(Jqp)
			qpPhy = geophy.interpolate_meshgrid_2d(*inputs, model.ctrlpts)
			u_interp = geophy.interpolate_meshgrid_2d(*inputs, u_tmp)
		if model.dim == 3:
			Jqp = geophy.eval_jacobien_3d(*inputs, model.ctrlpts)
			detJ, _ = geophy.eval_inverse_det(Jqp)
			qpPhy = geophy.interpolate_meshgrid_3d(*inputs, model.ctrlpts)
			u_interp = geophy.interpolate_meshgrid_3d(*inputs, u_tmp)
		if nr == 1: u_interp = np.ravel(u_interp)

		u_exact     = exactfun(qpPhy)
		ue_diff_uh2 = (u_exact - u_interp)**2 
		ue2         = u_exact**2
		if nr > 1: 
			ue_diff_uh2 = np.ravel(np.sum(ue_diff_uh2, axis=0))
			ue2         = np.ravel(np.sum(ue2, axis=0))
		ue_diff_uh2 = ue_diff_uh2 * detJ
		ue2         = ue2*detJ

		ue_diff_uh2 = np.reshape(ue_diff_uh2, tuple(nbqp), order='F')
		ue2         = np.reshape(ue2, tuple(nbqp), order='F')
		if model.dim == 2: 
			tmp1 = np.einsum('i,j,ij->', parweightsbyDir[0], parweightsbyDir[1], ue2)
			tmp2 = np.einsum('i,j,ij->', parweightsbyDir[0], parweightsbyDir[1], ue_diff_uh2)
		if model.dim == 3: 
			tmp1 = np.einsum('i,j,k,ijk->', parweightsbyDir[0], parweightsbyDir[1], parweightsbyDir[2], ue2)
			tmp2 = np.einsum('i,j,k,ijk->', parweightsbyDir[0], parweightsbyDir[1], parweightsbyDir[2], ue_diff_uh2)
		error = np.sqrt(tmp2/tmp1)

		return error
	
	def L2projectionCtrlptsSurf(self, nbFacePosition, funfield=None, datafield=None):
		" Given the solution field (function or scattered points), it computes the L2 projection, ie. the value at control points. "

		return

	def L2projectionCtrlptsVol(self, volfun=None, datafield=None):
		" Given the solution field (function or scattered points), it computes the L2 projection, ie. the value at control points. "
		if datafield is not None: prop = datafield
		if volfun is not None:    prop = volfun(self.part.qpPhy)
		prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.detJ, prop]
		if self.part.dim == 2: vector = geophy.get_forcevol_2d(*inputs)
		if self.part.dim == 3: vector = geophy.get_forcevol_3d(*inputs)

		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, 
				*self.part.weights, self.part.invJ, self.part.detJ, vector, self._nbIterPCG, self._thresholdPCG]
		if self.part.dim == 2: u_interp, _ = geophy.l2projection_ctrlpts_2d(*inputs)
		if self.part.dim == 3: u_interp, _ = geophy.l2projection_ctrlpts_3d(*inputs)
		if nr == 1: u_interp = np.ravel(u_interp)
		return u_interp

class heatproblem(problem):
	def __init__(self, material:thermomat, part:part, boundary:boundaryCondition, solverArgs={}):
		super().__init__(part, boundary, solverArgs)
		self.material = material
		return
	
	def eval_mfConductivity(self, array_in, args=None):
		if args is None: args = self.part.qpPhy
		prop = self.material.conductivity(args)
		inpts = [*super()._getInputs(), self.part.invJ, self.part.detJ, prop]
		if self.part.dim == 2: array_out = heatsolver.mf_get_ku_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = heatsolver.mf_get_ku_3d(*inpts, array_in)
		return array_out
	
	def eval_mfCapacity(self, array_in, args=None, isLumped=False): 
		if args is None: args = self.part.qpPhy
		prop = self.material.capacity(args)
		inpts = [*super()._getInputs(), isLumped, self.part.invJ, self.part.detJ, prop]
		if self.part.dim == 2: array_out = heatsolver.mf_get_cu_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = heatsolver.mf_get_cu_3d(*inpts, array_in)
		return array_out

	def interpolate_temperature(self, u_ctrlpts):
		" Computes the temperature at the quadrature points "
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, np.atleast_2d(u_ctrlpts)]
		if self.part.dim == 2:   uinterp = geophy.interpolate_meshgrid_2d(*inpts)
		elif self.part.dim == 3: uinterp = geophy.interpolate_meshgrid_3d(*inpts)
		uinterp = np.ravel(uinterp)
		return uinterp

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
		if self.part.dim == 2: sol, residue = heatsolver.solver_lineartransient_heat_2d(*inpts)
		if self.part.dim == 3: sol, residue = heatsolver.solver_lineartransient_heat_3d(*inpts)
		return sol, residue

	def solveNLTransientHeatProblemPy(self, Tinout, time_list, Fext_list, theta=1.0, isLumped=False):
		nbctrlpts_total = self.part.nbctrlpts_total
		nbSteps = len(time_list)
		dod = self.boundary.getThermalBoundaryConditionInfo()[0]

		VVn0 = np.zeros(nbctrlpts_total)
		if nbSteps == 2:
			dt = time_list[1] - time_list[0]
			VVn0[dod] = 1.0/dt*(Tinout[dod, 1] - Tinout[dod, 0])
		elif nbSteps > 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			VVn0[dod] = 1.0/(dt1*(factor - factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
		else: raise Warning('At least 2 steps')
		
		resPCG_list = []
		for i in range(1, nbSteps):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			TTn0 = np.copy(Tinout[:, i-1])

			# Get values of new step
			TTn1  = TTn0 + dt*(1-theta)*VVn0; TTn1[dod] = np.copy(Tinout[dod, i])
			TTn10 = np.copy(TTn1); VVn1 = np.zeros(nbctrlpts_total)
			VVn1[dod] = 1.0/theta*(1.0/dt*(Tinout[dod, i]-Tinout[dod, i-1]) - (1-theta)*VVn0[dod])
			Fstep = Fext_list[:, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR):

				# Compute temperature at each quadrature point
				TTinterp = self.interpolate_temperature(TTn1)
			
				# Compute internal force
				args = np.row_stack((self.part.qpPhy, TTinterp))
				CdT  = self.eval_mfCapacity(VVn1, args=args, isLumped=isLumped)
				KT   = self.eval_mfConductivity(TTn1, args=args)
				Fint = KT + CdT

				# Compute residue
				dF    = Fstep - Fint; dF[dod] = 0.0
				resNR = np.sqrt(np.dot(dF, dF))
				print('NR error: %.5e' %resNR)
				if resNR <= self._thresholdNR: break

				# Iterative solver
				resPCG = np.array([i, j+1])
				vtmp, resPCGt = self.solveLinearTransientHeatProblemFT(dF, theta*dt, args=args, isLumped=isLumped)
				resPCG = np.append(resPCG, resPCGt)
				resPCG_list.append(resPCG)

				# Update values
				VVn1 += vtmp
				TTn1 = TTn10 + theta*dt*VVn1

			Tinout[:, i] = np.copy(TTn1)
			VVn0 = np.copy(VVn1)

		return resPCG_list

class mechaproblem(problem):
	def __init__(self, material:mechamat, part:part, boundary:boundaryCondition, solverArgs={}):
		super().__init__(part, boundary, solverArgs)
		self.material = material
		return
	
	def eval_mfStiffness(self, array_in, mechArgs=None):
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
		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, self.part.invJ, displacement]
		if   self.part.dim == 2: strain = plasticitysolver.interpolate_strain_2d(*inputs)
		elif self.part.dim == 3: strain = plasticitysolver.interpolate_strain_3d(*inputs)
		return strain
	
	def compute_intForce(self, stress):
		"Compute internal force using sigma coefficients "
		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.invJ, self.part.detJ, stress]
		if   self.part.dim == 2: intForce = plasticitysolver.get_intforce_2d(*inputs)
		elif self.part.dim == 3: intForce = plasticitysolver.get_intforce_3d(*inputs)
		return intForce

	def solveElasticityProblemFT(self, Fext, mechArgs=None):
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
		inputs = [*super()._getInputs(), *dod, self.boundary.mchDirichletTable, 
				self.part.invJ, self.part.detJ, prop, mechArgs, Fext, self._nbIterPCG, self._thresholdPCG, self._methodPCG]
		if   self.part.dim == 2: displacement, residue = plasticitysolver.solver_elasticity_2d(*inputs)
		elif self.part.dim == 3: displacement, residue = plasticitysolver.solver_elasticity_3d(*inputs)

		strain = self.interpolate_strain(displacement)
		stress = self.material.evalElasticStress(strain)
		
		return displacement, residue, stress
	
	def solvePlasticityProblemPy(self, Fext_list): 

		if not self.material._isPlasticityPossible: raise Warning('Plasticity not defined')
		dimen  = self.part.dim
		nvoigt = int(dimen*(dimen+1)/2)
		nbqp_total = self.part.nbqp_total

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
		disp = np.zeros(np.shape(Fext_list))
		stress_r = np.zeros((nvoigt, nbqp_total, np.shape(Fext_list)[2]))
		resPCG_list = []

		for i in range(1, np.shape(Fext_list)[2]):
			
			# Get values of last step
			ddisp = np.zeros(np.shape(disp[:, :, i-1]))

			# Get values of new step
			Fstep = Fext_list[:, :, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR):

				# Compute strain at each quadrature point
				d_n1 = disp[:, :, i-1] + ddisp
				strain = self.interpolate_strain(d_n1)
	
				# Closest point projection in perfect plasticity
				output = self.material.returnMappingAlgorithm(strain, pls_n0, a_n0, b_n0)
				stress, pls_n1, a_n1 = output[:nvoigt, :], output[nvoigt:2*nvoigt, :], output[2*nvoigt, :]
				b_n1, mechArgs = output[2*nvoigt+1:3*nvoigt+1, :], output[3*nvoigt+1:, :]

				# Compute internal force 
				Fint = self.compute_intForce(stress)
				
				# Compute residue
				dF   = Fstep - Fint
				clean_dirichlet(dF, self.boundary.mchdod) 
				resNR = np.sqrt(block_dot_product(dimen, dF, dF))
				print('NR error: %.5e' %resNR)
				if resNR <= self._thresholdNR: break
				
				# Iterative solver
				resPCG = np.array([i, j+1])
				vtmp, resPCGt = self.solveElasticityProblemFT(Fext=dF, mechArgs=mechArgs)
				resPCG = np.append(resPCG, resPCGt)
				resPCG_list.append(resPCG)

				ddisp += vtmp

			disp[:, :, i] = d_n1
			stress_r[:, :, i] = stress			
			pls_n0 = np.copy(pls_n1)
			a_n0 = np.copy(a_n1)
			b_n0 = np.copy(b_n1)

		return disp, resPCG_list, stress_r

