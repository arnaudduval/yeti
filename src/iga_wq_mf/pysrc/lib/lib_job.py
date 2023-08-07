from .__init__ import *
from .lib_base import eraseRowsCSR
from .lib_quadrules import GaussQuadrature
from .lib_material import (thermomat, mechamat, 
							clean_dirichlet, block_dot_product)
from .lib_part import part
from .lib_boundary import boundaryCondition, get_INCTable

class problem():
	def __init__(self, part:part, boundary:boundaryCondition, solverArgs:dict):
		self.material = None
		self.part     = part
		self.boundary = boundary
		self.addSolverConstraints(solverArgs)
		return
	
	def addSolverConstraints(self, solverArgs:dict):
		self._nbIterPCG    = solverArgs.get('nbIterationsPCG', 50)
		self._nbIterNR     = solverArgs.get('nbIterationsNR', 30)
		self._thresholdPCG = solverArgs.get('PCGThreshold', 1e-10)
		self._thresholdNR  = solverArgs.get('NRThreshold', 1e-8)
		self._methodPCG    = solverArgs.get('PCGmethod', 'JMC')
		return
	
	def L2NormOfError(self, fun_exact, u_ctrlpts):
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

		u_exact     = fun_exact(qpPhy)
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
	
	def L2projectionCtrlpts(self, funfield=None, datafield=None):
		" Given the solution field (function or scattered points), it computes the L2 projection, ie. the value at control points. "
		coefs = None
		if datafield is not None: coefs = datafield*self.part.detJ
		if funfield is not None:  coefs = funfield(self.part.qpPhy)*self.part.detJ
		if coefs is None: raise Warning('Missing data')
		coefs = np.atleast_2d(coefs); nr = np.size(coefs, axis=0); u_interp = []

		for i in range(nr):
			inputs = [coefs[i, :], *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights]
			if self.part.dim == 2: vector = heatsolver.get_heatvol_2d(*inputs)
			if self.part.dim == 3: vector = heatsolver.get_heatvol_3d(*inputs)

			inputs = [self.part.detJ, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, 
					*self.part.weights, vector, self._nbIterPCG, self._thresholdPCG]
			if self.part.dim == 2: u_tmp, _ = geophy.l2projection_ctrlpts_2d(*inputs)
			if self.part.dim == 3: u_tmp, _ = geophy.l2projection_ctrlpts_3d(*inputs)
			u_interp.append(u_tmp)

		return np.array(u_interp)

class heatproblem(problem):
	def __init__(self, material:thermomat, part:part, boundary:boundaryCondition, solverArgs={}):
		super().__init__(part, boundary, solverArgs)
		self.material = material
		return
	
	def get_input4MatrixFree(self, table=None):
		" Returns necessary inputs to compute the product between a matrix and a vector "
		
		if table is None: table = self.boundary.thDirichletTable
		indices, basis, weights = [], [], []
		for i in range(self.part.dim):
			# Select data
			if np.array_equal(table[i, :], [0, 0]): rows2erase = []
			if np.array_equal(table[i, :], [0, 1]): rows2erase = [-1]
			if np.array_equal(table[i, :], [1, 0]): rows2erase = [0]
			if np.array_equal(table[i, :], [1, 1]): rows2erase = [0, -1]
			indi_t, indj_t, data_t = eraseRowsCSR(rows2erase, 
									self.part.indices[2*i], self.part.indices[2*i+1],  
									[self.part.basis[i], self.part.weights[i]])
			
			# Extract data and append to list
			[basist, weightst] = data_t
			indices.append(indi_t); indices.append(indj_t) 
			basis.append(basist); weights.append(weightst)

		inputs = [*self.part.nbqp[:self.part.dim], *indices, *basis, *weights]

		return inputs
	
	def eval_mfConductivity(self, array_in, coefs=None, table=None, args=None):
		inputs = self.get_input4MatrixFree(table=table)
		if coefs is None: coefs  = self.material.eval_conductivityCoefficients(self.part.invJ, self.part.detJ, args)
		if self.part.dim == 2: array_out = heatsolver.mf_get_ku_2d(coefs, *inputs, array_in)
		if self.part.dim == 3: array_out = heatsolver.mf_get_ku_3d(coefs, *inputs, array_in)
		return array_out
	
	def eval_mfCapacity(self, u, coefs=None, table=None, args=None): 
		inputs = self.get_input4MatrixFree(table=table)
		if coefs is None: coefs  = self.material.eval_capacityCoefficients(self.part.detJ, args)
		if self.part.dim == 2: result = heatsolver.mf_get_cu_2d(coefs, *inputs, u)
		if self.part.dim == 3: result = heatsolver.mf_get_cu_3d(coefs, *inputs, u)
		return result

	def eval_volForce(self, fun, indi=None): 
		if indi is None: indi = np.arange(self.part.nbctrlpts_total, dtype=int)
		coefs = self.material.eval_heatForceCoefficients(fun, self.part.detJ, self.part.qpPhy)
		inputs = [coefs, *self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights]
		if self.part.dim == 2: vector = heatsolver.get_heatvol_2d(*inputs)[indi]
		if self.part.dim == 3: vector = heatsolver.get_heatvol_3d(*inputs)[indi]
		return vector
	
	def eval_surfForce(self, fun, nbFacePosition):
		if self.part.dim != 3: raise Warning('Method only for 3D geometries')

		def get_faceInfo(nb):
			direction = int(np.floor(nb/2))
			if nb%2 == 1: side = 1
			else: side = 0
			return direction, side

		vector = np.zeros(self.part.nbctrlpts_total)
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

		coefs = fun(qpPhy)
		inpts = [coefs, JJ, *nnz, *indices, *weights]
		if   self.part.dim == 2: raise Warning('Not done yet')
		elif self.part.dim == 3: tmp = plasticitysolver.get_heatsurf_3d(*inpts)
		vector[CPList] = tmp

		return vector

	def interpolate_temperature(self, uctrlpts):
		" Computes the temperature at the quadrature points "
		basis   = self.part.basis
		indices = self.part.indices
		nbqp    = self.part.nbqp[:self.part.dim]
		inputs = [*nbqp, *indices, *basis, np.atleast_2d(uctrlpts)]
		if self.part.dim == 2:   uinterp = geophy.interpolate_meshgrid_2d(*inputs)
		elif self.part.dim == 3: uinterp = geophy.interpolate_meshgrid_3d(*inputs)
		uinterp = np.ravel(uinterp)
		return uinterp

	def solveSteadyHeatProblemFT(self, b, coefs=None):
		if coefs is None: coefs  = self.material.eval_conductivityCoefficients(self.part.invJ, self.part.detJ, self.part.qpPhy)
		tmp = self.get_input4MatrixFree(table=self.boundary.thDirichletTable)
		inputs = [coefs, *tmp, b, self._nbIterPCG, self._thresholdPCG, self._methodPCG]
		if self.part.dim == 2: sol, residue = heatsolver.solver_steady_heat_2d(*inputs)
		if self.part.dim == 3: sol, residue = heatsolver.solver_steady_heat_3d(*inputs)
		return sol, residue
	
	def solveLinearTransientHeatProblemFT(self, dt, b, theta=1.0, Ccoefs=None, Kcoefs=None, args={}):
		if Ccoefs is None: 
			temperature = args.get('temperature')
			inpt = self.part.qpPhy
			inpt = np.row_stack((inpt, temperature*np.ones(np.size(inpt, axis=1))))
			Ccoefs  = self.material.eval_capacityCoefficients(self.part.detJ, inpt)
			
		if Kcoefs is None:
			temperature = args.get('temperature')
			inpt = self.part.qpPhy
			inpt = np.row_stack((inpt, temperature*np.ones(np.size(inpt, axis=1))))
			Kcoefs  = self.material.eval_conductivityCoefficients(self.part.invJ, self.part.detJ, inpt)
		
		tmp = self.get_input4MatrixFree(table=self.boundary.thDirichletTable)
		inputs = [Ccoefs, Kcoefs, *tmp, b, theta*dt, self._nbIterPCG, self._thresholdPCG, self._methodPCG]
		if self.part.dim == 2: raise Warning('Until now not done')
		if self.part.dim == 3: sol, residue = heatsolver.solver_lineartransient_heat_3d(*inputs)
		return sol, residue

	def solveNLTransientHeatProblemPy(self, Tinout, time_list, Fext, theta=1.0):
		nbSteps = len(time_list)
		nbctrlpts_total = self.part.nbctrlpts_total
		dod, _, dof = self.boundary.getThermalBoundaryConditionInfo()

		VVn0 = np.zeros(nbctrlpts_total)
		resPCG_list = []
		if nbSteps == 2:
			dt = time_list[1] - time_list[0]
			VVn0[dod] = 1.0/dt*(Tinout[dod, 1] - Tinout[dod, 0])
		elif nbSteps > 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			VVn0[dod] = 1.0/(dt1*(factor - factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
		else: raise Warning('At least 2 steps')
		
		for i in range(1, nbSteps):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			TTn0 = np.copy(Tinout[:, i-1])

			# Get values of new step
			TTn1  = TTn0 + dt*(1-theta)*VVn0; TTn1[dod] = np.copy(Tinout[dod, i])
			TTn10 = np.copy(TTn1); VVn1 = np.zeros(np.shape(VVn0))
			VVn1[dod] = 1.0/theta*(1.0/dt*(Tinout[dod, i]-Tinout[dod, i-1]) - (1-theta)*VVn0[dod])
			Fstep = Fext[:, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR):

				# Compute temperature at each quadrature point
				TTinterp = self.interpolate_temperature(TTn1)
			
				# Compute internal force
				inpt = np.row_stack((self.part.qpPhy, TTinterp))
				Ccoefs = self.material.eval_capacityCoefficients(self.part.detJ, inpt)
				Kcoefs = self.material.eval_conductivityCoefficients(self.part.invJ, self.part.detJ, inpt)

				CdTemp = self.eval_mfCapacity(VVn1, coefs=Ccoefs, table=np.zeros((3, 2), dtype=bool))
				KTemp  = self.eval_mfConductivity(TTn1, coefs=Kcoefs, table=np.zeros((3, 2), dtype=bool))
				Fint   = KTemp + CdTemp

				# Compute residue
				dF    = Fstep[dof] - Fint[dof]
				resNR = np.sqrt(np.dot(dF, dF))
				print('NR error: %.5e' %resNR)
				if resNR <= self._thresholdNR: break

				# Iterative solver
				resPCG = np.array([i, j+1])
				vtmp, resPCGt = self.solveLinearTransientHeatProblemFT(dt, dF, Ccoefs=Ccoefs, Kcoefs=Kcoefs, theta=theta)
				resPCG = np.append(resPCG, resPCGt)
				resPCG_list.append(resPCG)

				# Update values
				VVn1[dof] += vtmp
				TTn1[dof] = TTn10[dof] + theta*dt*VVn1[dof]

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
		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, 
	    			*self.part.basis, *self.part.weights, self.part.invJ, self.part.detJ, mechArgs]
		if   self.part.dim == 2: array_out = plasticitysolver.mf_get_su_2d(*inputs, array_in)
		elif self.part.dim == 3: array_out = plasticitysolver.mf_get_su_3d(*inputs, array_in)
		return array_out
	
	def compute_strain(self, displacement):
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
	
	def eval_volForce(self, volfun):
		prop   = volfun(self.part.qpPhy)
		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.detJ, prop]
		if   self.part.dim == 2: volForce = plasticitysolver.get_forcevol_2d(*inputs)
		elif self.part.dim == 3: volForce = plasticitysolver.get_forcevol_3d(*inputs)
		return volForce
	
	def eval_surfForce(self, surffun, nbFacePosition):
		" Returns force vector at the surface. In 3D: surface integrals. "

		def get_faceInfo(nb):
			direction = int(np.floor(nb/2))
			if nb%2 == 1: side = 1
			else: side = 0
			return direction, side

		surfForce = np.zeros((self.part.dim, self.part.nbctrlpts_total))
		INC_ctrlpts = get_INCTable(self.part.nbctrlpts)
		direction, side = get_faceInfo(nbFacePosition)
		
		if direction>=2*self.part.dim: raise Warning('Not possible')
		valrange = [i for i in range(self.part.dim)]
		valrange.pop(direction)

		# Get control points and quadrature points list
		if side == 0: CPList = np.where(INC_ctrlpts[:, direction] == 0)[0]
		elif side == 1: CPList = np.where(INC_ctrlpts[:, direction] == self.part.nbctrlpts[direction]-1)[0]
		
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

		# Get force values at quadrature points
		prop  = surffun(qpPhy)
		inpts = [*nnz, *indices, *weights, JJ, prop]
		if   self.part.dim == 2: tmp = plasticitysolver.get_forcesurf_2d(*inpts)
		elif self.part.dim == 3: tmp = plasticitysolver.get_forcesurf_3d(*inpts)
		surfForce[:, CPList] = tmp

		for i in range(self.part.dim): surfForce[i, self.boundary.mchdod[i]] = 0.0

		return surfForce

	def solveElasticityProblemFT(self, Fext, mechArgs=None):
		dod_total = deepcopy(self.boundary.mchdod)
		for i, dod in enumerate(dod_total):
			tmp = dod + 1; dod_total[i] = tmp

		dimen  = self.part.dim
		nvoigt = int(dimen*(dimen+1)/2)
		prop = [self.material.elasticmodulus, self.material.poissonratio, self.material.elasticlimit]
		if mechArgs is None:
			mechArgs = np.zeros((nvoigt+3, self.part.nbqp_total))
			mechArgs[0, :] = self.material.lame_lambda
			mechArgs[1, :] = self.material.lame_mu
		inputs = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, 
				*self.part.weights, Fext, *dod_total, self.boundary.mchDirichletTable, 
				self.part.invJ, self.part.detJ, prop, mechArgs, self._nbIterPCG, self._thresholdPCG, self._methodPCG]
		if   self.part.dim == 2: displacement, residue = plasticitysolver.solver_elasticity_2d(*inputs)
		elif self.part.dim == 3: displacement, residue = plasticitysolver.solver_elasticity_3d(*inputs)

		strain = self.compute_strain(displacement)
		stress = self.material.evalElasticStress(strain)
		
		return displacement, residue, stress
	
	def solvePlasticityProblemPy(self, Fext): 

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
		disp = np.zeros(np.shape(Fext))
		stress_r = np.zeros((nvoigt, nbqp_total, np.shape(Fext)[2]))
		resPCG_list = []

		for i in range(1, np.shape(Fext)[2]):
			
			# Get values of last step
			ddisp = np.zeros(np.shape(disp[:, :, i-1]))

			# Get values of new step
			Fstep = Fext[:, :, i]

			print('Step: %d' %i)
			for j in range(self._nbIterNR):

				# Compute strain at each quadrature point
				d_n1 = disp[:, :, i-1] + ddisp
				strain = self.compute_strain(d_n1)
	
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

