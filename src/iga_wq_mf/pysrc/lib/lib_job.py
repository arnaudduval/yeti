from .__init__ import *
from .lib_base import get_faceInfo, get_INCTable, evalDersBasisFortran
from .lib_quadrules import GaussQuadrature
from .lib_material import (heatmat, mechamat, clean_dirichlet, block_dot_product)
from .lib_part import part
from .lib_boundary import boundaryCondition

class problem():
	def __init__(self, part:part, boundary:boundaryCondition, solverArgs:dict):
		self.part     = part
		self.boundary = boundary
		self.addSolverConstraints(solverArgs)
		return
	
	def _getInputs(self):
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, *self.part.weights]
		return inpts

	def addSolverConstraints(self, solverArgs:dict):
		self._linPreCond = solverArgs.get('preconditioner', 'JMC')
		self._itersNL = solverArgs.get('iters_nonlinear', 20)
		self._thresNL = solverArgs.get('thres_nonlinear', 1e-8)
		self._itersLin = solverArgs.get('iters_linear', 100)
		self._thresLin = solverArgs.get('thres_linear', 1e-8)
		self._safeguard = 1e-14
		return
	
	def __getInfo4surfForce(self, nbFacePosition):
		INC_ctrlpts = get_INCTable(self.part.nbctrlpts)
		direction, side = get_faceInfo(nbFacePosition)
		assert direction < 2*self.part.dim, 'Not possible'
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
			The surffun is a Neumann like function, ie, in transfer heat q = -(k grad(T)).normal
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
	
	def normOfError(self, u_ctrlpts, normArgs:dict):
		""" Computes the norm L2 or H1 of the error. The exactfun is the function of the exact solution. 
			and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
			whether the default quadrature is weighted quadrature. 
		"""	
		typeNorm = normArgs.get('type', 'l2').lower()
		if all(norm != typeNorm  for norm in ['l2', 'h1', 'semih1']): raise Warning('Unknown norm')

		# Compute u interpolation
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
			detJ, invJ = geophy.eval_inverse_det(Jqp)
			qpPhy = geophy.interpolate_meshgrid_2d(*inpts, self.part.ctrlpts)
			u_interp = geophy.interpolate_meshgrid_2d(*inpts, np.atleast_2d(u_ctrlpts))
			derstemp = geophy.eval_jacobien_2d(*inpts, np.atleast_2d(u_ctrlpts))

		elif self.part.dim == 3:
			Jqp = geophy.eval_jacobien_3d(*inpts, self.part.ctrlpts)
			detJ, invJ = geophy.eval_inverse_det(Jqp)
			qpPhy = geophy.interpolate_meshgrid_3d(*inpts, self.part.ctrlpts)
			u_interp = geophy.interpolate_meshgrid_3d(*inpts, np.atleast_2d(u_ctrlpts))
			derstemp = geophy.eval_jacobien_3d(*inpts, np.atleast_2d(u_ctrlpts))
	
		u_interp = np.atleast_2d(u_interp); uders_interp = None
		uders_interp = np.atleast_3d(np.einsum('ijl,jkl->ikl', derstemp, invJ))

		# Compute u exact
		u_exact, uders_exact = None, None
		exactfun = normArgs.get('exactFunction', None)
		exactfunders = normArgs.get('exactFunctionDers', None)
		if callable(exactfun): u_exact = np.atleast_2d(exactfun(qpPhy))
		if callable(exactfunders): uders_exact = np.atleast_3d(exactfunders(qpPhy))

		part_ref = normArgs.get('part_ref', None); u_ref = normArgs.get('u_ref', None)
		if isinstance(part_ref, part) and isinstance(u_ref, np.ndarray):
			basisExact, indicesExact = [], []
			for i in range(self.part.dim):
				basis, indi, indj = evalDersBasisFortran(part_ref.degree[i], part_ref.knotvector[i], quadPts[i])
				basisExact.append(basis); indicesExact.append(indi); indicesExact.append(indj)
			inpts = [*nbqp, *indicesExact, *basisExact]
			if self.part.dim == 2:   
				u_exact = geophy.interpolate_meshgrid_2d(*inpts, np.atleast_2d(u_ref))    
				JqpExact = geophy.eval_jacobien_2d(*inpts, part_ref.ctrlpts)
				_, invJExact = geophy.eval_inverse_det(JqpExact) 
				derstemp = geophy.eval_jacobien_2d(*inpts, np.atleast_2d(u_ref))
			elif self.part.dim == 3: 
				u_exact = geophy.interpolate_meshgrid_3d(*inpts, np.atleast_2d(u_ref))
				JqpExact = geophy.eval_jacobien_3d(*inpts, part_ref.ctrlpts)
				_, invJExact = geophy.eval_inverse_det(JqpExact)
				derstemp = geophy.eval_jacobien_3d(*inpts, np.atleast_2d(u_ref))

			u_exact = np.atleast_2d(u_exact)
			uders_exact = np.atleast_3d(np.einsum('ijl,jkl->ikl', derstemp, invJExact))

		# Compute error
		uedfuf2_l2, uedfuf2_sh1 = 0., 0.
		ue2_l2, ue2_sh1 = 0., 0.

		if typeNorm == 'l2' or typeNorm == 'h1':
			uedfuf2_l2 += np.einsum('il->l', (u_exact - u_interp)**2)
			ue2_l2     += np.einsum('il->l', u_exact**2)

		if typeNorm == 'h1' or typeNorm == 'semih1':
			uedfuf2_sh1 += np.einsum('ijl->l', (uders_exact - uders_interp)**2)
			ue2_sh1     += np.einsum('ijl->l', uders_exact**2)

		norm1 = (uedfuf2_l2 + uedfuf2_sh1)*detJ
		norm2 = (ue2_l2 + ue2_sh1)*detJ

		norm1 = np.reshape(norm1, tuple(nbqp), order='F')
		norm2 = np.reshape(norm2, tuple(nbqp), order='F')
		if self.part.dim == 2: 
			tmp1 = np.einsum('i,j,ij->', parametricWeights[0], parametricWeights[1], norm1)
			tmp2 = np.einsum('i,j,ij->', parametricWeights[0], parametricWeights[1], norm2)
		if self.part.dim == 3: 
			tmp1 = np.einsum('i,j,k,ijk->', parametricWeights[0], parametricWeights[1], parametricWeights[2], norm1)
			tmp2 = np.einsum('i,j,k,ijk->', parametricWeights[0], parametricWeights[1], parametricWeights[2], norm2)
			
		abserror = np.sqrt(tmp1)
		relerror = np.sqrt(tmp1/tmp2)

		return abserror, relerror

	def L2projectionCtrlpts(self, u_atqp):
		" Given the solution field (function) over a physical space, it computes the L2 projection, ie. the value at control points. "
		prop = np.atleast_2d(u_atqp); nr = np.size(prop, axis=0)
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.detJ, prop]
		if self.part.dim == 2: volForce = geophy.get_forcevol_2d(*inpts)
		if self.part.dim == 3: volForce = geophy.get_forcevol_3d(*inpts)

		volForce = np.atleast_2d(volForce)
		inpts = [*self._getInputs(), self.part.detJ, volForce, self._itersLin, self._thresLin]
		if self.part.dim == 2: u_interp, _ = geophy.l2projection_ctrlpts_2d(*inpts)
		if self.part.dim == 3: u_interp, _ = geophy.l2projection_ctrlpts_3d(*inpts)
		if nr == 1: u_interp = np.ravel(u_interp)
		return u_interp

	def fastDiagonalization(self, array_in, fdtype='heat'):
		if fdtype == 'heat':
			inpts = [*self._getInputs(), array_in]
			if self.part.dim == 2: array_out = eigensolver.fastdiagonalization_2d(*inpts)
			if self.part.dim == 3: array_out = eigensolver.fastdiagonalization_3d(*inpts)
		else: raise Warning('Not coded')
		return array_out
		
class heatproblem(problem):
	def __init__(self, heat_material:heatmat, part:part, boundary:boundaryCondition, solverArgs={}):
		problem.__init__(self, part, boundary, solverArgs)
		self.heatmaterial = heat_material
		if self.heatmaterial.density is None: self.heatmaterial.addDensity(inpt=1.0, isIsotropic=True)
		return
	
	def compute_mfConductivity(self, array_in, args=None):
		if args is None: args = self.part.qpPhy
		prop = self.heatmaterial.conductivity(args)
		inpts = [*self._getInputs(), self.part.invJ, self.part.detJ, prop]
		if self.part.dim == 2: array_out = heatsolver.mf_conductivity_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = heatsolver.mf_conductivity_3d(*inpts, array_in)
		return array_out
	
	def compute_mfCapacity(self, array_in, args=None, isLumped=False): 
		if args is None: args = self.part.qpPhy
		prop = self.heatmaterial.capacity(args)*self.heatmaterial.density(args)
		inpts = [*self._getInputs(), isLumped, self.part.invJ, self.part.detJ, prop]
		if self.part.dim == 2: array_out = heatsolver.mf_capacity_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = heatsolver.mf_capacity_3d(*inpts, array_in)
		return array_out

	def interpolate_temperature(self, T_ctrlpts):
		" Computes the temperature at the quadrature points "
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, np.atleast_2d(T_ctrlpts)]
		if self.part.dim == 2:   T_interp = geophy.interpolate_meshgrid_2d(*inpts)
		elif self.part.dim == 3: T_interp = geophy.interpolate_meshgrid_3d(*inpts)
		T_interp = np.ravel(T_interp)
		return T_interp

	def compute_HeatIntForce(self, temp, flux, args=None, isLumped=False):
		if args is None: args = self.part.qpPhy
		intForce = self.compute_mfCapacity(flux, args=args, isLumped=isLumped) + self.compute_mfConductivity(temp, args=args)
		return intForce
	
	def compute_eigs_LOBPCG(self, ishigher=False, args=None):
		if args is None: args = self.part.qpPhy
		Cprop = self.heatmaterial.capacity(args)*self.heatmaterial.density(args)
		Kprop = self.heatmaterial.conductivity(args)
		inpts = [*self._getInputs(), self.boundary.thDirichletTable, self.part.invJ, self.part.detJ, 
				Cprop, Kprop, ishigher, self._itersLin, self._thresLin]
		if self.part.dim == 2: eigenval, eigenvec = eigensolver.solver_lobpcg_heat_2d(*inpts)
		if self.part.dim == 3: eigenval, eigenvec = eigensolver.solver_lobpcg_heat_3d(*inpts)
		return eigenval, eigenvec

	def solveSteadyHeatProblem(self, Fext, args=None):
		if args is None: args = self.part.qpPhy
		prop = self.heatmaterial.conductivity(args)
		inpts = [*self._getInputs(), self.boundary.thDirichletTable, self.part.invJ, self.part.detJ, 
				prop, Fext,  self._itersLin, self._thresLin, self._linPreCond]
		if self.part.dim == 2: temperature, residue = heatsolver.solver_linearsteady_heat_2d(*inpts)
		if self.part.dim == 3: temperature, residue = heatsolver.solver_linearsteady_heat_3d(*inpts)
		return temperature, residue
	
	def _solveLinearizedTransientProblem(self, Fext, tsfactor, args=None, isLumped=False):
		if args is None: args = self.part.qpPhy
		Cprop = self.heatmaterial.capacity(args)*self.heatmaterial.density(args)
		Kprop = self.heatmaterial.conductivity(args)
		inpts = [*self._getInputs(), isLumped, self.boundary.thDirichletTable, self.part.invJ, self.part.detJ,
				Cprop, Kprop, tsfactor, Fext, self._itersLin, self._thresLin, self._linPreCond]
		if self.part.dim == 2: temperature, residue = heatsolver.solver_lineartransient_heat_2d(*inpts)
		if self.part.dim == 3: temperature, residue = heatsolver.solver_lineartransient_heat_3d(*inpts)
		return temperature, residue

	def solveFourierTransientProblem(self, Tinout, Fext_list, time_list, alpha=1.0, isLumped=False):
		nbctrlpts_total = self.part.nbctrlpts_total; nsteps = len(time_list)
		dod = self.boundary.getThermalBoundaryConditionInfo()[0]

		# Compute inital velocity using interpolation
		assert nsteps > 2, 'At least 2 steps'
		V_n0 = np.zeros(nbctrlpts_total)
		dt1 = time_list[1] - time_list[0]
		dt2 = time_list[2] - time_list[0]
		factor = dt2/dt1
		V_n0[dod] = 1.0/(dt1*(factor - factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
		# ADD A PROCESS TO COMPUTE THE VELOCITY IN DOF

		AllresLin = []
		for i in range(1, nsteps):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			d_n0 = np.copy(Tinout[:, i-1])

			# Predict values of new step
			dj_n1 = d_n0 + (1 - alpha)*dt*V_n0
			Vj_n1 = np.zeros(self.part.nbctrlpts_total)
			
			# Overwrite inactive control points
			Vj_n1[dod] = 1.0/(alpha*dt)*(Tinout[dod, i] - dj_n1[dod])
			dj_n1[dod] = Tinout[dod, i]

			Fext_n1 = np.copy(Fext_list[:, i])

			print('Step: %d' %i)
			for j in range(self._itersNL):

				# Compute temperature at each quadrature point
				temperature = self.interpolate_temperature(dj_n1)
			
				# Compute internal force
				args = np.row_stack((self.part.qpPhy, temperature))
				Fint_dj = self.compute_HeatIntForce(dj_n1, Vj_n1, args=args, isLumped=isLumped)

				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				r_dj[dod] = 0.0

				resNLj = np.sqrt(np.dot(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('NonLinear error: %.5e' %resNLj)
				if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break

				# Solve for active control points
				resLinj = np.array([i, j+1])
				deltaV, resLin = self._solveLinearizedTransientProblem(r_dj, alpha*dt, args=args, isLumped=isLumped)
				resLinj = np.append(resLinj, resLin)
				
				# Update active control points
				dj_n1 += alpha*dt*deltaV
				Vj_n1 += deltaV
				AllresLin.append(resLinj)

			Tinout[:, i] = np.copy(dj_n1)
			V_n0 = np.copy(Vj_n1)

		return AllresLin

class mechaproblem(problem):
	def __init__(self, mechanical_material:mechamat, part:part, boundary:boundaryCondition, solverArgs={}):
		problem.__init__(self, part, boundary, solverArgs)
		self.mechamaterial = mechanical_material
		if self.mechamaterial.density is None: self.mechamaterial.addDensity(inpt=1.0, isIsotropic=True)
		return
	
	def compute_mfMass(self, array_in, args=None, isLumped=False):
		if args is None: args = self.part.qpPhy
		prop = self.mechamaterial.density(args)
		inpts = [*self._getInputs(), isLumped, self.part.invJ, self.part.detJ, prop]
		if self.part.dim == 2: array_out = plasticitysolver.mf_mass_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = plasticitysolver.mf_mass_3d(*inpts, array_in)
		return array_out
	
	def compute_mfStiffness(self, array_in, mechArgs=None):
		if mechArgs is None: 
			dimen = self.part.dim; nvoigt = int(dimen*(dimen+1)/2)
			mechArgs = np.zeros((2, self.part.nbqp_total))
			mechArgs[0, :] = self.mechamaterial.lame_lambda
			mechArgs[1, :] = self.mechamaterial.lame_mu
		inpts = [*self._getInputs(), self.part.invJ, self.part.detJ, mechArgs]
		if   self.part.dim == 2: array_out = plasticitysolver.mf_stiffness_2d(*inpts, array_in)
		elif self.part.dim == 3: array_out = plasticitysolver.mf_stiffness_3d(*inpts, array_in)
		return array_out
	
	def interpolate_strain(self, displacement):
		" Compute strain field from displacement field "
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.basis, self.part.invJ, displacement]
		if   self.part.dim == 2: strain = plasticitysolver.interpolate_strain_2d(*inpts)
		elif self.part.dim == 3: strain = plasticitysolver.interpolate_strain_3d(*inpts)
		return strain
	
	def compute_MechStaticIntForce(self, stress):
		"Compute internal force using sigma coefficients "
		inpts = [*self.part.nbqp[:self.part.dim], *self.part.indices, *self.part.weights, self.part.invJ, self.part.detJ, stress]
		if   self.part.dim == 2: intForce = plasticitysolver.get_intforce_2d(*inpts)
		elif self.part.dim == 3: intForce = plasticitysolver.get_intforce_3d(*inpts)
		return intForce

	def compute_MechDynamicIntForce(self, stress, accel, args=None, isLumped=False):
		if args is None: args = self.part.qpPhy
		intForce = self.compute_MechStaticIntForce(stress) + self.compute_mfMass(accel, args=args, isLumped=isLumped)
		return intForce
	
	def compute_eigs_LOBPCG(self, ishigher=False, mechArgs=None, args=None):
		if mechArgs is None:
			mechArgs = np.zeros((2, self.part.nbqp_total))
			mechArgs[0, :] = self.mechamaterial.lame_lambda
			mechArgs[1, :] = self.mechamaterial.lame_mu
		if args is None: args = self.part.qpPhy
		massProp = self.mechamaterial.density(args)
		inpts = [*self._getInputs(), self.boundary.mchDirichletTable, self.part.invJ, self.part.detJ, 
				mechArgs, massProp, ishigher, self._itersLin, self._thresLin]
		if self.part.dim == 2: eigenval, eigenvec = eigensolver.solver_lobpcg_elasticity_2d(*inpts)
		if self.part.dim == 3: eigenval, eigenvec = eigensolver.solver_lobpcg_elasticity_3d(*inpts)
		return eigenval, eigenvec

	def solveElasticityProblem(self, Fext, mechArgs=None):
		if mechArgs is None:
			mechArgs = np.zeros((2, self.part.nbqp_total))
			mechArgs[0, :] = self.mechamaterial.lame_lambda
			mechArgs[1, :] = self.mechamaterial.lame_mu
		inpts = [*self._getInputs(), self.boundary.mchDirichletTable, self.part.invJ, self.part.detJ, 
				mechArgs, Fext, self._itersLin, self._thresLin, self._linPreCond]
		if   self.part.dim == 2: displacement, residual = plasticitysolver.solver_linearelasticity_2d(*inpts)
		elif self.part.dim == 3: displacement, residual = plasticitysolver.solver_linearelasticity_3d(*inpts)
		return displacement, residual
	
	def solvePlasticityProblem(self, dispinout, Fext_list): 

		dimen  = self.part.dim
		nvoigt = int(dimen*(dimen+1)/2)
		nbChaboche = self.mechamaterial._chabocheNBparameters
		nbqp_total = self.part.nbqp_total
		nsteps = np.shape(Fext_list)[2]

		# Internal variables
		pls_n0 = np.zeros((nvoigt, nbqp_total))
		a_n0   = np.zeros((1, nbqp_total))
		b_n0   = np.zeros((nbChaboche, nvoigt, nbqp_total))
		pls_n1 = np.zeros((nvoigt, nbqp_total))
		a_n1   = np.zeros((1, nbqp_total))
		b_n1   = np.zeros((nbChaboche, nvoigt, nbqp_total))
		
		# Output variables
		Allstress  	 = np.zeros((nvoigt, nbqp_total, nsteps))
		Allstrain 	 = np.zeros((nvoigt, nbqp_total, nsteps))
		Allhardening = np.zeros((1, nbqp_total, nsteps))
		AllresLin 	 = []

		for i in range(1, nsteps):
			
			# Get values of last step
			d_n0 = np.copy(dispinout[:, :, i-1])

			# Predict values of new step
			dj_n1 = np.copy(d_n0)

			# Overwrite inactive control points
			for k in range(self.part.dim):
				dod = self.boundary.mchdod[k]
				dj_n1[k, dod] = dispinout[k, dod, i]
			
			Fext_n1 = np.copy(Fext_list[:, :, i])

			print('Step: %d' %i)
			for j in range(self._itersNL):

				# Compute strain at each quadrature point
				strain = self.interpolate_strain(dj_n1)
	
				# Closest point projection in perfect plasticity
				output, isElasticLoad = self.mechamaterial.J2returnMappingAlgorithm3D(strain, pls_n0, a_n0, b_n0)
				stress = output['stress']; pls_n1 = output['pls']; a_n1 = output['alpha']
				b_n1 = output['beta']; mechArgs = output['mechArgs']

				# Compute internal force 
				Fint_dj = self.compute_MechStaticIntForce(stress)
				
				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				clean_dirichlet(r_dj, self.boundary.mchdod) 

				resNLj = np.sqrt(block_dot_product(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('NonLinear error: %.5e' %resNLj)
				if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break
				if j > 0 and isElasticLoad: break
				
				# Solver for active control points
				resLinj = np.array([i, j+1])
				deltaD, resLin = self.solveElasticityProblem(Fext=r_dj, mechArgs=mechArgs)
				resLinj = np.append(resLinj, resLin)
				
				# Update active control points
				dj_n1 += deltaD
				AllresLin.append(resLinj)

			dispinout[:, :, i] = dj_n1
			Allstress[:, :, i] = stress	
			Allstrain[:, :, i] = strain
			Allhardening[0, :, i] = a_n1[0, :]

			pls_n0, a_n0, b_n0 = np.copy(pls_n1), np.copy(a_n1), np.copy(b_n1)

		return AllresLin, {'stress': Allstress, 'totalstrain': Allstrain, 'hardening':Allhardening}

	def _solveLinearizedElastoDynamicProblem(self, Fext, tsfactor, mechArgs=None, args=None, isLumped=False):
		if mechArgs is None:
			mechArgs = np.zeros((2, self.part.nbqp_total))
			mechArgs[0, :] = self.mechamaterial.lame_lambda
			mechArgs[1, :] = self.mechamaterial.lame_mu
		if args is None: args = self.part.qpPhy
		massProp = self.mechamaterial.density(args)
		inpts = [*self._getInputs(), isLumped, self.boundary.mchDirichletTable, 
				self.part.invJ, self.part.detJ, mechArgs, massProp, tsfactor,
				Fext, self._itersLin, self._thresLin, self._linPreCond]
		if   self.part.dim == 2: displacement, residual = plasticitysolver.solver_lineardynamics_2d(*inpts)
		elif self.part.dim == 3: displacement, residual = plasticitysolver.solver_lineardynamics_3d(*inpts)
		return displacement, residual

	def solveDynamicsProblem(self, dispinout, Fext_list, time_list, beta=0.25, gamma=0.5, isLumped=False):
		nbctrlpts_total = self.part.nbctrlpts_total; nsteps = len(time_list)
		V_n0 = np.zeros((self.part.dim, nbctrlpts_total))
		A_n0 = np.zeros((self.part.dim, nbctrlpts_total))

		# Compute initial velocity using interpolation
		if np.shape(dispinout)[2] > 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			for k in range(self.part.dim):
				dod = self.boundary.mchdod[k]
				V_n0[k, dod] = 1.0/(dt1*(factor-factor**2))*(dispinout[k, dod, 2] - (factor**2)*dispinout[k, dod, 1] - (1 - factor**2)*dispinout[k, dod, 0])
				A_n0[k, dod] = 2.0/(dt1*dt2)*((dispinout[k, dod, 2] - factor*dispinout[k, dod, 1])/(factor - 1) + dispinout[k, dod, 0])
		else: raise Warning('We need more than 2 steps')

		AllresLin = []
		for i in range(1, nsteps):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			d_n0 = np.copy(dispinout[:, :, i-1])

			# Predict values of new step
			dj_n1 = d_n0 + dt*V_n0 + 0.5*dt**2*(1 - 2*beta)*A_n0
			Vj_n1 = V_n0 + (1 - gamma)*dt*A_n0
			Aj_n1 = np.zeros((self.part.dim, nbctrlpts_total))
			
			# Overwrite inactive control points
			for k in range(self.part.dim):
				dod = self.boundary.mchdod[k]
				Aj_n1[k, dod] = (dispinout[k, dod, i] - dj_n1[k, dod])/(beta*dt**2)
				Vj_n1[k, dod] = Vj_n1[k, dod] + gamma*dt*Aj_n1[k, dod]
				dj_n1[k, dod] = dispinout[k, dod, i]

			Fext_n1 = np.copy(Fext_list[:, :, i])

			print('Step: %d' %i)
			for j in range(self._itersNL):

				# Compute strain and stress at each quadrature point
				strain = self.interpolate_strain(dj_n1)
				stress = self.mechamaterial.evalElasticStress(strain, self.part.dim)

				# Compute internal force 
				args = self.part.qpPhy
				Fint_dj = self.compute_MechDynamicIntForce(stress, Aj_n1, args=args, isLumped=isLumped)
				
				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				clean_dirichlet(r_dj, self.boundary.mchdod) 

				resNLj = np.sqrt(block_dot_product(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('NonLinear error: %.5e' %resNLj)
				if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break

				# Solve for active control points
				resLinj = np.array([i, j+1])
				deltaA, resLin = self._solveLinearizedElastoDynamicProblem(Fext=r_dj, tsfactor=beta*dt**2, isLumped=isLumped)
				resLinj = np.append(resLinj, resLin)
				
				# Update active control points
				dj_n1 += beta*dt**2*deltaA
				Vj_n1 += gamma*dt*deltaA
				Aj_n1 += deltaA
				AllresLin.append(resLinj)

			dispinout[:, :, i] = np.copy(dj_n1)
			V_n0 = np.copy(Vj_n1)
			A_n0 = np.copy(Aj_n1)

		return AllresLin

class thermomechaproblem(heatproblem, mechaproblem):

	def __init__(self, heatmaterial:heatmat, mechamaterial:mechamat, part:part, boundary:boundaryCondition, solverArgs={}):
		heatproblem.__init__(self, heat_material=heatmaterial, part=part, 
				boundary=boundary, solverArgs=solverArgs)
		mechaproblem.__init__(self, mechanical_material=mechamaterial, part=part, 
				boundary=boundary, solverArgs=solverArgs)
		return
	
	def addDensity(self, inpt, isIsotropic):
		self.mechamaterial.density = self.mechamaterial.setScalarProperty(inpt, isIsotropic=isIsotropic)
		self.heatmaterial.density  = self.heatmaterial.setScalarProperty(inpt, isIsotropic=isIsotropic)
		return

	def compute_mfCoupled(self, array_in, args=None, isThermal=True):
		if args is None: args = self.part.qpPhy
		prop = 3*self.mechamaterial.thermalExpansion*self.mechamaterial.lame_bulk*np.ones(self.part.nbqp_total)
		inpts = [*self._getInputs(), self.part.invJ, self.part.detJ, prop]
		if isThermal:
			if self.part.dim == 2: array_out = heatsolver.mf_thmchcoupled_2d(*inpts, array_in)
			if self.part.dim == 3: array_out = heatsolver.mf_thmchcoupled_3d(*inpts, array_in)
		else:
			if self.part.dim == 2: array_out = plasticitysolver.mf_thmchcoupled_2d(*inpts, array_in)
			if self.part.dim == 3: array_out = plasticitysolver.mf_thmchcoupled_3d(*inpts, array_in)
		return array_out
	
	def compute_ThermomechIntForce(self, disp, vel, accel, stress, args=None, isLumped=False):
		if args is None: args = self.part.qpPhy
		intForce = np.zeros((self.part.dim+1, self.part.nbctrlpts_total))
		intForce[:-1, :] = (self.compute_MechDynamicIntForce(stress, accel[:-1, :], args=args, isLumped=isLumped) 
							- self.compute_mfCoupled(disp[-1, :], args=args, isThermal=False))
		intForce[-1, :]  = (self.compute_HeatIntForce(disp[-1, :], vel[-1, :], args=args, isLumped=isLumped)
							+ self.heatmaterial.refTemp*self.compute_mfCoupled(vel[:-1, :], args=args, isThermal=True))
		return intForce

	def _solveLinearizedThermoElasticityProblem(self, Fext, tsfactor1, tsfactor2, mechArgs=None, args=None, isLumped=False):
		displacement = np.zeros((self.part.dim+1, self.part.nbctrlpts_total))
		displacement[:-1, :] = self._solveLinearizedElastoDynamicProblem(Fext[:-1, :], tsfactor2, mechArgs=mechArgs, args=args, isLumped=isLumped)[0]
		displacement[-1, :] = self._solveLinearizedTransientProblem(Fext[-1, :], tsfactor2/tsfactor1, args=args, isLumped=isLumped)[0]/tsfactor1
		return displacement
	
	def solveThermoElasticityProblem(self, dispinout, Tinout, Fmech_list, Fheat_list, time_list, beta=0.25, gamma=0.5, isLumped=False):
		nbctrlpts_total = self.part.nbctrlpts_total; nsteps = len(time_list)
		d_n0 = np.zeros((self.part.dim+1, nbctrlpts_total))
		V_n0 = np.zeros((self.part.dim+1, nbctrlpts_total))
		A_n0 = np.zeros((self.part.dim+1, nbctrlpts_total))

		# Compute initial velocity using interpolation
		if np.shape(dispinout)[2] > 2:
			dt1 = time_list[1] - time_list[0]
			dt2 = time_list[2] - time_list[0]
			factor = dt2/dt1
			for k in range(self.part.dim):
				dod = self.boundary.mchdod[k]
				V_n0[k, dod] = 1.0/(dt1*(factor-factor**2))*(dispinout[k, dod, 2] - (factor**2)*dispinout[k, dod, 1] - (1 - factor**2)*dispinout[k, dod, 0])
				A_n0[k, dod] = 2.0/(dt1*dt2)*((dispinout[k, dod, 2] - factor*dispinout[k, dod, 1])/(factor - 1) + dispinout[k, dod, 0])
			dod = self.boundary.thdod
			V_n0[-1, dod] = 1.0/(dt1*(factor-factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
			A_n0[-1, dod] = 2.0/(dt1*dt2)*((Tinout[dod, 2] - factor*Tinout[dod, 1])/(factor - 1) + Tinout[dod, 0])
		else: raise Warning('We need more than 2 steps')

		for i in range(1, nsteps):
			
			# Get delta time
			dt = time_list[i] - time_list[i-1]
			
			# Get values of last step
			d_n0[:-1, :] = np.copy(dispinout[:, :, i-1])
			d_n0[-1, :] = np.copy(Tinout[:, i-1])

			# Predict values of new step
			dj_n1 = d_n0 + dt*V_n0 + 0.5*dt**2*(1 - 2*beta)*A_n0
			Vj_n1 = V_n0 + (1 - gamma)*dt*A_n0
			Aj_n1 = np.zeros((self.part.dim+1, self.part.nbctrlpts_total))
			
			# Overwrite inactive control points
			for k in range(self.part.dim):
				dod = self.boundary.mchdod[k]
				Aj_n1[k, dod] = (dispinout[k, dod, i] - dj_n1[k, dod])/(beta*dt**2)
				Vj_n1[k, dod] = Vj_n1[k, dod] + gamma*dt*Aj_n1[k, dod]
				dj_n1[k, dod] = dispinout[k, dod, i]

			dod = self.boundary.getThermalBoundaryConditionInfo()[0]
			Aj_n1[-1, dod] = (Tinout[dod, i] - dj_n1[-1, dod])/(beta*dt**2)
			Vj_n1[-1, dod] = Vj_n1[-1, dod] + gamma*dt*Aj_n1[-1, dod]
			dj_n1[-1, dod] = Tinout[dod, i]
			
			Fext_n1 = np.zeros((self.part.dim+1, nbctrlpts_total))
			Fext_n1[:-1, :] = np.copy(Fmech_list[:, :, i])
			Fext_n1[-1, :] = np.copy(Fheat_list[:, i])

			print('Step: %d' %i)
			for j in range(self._itersNL):
				
				# Compute strain and stress at each quadrature point
				strain = self.interpolate_strain(dj_n1[:-1, :])
				stress = self.mechamaterial.evalElasticStress(strain, self.part.dim)
				temperature = self.interpolate_temperature(dj_n1[-1, :])

				# Compute internal force 
				args = np.row_stack((self.part.qpPhy, temperature))
				Fint_dj = self.compute_ThermomechIntForce(dj_n1, Vj_n1, Aj_n1, stress, 
														args=args, isLumped=isLumped)
				
				# Compute residue
				r_dj = Fext_n1 - Fint_dj
				dod = [*self.boundary.mchdod, self.boundary.thdod]
				clean_dirichlet(r_dj, dod) 

				resNLj = np.sqrt(block_dot_product(r_dj, r_dj))
				if j == 0: resNL0 = resNLj
				print('Nonlinear error: %.5e' %resNLj)
				if resNLj <= max([self._safeguard, self._thresNL*resNL0]): break

				# Solve for active control points
				deltaA = self._solveLinearizedThermoElasticityProblem(Fext=r_dj, tsfactor1=gamma*dt, tsfactor2=beta*dt**2, args=args, isLumped=isLumped)
				
				# Update active control points
				dj_n1 += beta*dt**2*deltaA
				Vj_n1 += gamma*dt*deltaA
				Aj_n1 += deltaA

			dispinout[:, :, i] = np.copy(dj_n1[:-1, :])
			Tinout[:, i] = np.copy(dj_n1[-1, :])
			V_n0 = np.copy(Vj_n1)
			A_n0 = np.copy(Aj_n1)
		
		return 