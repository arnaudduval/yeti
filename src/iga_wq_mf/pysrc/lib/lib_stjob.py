"SPACE TIME GALERKIN"
from .__init__ import *
from .lib_base import get_faceInfo, get_INCTable, evalDersBasisFortran
from .lib_quadrules import GaussQuadrature
from .lib_material import (heatmat)
from .lib_part import part, part1D
from .lib_boundary import boundaryCondition

class stproblem():
	def __init__(self, part:part, tspan:part1D, boundary:boundaryCondition, solverArgs:dict):
		self.part     = part
		self.time     = tspan
		self.boundary = boundary
		self.addSolverConstraints(solverArgs)
		return
	
	def _getInputs(self):
		inpts = [*self.part.nbqp[:self.part.dim], self.time.nbqp, *self.part.indices, *self.time.quadRule.dersIndices, 
		   		*self.part.basis, self.time.quadRule.dersBasis, *self.part.weights, self.time.quadRule.dersWeights]
		return inpts

	def addSolverConstraints(self, solverArgs:dict):
		self._linSolv = solverArgs.get('linearsolver', 'GMRES')
		self._itersLin = solverArgs.get('iters_linear', 100)
		self._thresLin = solverArgs.get('thres_linear', 1e-8)
		self._linPreCond = solverArgs.get('preconditioner', 'JMC')
		self._itersNL = solverArgs.get('iters_nonlinear', 15)
		self._thresNL = solverArgs.get('thres_nonlinear', 1e-8)
		self._safeguard = 1e-14
		return
	
	def compute_volForce(self, volfun, args=None): 
		" Computes the volume force over a geometry "
		if args is None: args={'position': self.part.qpPhy}
		prop = volfun(args)
		prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
		inpts = [*self.part.nbqp[:self.part.dim], self.time.nbqp, *self.part.indices, *self.time.quadRule.dersIndices, 
				*self.part.weights, self.time.quadRule.dersWeights, self.part.detJ, self.time.detJ, prop]
		if self.part.dim == 2: volForce = geophy.get_forcevol_st2d(*inpts)
		if self.part.dim == 3: volForce = geophy.get_forcevol_st3d(*inpts)
		if nr == 1: volForce = np.ravel(volForce)
		return volForce	

	def normOfError(self, u_ctrlpts, normArgs:dict):
		""" Computes the norm L2 or H1 of the error. The exactfun is the function of the exact solution. 
			and u_ctrlpts is the field at the control points. We compute the integral using Gauss Quadrature
			whether or not the default quadrature is weighted quadrature. 
		"""	
		typeNorm = normArgs.get('type', 'l2').lower()
		if all(norm != typeNorm  for norm in ['l2', 'h1', 'semih1']): raise Warning('Unknown norm')

		# Compute u interpolation in space
		nbqp, quadPts, indices, basis, parametricWeights = [], [], [], [], []
		for i in range(self.part.dim):
			quadRule = GaussQuadrature(self.part.degree[i], self.part.knotvector[i], quadArgs={'type':'leg'})
			quadPtsByDir, indicesByDir, basisByDir, _ = quadRule.getQuadratureRulesInfo()
			indi, indj = indicesByDir; parweightsByDir = quadRule._parametricWeights
			
			nbqp.append(quadRule.nbqp); quadPts.append(quadPtsByDir); indices.append(indi); indices.append(indj)
			basis.append(basisByDir); parametricWeights.append(parweightsByDir)

		# Compute u interpolation in time
		quadRule = GaussQuadrature(self.time.degree, self.time.knotvector, quadArgs={'type':'leg'})
		quadPtsByDir, indicesByDir, basisByDir, _ = quadRule.getQuadratureRulesInfo()
		indi, indj = indicesByDir; parweightsByDir = quadRule._parametricWeights
		nbqp.append(quadRule.nbqp); quadPts.append(quadPtsByDir); indices.append(indi); indices.append(indj)
		basis.append(basisByDir); parametricWeights.append(parweightsByDir)

		sptimectrlpts = np.zeros((self.part.dim+1, self.part.nbctrlpts_total*self.time.nbctrlpts_total))
		for i in range(self.time.nbctrlpts_total):
			iold = i*self.part.nbctrlpts_total; inew = (i + 1)*self.part.nbctrlpts_total
			sptimectrlpts[:, iold:inew] = np.vstack([self.part.ctrlpts, self.time.ctrlpts[i]*np.ones(self.part.nbctrlpts_total)])

		inpts = [*nbqp, *indices, *basis]
		if self.part.dim == 2:
			Jqp = geophy.eval_jacobien_3d(*inpts, sptimectrlpts)
			detJ, invJ = geophy.eval_inverse_det(Jqp)
			qpPhy = geophy.interpolate_meshgrid_3d(*inpts, sptimectrlpts)
			u_interp = geophy.interpolate_meshgrid_3d(*inpts, np.atleast_2d(u_ctrlpts))
			derstemp = geophy.eval_jacobien_3d(*inpts, np.atleast_2d(u_ctrlpts))
	
		elif self.part.dim == 3:
			Jqp = geophy.eval_jacobien_4d(*inpts, sptimectrlpts)
			detJ, invJ = geophy.eval_inverse_det(Jqp)
			qpPhy = geophy.interpolate_meshgrid_4d(*inpts, sptimectrlpts)
			u_interp = geophy.interpolate_meshgrid_4d(*inpts, np.atleast_2d(u_ctrlpts))
			derstemp = geophy.eval_jacobien_4d(*inpts, np.atleast_2d(u_ctrlpts))

		u_interp = np.atleast_2d(u_interp)
		uders_interp = np.atleast_3d(np.einsum('ijl,jkl->ikl', derstemp, invJ))

		# Compute u exact
		u_exact, uders_exact = None, None
		exactfun = normArgs.get('exactFunction', None)
		exactfunders = normArgs.get('exactFunctionDers', None)
		exactextraArgs = normArgs.get('exactExtraArgs', None)
		if exactextraArgs is not None:
			assert isinstance(exactextraArgs, dict), 'Error type of extra args'
			if not exactextraArgs.has_key('position'): exactextraArgs['position'] = qpPhy
			if callable(exactfun): u_exact = np.atleast_2d(exactfun(exactextraArgs))
			if callable(exactfunders): uders_exact = np.atleast_3d(exactfunders(exactextraArgs))
		else:
			if callable(exactfun): u_exact = np.atleast_2d(exactfun(qpPhy))
			if callable(exactfunders): uders_exact = np.atleast_3d(exactfunders(qpPhy))

		part_ref = normArgs.get('part_ref', None); time_ref = normArgs.get('time_ref', None); u_ref = normArgs.get('u_ref', None)
		if isinstance(part_ref, part) and isinstance(time_ref, part1D) and isinstance(u_ref, np.ndarray):
			basisExact, indicesExact = [], [], []
			for i in range(self.part.dim):
				basis, indi, indj = evalDersBasisFortran(part_ref.degree[i], part_ref.knotvector[i], quadPts[i])
				basisExact.append(basis); indicesExact.append(indi); indicesExact.append(indj)
			basis, indi, indj = evalDersBasisFortran(time_ref.degree, time_ref.knotvector, quadPts[-1])
			basisExact.append(basis); indicesExact.append(indi); indicesExact.append(indj)

			sptimectrlpts = np.zeros((part_ref.dim+1, part_ref.nbctrlpts_total*time_ref.nbctrlpts))
			for i in range(time_ref.nbctrlpts):
				iold = i*part_ref.nbctrlpts_total; inew = (i + 1)*part_ref.nbctrlpts_total
				sptimectrlpts[:, iold:inew] = np.vstack([part_ref.ctrlpts, time_ref.ctrlpts[i]*np.ones(part_ref.nbctrlpts_total)])

			inpts = [*nbqp, *indicesExact, *basisExact]
			if self.part.dim == 2:   
				u_exact = geophy.interpolate_meshgrid_3d(*inpts, np.atleast_2d(u_ref))    
				JqpExact = geophy.eval_jacobien_3d(*inpts, sptimectrlpts)
				_, invJExact = geophy.eval_inverse_det(JqpExact) 
				derstemp = geophy.eval_jacobien_3d(*inpts, np.atleast_2d(u_ref))
			elif self.part.dim == 3: 
				u_exact = geophy.interpolate_meshgrid_4d(*inpts, np.atleast_2d(u_ref))
				JqpExact = geophy.eval_jacobien_4d(*inpts, sptimectrlpts)
				_, invJExact = geophy.eval_inverse_det(JqpExact)
				derstemp = geophy.eval_jacobien_4d(*inpts, np.atleast_2d(u_ref))

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
			tmp1 = np.einsum('i,j,k,ijk->', parametricWeights[0], parametricWeights[1], parametricWeights[2], norm1)
			tmp2 = np.einsum('i,j,k,ijk->', parametricWeights[0], parametricWeights[1], parametricWeights[2], norm2)
		if self.part.dim == 3: 
			tmp1 = np.einsum('i,j,k,l,ijkl->', parametricWeights[0], parametricWeights[1], parametricWeights[2], parametricWeights[3], norm1)
			tmp2 = np.einsum('i,j,k,l,ijkl->', parametricWeights[0], parametricWeights[1], parametricWeights[2], parametricWeights[3], norm2)
	
		abserror = np.sqrt(tmp1)
		relerror = np.sqrt(tmp1/tmp2) if tmp2!=0 else np.sqrt(tmp1)
		if tmp2 == 0: print('Warning: Dividing by zero')

		return abserror, relerror

class stheatproblem(stproblem):
	def __init__(self, heat_material:heatmat, part:part, tspan:part1D, boundary:boundaryCondition, solverArgs={}):
		stproblem.__init__(self, part, tspan, boundary, solverArgs)
		self.heatmaterial = heat_material
		if self.heatmaterial.density is None: self.heatmaterial.density = lambda x: np.ones(self.part.nbqp_total)
		return
	
	def compute_mfSTConductivity(self, array_in, args=None):
		if args is None: args = {'position': self.part.qpPhy}
		prop = self.heatmaterial.conductivity(args)
		inpts = [*self._getInputs(), self.part.invJ, self.part.detJ, self.time.detJ, prop]
		if self.part.dim == 2: array_out = stheatsolver.mf_stconductivity_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = stheatsolver.mf_stconductivity_3d(*inpts, array_in)
		return array_out
	
	def compute_mfSTCapacity(self, array_in, args=None):
		if args is None: args = {'position': self.part.qpPhy}
		prop = self.heatmaterial.capacity(args)
		inpts = [*self._getInputs(), self.part.invJ, self.part.detJ, self.time.detJ, prop]
		if self.part.dim == 2: array_out = stheatsolver.mf_stcapacity_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = stheatsolver.mf_stcapacity_3d(*inpts, array_in)
		return array_out

	def interpolate_STtemperature_gradients(self, u_ctrlpts):
		inpts = [*self.part.nbqp[:self.part.dim], self.time.nbqp, *self.part.indices, 
				*self.time.quadRule.dersIndices, *self.part.basis, self.time.quadRule.dersBasis]
		
		sptimectrlpts = np.zeros((self.part.dim+1, self.part.nbctrlpts_total*self.time.nbctrlpts_total))
		for i in range(self.time.nbctrlpts_total):
			iold = i*self.part.nbctrlpts_total; inew = (i + 1)*self.part.nbctrlpts_total
			sptimectrlpts[:, iold:inew] = np.vstack([self.part.ctrlpts, self.time.ctrlpts[i]*np.ones(self.part.nbctrlpts_total)])
		
		if self.part.dim == 2:
			Jqp = geophy.eval_jacobien_3d(*inpts, sptimectrlpts)
			_, invJ = geophy.eval_inverse_det(Jqp)
			u_interp = geophy.interpolate_meshgrid_3d(*inpts, np.atleast_2d(u_ctrlpts))
			derstemp = geophy.eval_jacobien_3d(*inpts, np.atleast_2d(u_ctrlpts))
	
		elif self.part.dim == 3:
			Jqp = geophy.eval_jacobien_4d(*inpts, sptimectrlpts)
			_, invJ = geophy.eval_inverse_det(Jqp)
			u_interp = geophy.interpolate_meshgrid_4d(*inpts, np.atleast_2d(u_ctrlpts))
			derstemp = geophy.eval_jacobien_4d(*inpts, np.atleast_2d(u_ctrlpts))

		u_interp = np.atleast_2d(u_interp); u_interp = np.ravel(u_interp)
		uders_interp = np.atleast_3d(np.einsum('ijl,jkl->ikl', derstemp, invJ)); uders_interp = uders_interp[0, :, :]
		return u_interp, uders_interp

	def _compute_STHeatIntForce(self, array_in, args=None):
		assert args is not None, 'Please enter a valid argument'
		intForce = self.compute_mfSTCapacity(array_in, args) + self.compute_mfSTConductivity(array_in, args)
		return intForce
	
	def _solveLinearizedSTHeatProblem(self, Fext, args=None, isfull=False, threshold=None):
		assert args is not None, 'Please enter a valid argument'
		if threshold is None: threshold = self._thresLin
		Cprop = self.heatmaterial.capacity(args)
		Kprop = self.heatmaterial.conductivity(args)
		if isfull:
			gradTemperature = args['gradients']
			Cdersprop = self.heatmaterial.capacityDers(args)*gradTemperature[-1, :]
			Kdersprop = np.einsum('ijk,jk->ik', self.heatmaterial.conductivityDers(args), gradTemperature[:self.part.dim, :])
			inpts = [*self._getInputs(), self.boundary.thDirichletTable, self.part.invJ, self.part.detJ,
					self.time.detJ, Cprop, Cdersprop, Kprop, Kdersprop, Fext, self._itersLin, 
					threshold, self._linPreCond, self._linSolv]
			if self.part.dim == 2: temperature, residue = stheatsolver.solver_newtonspacetime_heat_2d(*inpts)
			if self.part.dim == 3: temperature, residue = stheatsolver.solver_newtonspacetime_heat_3d(*inpts)
		else:
			inpts = [*self._getInputs(), self.boundary.thDirichletTable, self.part.invJ, self.part.detJ,
					self.time.detJ, Cprop, Kprop, Fext, self._itersLin, threshold, self._linPreCond, self._linSolv]
			if self.part.dim == 2: temperature, residue = stheatsolver.solver_picardspacetime_heat_2d(*inpts)
			if self.part.dim == 3: temperature, residue = stheatsolver.solver_picardspacetime_heat_3d(*inpts)
		return temperature, residue
	
	def solveFourierSTHeatProblem(self, Tinout, Fext, isfull=False, isadaptive=True, solvArgs={}):
		eps_kr0  = solvArgs.get('initial', .5)
		gamma_kr = solvArgs.get('coefficient', 1.0)
		omega_kr = solvArgs.get('exponential', 2.0)

		dod = self.boundary.getThermalBoundaryConditionInfo()[0]
		AllresLin, AllresNewton, Allsol, Allthres, Alldelta = [], [], [], [], []
		threshold_inner = None
		for j in range(self._itersNL):
			# enablePrint()

			# Compute temperature at each quadrature point
			temperature, gradtemperature = self.interpolate_STtemperature_gradients(Tinout)

			# Compute internal force
			Fint_dj = self._compute_STHeatIntForce(Tinout, args={'temperature':temperature})

			# Compute residue
			r_dj = Fext - Fint_dj
			r_dj[dod] = 0.0

			resNLj1 = np.sqrt(np.dot(r_dj, r_dj))
			print('Nonlinear error: %.3e' %resNLj1)

			if j == 0: resNL0 = np.copy(resNLj1)
			if resNLj1 <= max([self._safeguard, self._thresNL*resNL0]): break
			AllresNewton.append(resNLj1); Allsol.append(np.copy(Tinout))

			# Update inner threshold
			if isadaptive: 
				if j == 0: 
					threshold_ref = np.copy(eps_kr0)

				else:
					ratio = resNLj1/resNLj0
					eps_kr_k = gamma_kr*np.power(ratio, omega_kr)
					eps_kr_r = gamma_kr*np.power(threshold_inner, omega_kr)
					print('Ratio:%.3e' %ratio)
					print('%.3e, %.3e' %(eps_kr_r, eps_kr_k))
					if eps_kr_r <= 0.1: 
						threshold_ref = np.copy(eps_kr_k)
						print('Method 1')
					else: 
						threshold_ref = max([eps_kr_k, eps_kr_r])
						print('Method 2')

				print('Thres_ref: %.3e' %threshold_ref)
				threshold_inner = min([eps_kr0, max([self._thresLin, threshold_ref])])

			else: 
				threshold_inner = np.copy(self._thresLin)
			Allthres.append(threshold_inner)

			# Solve for active control points
			deltaD, resLinj = self._solveLinearizedSTHeatProblem(r_dj, 
										args={'temperature':temperature, 
											'gradients':gradtemperature}, 
										isfull=isfull, threshold=threshold_inner)			

			# Update active control points
			Tinout += deltaD
			resNLj0 = np.copy(resNLj1)
			AllresLin.append(resLinj); Alldelta.append(deltaD)

			blockPrint()

			
		output = {'KrylovRes': AllresLin, 'NewtonRes':AllresNewton, 'Solution':Allsol, 'Threshold':Allthres, 'Delta':Alldelta}
		return output
