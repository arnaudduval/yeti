"SPACE TIME GALERKIN"
from .__init__ import *
from .lib_base import get_faceInfo, get_INCTable, evalDersBasisFortran
from .lib_quadrules import GaussQuadrature
from .lib_material import (heatmat, clean_dirichlet, block_dot_product)
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
		self._nbIterPCG    = solverArgs.get('nbIterationsPCG', 50)
		self._nbIterNR     = solverArgs.get('nbIterationsNR', 50)
		self._thresholdPCG = solverArgs.get('PCGThreshold', 1e-12)
		self._thresholdNR  = solverArgs.get('NRThreshold', 1e-4)
		self._methodPCG    = solverArgs.get('PCGmethod', 'JMC')
		return
	
	# def __getInfo4surfForce(self, nbFacePosition):
	# 	# INC_ctrlpts = get_INCTable(self.part.nbctrlpts)
	# 	# direction, side = get_faceInfo(nbFacePosition)
	# 	# if direction>=2*self.part.dim: raise Warning('Not possible')
	# 	# valrange = [i for i in range(self.part.dim)]
	# 	# valrange.pop(direction)
	# 	# if side == 0:   CPList = np.where(INC_ctrlpts[:, direction] == 0)[0]
	# 	# elif side == 1: CPList = np.where(INC_ctrlpts[:, direction] == self.part.nbctrlpts[direction]-1)[0]
	# 	# CPList = list(np.sort(CPList))

	# 	# nnz, indices, basis, weights = [], [], [], []
	# 	# for _ in valrange:
	# 	# 	nnz.append(self.part.nbqp[_]); basis.append(self.part.basis[_]); weights.append(self.part.weights[_])
	# 	# 	indices.append(self.part.indices[2*_]); indices.append(self.part.indices[2*_+1]) 
	# 	# inpts = [*nnz, *indices, *basis, self.part.ctrlpts[:, CPList]]
	# 	# if self.part.dim == 2: 
	# 	# 	Jqp = geophy.eval_jacobien_1d(*inpts)
	# 	# 	qpPhy = geophy.interpolate_meshgrid_1d(*inpts)
	# 	# elif self.part.dim == 3:
	# 	# 	Jqp = geophy.eval_jacobien_2d(*inpts)
	# 	# 	qpPhy = geophy.interpolate_meshgrid_2d(*inpts)
	# 	return # nnz, indices, basis, weights, Jqp, qpPhy, CPList

	# def compute_surfForce(self, surffun, nbFacePosition):
	# 	""" Computes the surface foce over the boundary of a geometry. 
	# 		The surffun is a Neumann like function, ie, in transfer heat q = -(k grad(T)).normal
	# 		and in elasticity t = sigma.normal
	# 	"""
	# 	# nnz, indices, _, weights, Jqp, qpPhy, CPList = self.__getInfo4surfForce(nbFacePosition)
	# 	# prop = surffun(qpPhy); prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
	# 	# inpts = [*nnz, *indices, *weights, Jqp, prop]
	# 	# if   self.part.dim == 2: tmp = geophy.get_forcesurf_2d(*inpts)
	# 	# elif self.part.dim == 3: tmp = geophy.get_forcesurf_3d(*inpts)
	# 	# surfForce = np.zeros((nr, self.part.nbctrlpts_total))
	# 	# surfForce[:, CPList] = tmp
	# 	# if nr == 1: surfForce = np.ravel(surfForce)
	# 	return # surfForce, CPList
	
	# def compute_surfForce_woNormal(self, gradfun, nbFacePosition):
	# 	# nnz, indices, _, weights, Jqp, qpPhy, CPList = self.__getInfo4surfForce(nbFacePosition)

	# 	# normal_qpPhy = geophy.eval_normal(Jqp)
	# 	# if np.any(np.array([0, 3, 4], dtype=int) == nbFacePosition): normal_qpPhy = -normal_qpPhy 
	# 	# # By the moment I do not why in this surfaces the normal has the wrong sign (the same for 2 and 3 dimensions)
	# 	# grad_qpPhy = gradfun(qpPhy)
	# 	# if np.array(grad_qpPhy).ndim == 3:
	# 	# 	prop = np.einsum('ijl,jl->il', grad_qpPhy, normal_qpPhy)
	# 	# elif np.array(grad_qpPhy).ndim == 2:
	# 	# 	prop = np.einsum('jl,jl->l', grad_qpPhy, normal_qpPhy)
	# 	# else: raise Warning('Size problem')

	# 	# prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
	# 	# inpts = [*nnz, *indices, *weights, Jqp, prop]
	# 	# if   self.part.dim == 2: tmp = geophy.get_forcesurf_2d(*inpts)
	# 	# elif self.part.dim == 3: tmp = geophy.get_forcesurf_3d(*inpts)
	# 	# surfForce = np.zeros((nr, self.part.nbctrlpts_total))
	# 	# surfForce[:, CPList] = tmp
	# 	# if nr == 1: surfForce = np.ravel(surfForce)
	# 	return # surfForce, CPList
	
	def compute_volForce(self, volfun, args=None): 
		" Computes the volume force over a geometry "
		assert args is not None, 'Please enter a valid argument'
		prop = volfun(args)
		prop = np.atleast_2d(prop); nr = np.size(prop, axis=0)
		inpts = [*self.part.nbqp[:self.part.dim], self.time.nbqp, *self.part.indices, *self.time.quadRule.dersIndices, 
				*self.part.weights, self.time.quadRule.dersWeights, self.part.detJ, self.time.detJ, prop]
		if self.part.dim == 2: volForce = geophy.get_forcevol_st2d(*inpts)
		if self.part.dim == 3: volForce = geophy.get_forcevol_st3d(*inpts)
		if nr == 1: volForce = np.ravel(volForce)
		return volForce	
	
	def normOfError(self, u_ctrlpts, normArgs:dict, isRelative=True):
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

		sptimectrlpts = np.zeros((4,self.part.nbctrlpts_total*self.time.nbctrlpts))
		iold = 0
		for i in range(self.time.nbctrlpts):
			inew = (i + 1)*self.part.nbctrlpts_total
			sptimectrlpts[:, iold:inew] = np.stack([self.part.ctrlpts, self.time.ctrlpts[i]*np.ones(self.part.nbctrlpts_total)])
			iold = np.copy(inew)

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

		u_interp = np.atleast_2d(u_interp); uders_interp = None
		uders_interp = np.atleast_3d(np.einsum('ijl,jkl->ikl', derstemp, invJ))

		# Compute u exact
		u_exact, uders_exact = None, None
		exactfun = normArgs.get('exactFunction', None)
		exactfunders = normArgs.get('exactFunctionDers', None)
		if callable(exactfun): u_exact = np.atleast_2d(exactfun(qpPhy))
		if callable(exactfunders): uders_exact = np.atleast_3d(exactfunders(qpPhy))

		part_ref = normArgs.get('part_ref', None); time_ref = normArgs.get('time_ref', None); u_ref = normArgs.get('u_ref', None)
		if isinstance(part_ref, part) and isinstance(time_ref, part1D) and isinstance(u_ref, np.ndarray):
			nbqpExact, basisExact, indicesExact = [], [], []
			for i in range(self.part.dim):
				basis, indi, indj = evalDersBasisFortran(part_ref.degree[i], part_ref.knotvector[i], quadPts[i])
				nbqpExact.append(len(quadPts[i])); basisExact.append(basis); indicesExact.append(indi); indicesExact.append(indj)
			basis, indi, indj = evalDersBasisFortran(time_ref.degree, time_ref.knotvector, time_ref.quadRule.quadPtsPos)
			nbqpExact.append(len(quadPts[-1])); basisExact.append(basis); indicesExact.append(indi); indicesExact.append(indj)
			
			inpts = [*nbqpExact, *indicesExact, *basisExact]
			if self.part.dim == 2:   
				u_exact = geophy.interpolate_meshgrid_3d(*inpts, np.atleast_2d(u_ref))    
				JqpExact = geophy.eval_jacobien_3d(*inpts, part_ref.ctrlpts)
				_, invJExact = geophy.eval_inverse_det(JqpExact) 
				derstemp = geophy.eval_jacobien_3d(*inpts, np.atleast_2d(u_ref))
			elif self.part.dim == 3: 
				u_exact = geophy.interpolate_meshgrid_4d(*inpts, np.atleast_2d(u_ref))
				JqpExact = geophy.eval_jacobien_4d(*inpts, part_ref.ctrlpts)
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
			
		if isRelative: error = np.sqrt(tmp1/tmp2)
		else:          error = np.sqrt(tmp1)

		return


class stheatproblem(stproblem):
	def __init__(self, heat_material:heatmat, part:part, tspan:part1D, boundary:boundaryCondition, solverArgs={}):
		stproblem.__init__(self, part, tspan, boundary, solverArgs)
		self.heatmaterial = heat_material
		if self.heatmaterial.density is None: self.heatmaterial.addDensity(inpt=1.0, isIsotropic=True)
		return
	
	def compute_mfSTConductivity(self, array_in, args=None):
		assert args is not None, 'Please enter a valid argument'
		prop = self.heatmaterial.conductivity(args)*self.heatmaterial.density(args)
		inpts = [*self._getInputs(), self.part.invJ, self.part.detJ, self.time.detJ, prop]
		if self.part.dim == 2: array_out = stheatsolver.mf_stconductivity_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = stheatsolver.mf_stconductivity_3d(*inpts, array_in)
		return array_out
	
	def compute_mfSTCapacity(self, array_in, args=None):
		assert args is not None, 'Please enter a valid argument'
		prop = self.heatmaterial.capacity(args)*self.heatmaterial.density(args)
		inpts = [*self._getInputs(), self.part.invJ, self.part.detJ, self.time.detJ, prop]
		if self.part.dim == 2: array_out = stheatsolver.mf_stcapacity_2d(*inpts, array_in)
		if self.part.dim == 3: array_out = stheatsolver.mf_stcapacity_3d(*inpts, array_in)
		return array_out

	def interpolate_STtemperature(self, T_ctrlpts):
		inpts = [*self.part.nbqp[:self.part.dim], self.time.nbqp, *self.part.indices, 
			*self.time.quadRule.dersIndices, *self.part.basis, self.time.quadRule.dersBasis, 
			np.atleast_2d(T_ctrlpts)]
		if self.part.dim == 2:   T_interp = geophy.interpolate_meshgrid_3d(*inpts)
		elif self.part.dim == 3: T_interp = geophy.interpolate_meshgrid_4d(*inpts)
		T_interp = np.ravel(T_interp)
		return T_interp

	def compute_STHeatIntForce(self, array_in, args=None):
		assert args is not None, 'Please enter a valid argument'
		intForce = self.compute_mfSTCapacity(array_in, args) + self.compute_mfSTConductivity(array_in, args)
		return intForce
	
	def _solveLinearizedSTHeatProblem(self, Fext, args=None):
		assert args is not None, 'Please enter a valid argument'
		Cprop = self.heatmaterial.capacity(args)*self.heatmaterial.density(args)
		Kprop = self.heatmaterial.conductivity(args)
		inpts = [*self._getInputs(), self.boundary.thDirichletTable, self.part.invJ, self.part.detJ,
				self.time.detJ, Cprop, Kprop, Fext, self._nbIterPCG, self._thresholdPCG, self._methodPCG]
		if self.part.dim == 2: temperature, residue = stheatsolver.solver_linearspacetime_heat_2d(*inpts)
		if self.part.dim == 3: temperature, residue = stheatsolver.solver_linearspacetime_heat_3d(*inpts)
		return temperature, residue
	
	def solveFourierSTHeatProblem(self, Tguess, Fext):
		dod = self.boundary.getThermalBoundaryConditionInfo()[0]
		dj_n1 = np.copy(Tguess); d_n1ref = np.copy(Tguess)
		
		AllresPCG = []
		for j in range(self._nbIterNR):

			# Compute temperature at each quadrature point
			temperature = self.interpolate_STtemperature(Tguess)
		
			# Compute internal force
			Fint_dj = self.compute_STHeatIntForce(dj_n1, args=temperature)

			# Compute residue
			r_dj = Fext - Fint_dj
			r_dj[dod] = 0.0

			# Solve for active control points
			deltaD, resPCGj = self._solveLinearizedSTHeatProblem(r_dj, args=temperature)
			d_n1ref += deltaD

			# Compute residue of Newton Raphson using an energetic approach
			resNRj = abs(np.dot(d_n1ref, r_dj))
			if j == 0: resNR0 = resNRj
			print('NR error: %.5e' %resNRj)
			if resNRj <= self._thresholdNR*resNR0: break

			# Update active control points
			dj_n1 += deltaD
			AllresPCG.append(resPCGj)

		return dj_n1, AllresPCG
