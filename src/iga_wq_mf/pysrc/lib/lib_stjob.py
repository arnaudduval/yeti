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
		self._nbIterPCG    = solverArgs.get('nbIterationsPCG', 100)
		self._nbIterNR     = solverArgs.get('nbIterationsNR', 50)
		self._thresholdPCG = solverArgs.get('PCGThreshold', 1e-12)
		self._thresholdNR  = solverArgs.get('NRThreshold', 1e-10)
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
	
	def solveFourierSTHeatProblem(self, Tinout, Fext):
		dod = self.boundary.getThermalBoundaryConditionInfo()[0]
		dj_n1 = np.copy(Tinout); d_n1ref = np.copy(Tinout)
		
		AllresPCG = []
		for j in range(self._nbIterNR):

			# Compute temperature at each quadrature point
			temperature = self.interpolate_STtemperature(Tinout)
		
			# Compute internal force
			Fint_dj = self.compute_STHeatIntForce(dj_n1, args=temperature)

			# Compute residue
			r_dj = Fext - Fint_dj
			r_dj[dod] = 0.0

			# Solve for active control points
			resPCGj = np.array([j+1])
			deltaD, resPCG = self._solveLinearizedSTHeatProblem(r_dj, args=temperature)
			resPCGj = np.append(resPCGj, resPCG)
			d_n1ref += deltaD

			# Compute residue of Newton Raphson using an energetic approach
			resNRj = abs(np.dot(d_n1ref, r_dj))
			if j == 0: resNR0 = resNRj
			print('NR error: %.5e' %resNRj)
			if resNRj <= self._thresholdNR*resNR0: break

			# Update active control points
			dj_n1 += deltaD
			AllresPCG.append(resPCGj)

		Tinout = np.copy(dj_n1)

		return AllresPCG
