from lib.lib_job import *
from lib.lib_quadrules import GaussQuadrature, WeightedQuadrature

class heatproblemSpTm(heatproblem):

	def __init__(self, material:thermomat, part:part, boundary:boundaryCondition, timeArgs:dict):
		super().__init__(material, part, boundary)
		self.timeDisc = self.timeDiscretisation()
		self.__getInfo(timeArgs)
		self.__update_thDirichletBound()
		return
	
	def __getInfo(self, timeArgs:dict):
		tmspan       = timeArgs.get('timespan', 1.)
		tmdegree     = timeArgs.get('degree', None)
		tmknotvector = timeArgs.get('knotvector', None)
		if any(el for el in [tmdegree, tmknotvector]): raise Warning('Not possible')
		tmquadrule   = timeArgs.get('quadrule', 'wq').lower()
		if   tmquadrule == 'iga': quadRule = GaussQuadrature(tmdegree, tmknotvector, quadArgs=timeArgs)
		elif tmquadrule == 'wq' : quadRule = WeightedQuadrature(tmdegree, tmknotvector, quadArgs=timeArgs)
		else: raise Warning('Not found')
		quadRule.getQuadratureRulesInfo()
		self.timeDisc.nbctrlpts = quadRule.nbctrlpts
		self.timeDisc.indices = quadRule.dersIndices
		self.timeDisc.basis   = quadRule.dersBasis
		self.timeDisc.weights = quadRule.dersWeights
		self.timeDisc.nbqp    = len(quadRule.quadPtsPos)
		self.timeDisc.detG    = tmspan*np.ones(self.timeDisc.nbqp)
		self.timeDisc.qpPhy   = quadRule.quadPtsPos*self.timeDisc.detG
		return
	
	def __update_thDirichletBound(self):
		nbctrlpts_total = np.product(self.boundary._nbctrlpts)*self.timeDisc.nbctrlpts
		spthDirichletBound = self.boundary.thDirichletBound
		spdod = self.boundary.thdod
		dod   = np.arange(self.boundary._nbctrlpts, dtype=int)
		thDirichletBound = np.zeros(nbctrlpts_total)
		for i in range(self.timeDisc.nbctrlpts):
			tmp = spdod + i*self.boundary._nbctrlpts
			dod = np.append(dod, tmp)
			thDirichletBound[tmp] = spthDirichletBound
		dod = set(dod)
		dof = set(np.arange(nbctrlpts_total, dtype=int)).difference(spdod)
		self.boundary.thdod = np.sort(np.array(list(spdod), dtype=int))
		self.boundary.thdof = np.sort(np.array(list(dof), dtype=int))
		self.boundary.thDirichletBound = thDirichletBound
		return

	def add_initialTemperature(self, temperature=0.):
		dod = np.arange(self.boundary._nbctrlpts, dtype=int)
		self.boundary.thDirichletBound[dod] = temperature
		return

	def get_input4MatrixFree(self, sptable=None):
		" Returns necessary inputs to compute the product between a matrix and a vector "
		
		if sptable is None: sptable = self.boundary.thDirichletTable
		indices, basis, weights = [], [], []
		for i in range(self.part.dim):
			# Select data
			if np.array_equal(sptable[i, :], [0, 0]): rows2erase = []
			if np.array_equal(sptable[i, :], [0, 1]): rows2erase = [-1]
			if np.array_equal(sptable[i, :], [1, 0]): rows2erase = [0]
			if np.array_equal(sptable[i, :], [1, 1]): rows2erase = [0, -1]
			indi_t, indj_t, data_t = eraseRowsCSR(rows2erase, 
									self.part.indices[2*i], self.part.indices[2*i+1],  
									[self.part.basis[i], self.part.weights[i]])
			
			# Extract data and append to list
			[basist, weightst] = data_t
			indices.append(indi_t); indices.append(indj_t) 
			basis.append(basist); weights.append(weightst)

		rows2erase = [0]
		indi_t, indj_t, data_t = eraseRowsCSR(rows2erase, 
								self.timeDisc.indices[0], self.timeDisc.indices[1],  
								[self.timeDisc.basis, self.timeDisc.weights])
		
		# Extract data and append to list
		[basist, weightst] = data_t
		indices.append(indi_t); indices.append(indj_t) 
		basis.append(basist); weights.append(weightst)

		inputs = [*self.part.nbqp, self.timeDisc.nbqp, *indices, *basis, *weights]

		return inputs
	
	def eval_mfSpTmMatrix(self, u, table=None):
		inputs  = self.get_input4MatrixFree(table=table)
		Kcoefs  = self.material.eval_conductivityCoefficients(self.part.invJ, self.part.detJ, self.part.qpPhy)
		Ccoefs  = self.material.eval_capacityCoefficients(self.part.detJ, self.part.qpPhy)

		if self.part.dim == 2: raise Warning('Until now not done')
		if self.part.dim == 3: result = heatsolver.mf_wq_get_spacetimeheat_3d(Ccoefs, Kcoefs, self.timeDisc.detG, *inputs, u)
		return result
	
	class timeDiscretisation():
		def __init__(self):
			self.nbctrlpts = None
			self.indices = None
			self.basis   = None
			self.weights = None
			self.detG    = None 
			self.qpPhy   = None
			self.nbqp    = None
			return