from .__init__ import *
from .lib_base import get_INCTable

class boundaryCondition():

	def __init__(self, nbctrlpts=np.array([1, 1, 1])):
		self._nbctrlpts = nbctrlpts
		self._dim = np.count_nonzero(nbctrlpts>1)
		if self._dim == 1: raise Warning('Not possible')	
		self.clear_Dirichlet()
		return
	
	def __activate_DirichletThermal(self):
		nbctrlpts_total = np.product(self._nbctrlpts)
		if self.thDirichletBound is None:
			self.thDirichletBound = np.zeros(nbctrlpts_total)
		return
	
	def __activate_DirichletMechanical(self):
		nbctrlpts_total = np.product(self._nbctrlpts)
		if self.mchDirichletBound is None:
			self.mchDirichletBound = np.zeros((nbctrlpts_total, self._dim))
		return
	
	def __get_boundaryNodes(self, table, nbctrlpts, dimen=3): 
		" Gets the indices of the blocked and free control points from table"

		# The table of dirichlet boundaries must be at least 3D
		tablecopy = np.atleast_3d(np.copy(table))

		# Get number of degree of freedom (DOF) per node
		nbDOF = np.shape(tablecopy)[2]
		if np.shape(tablecopy)[0] < dimen or np.shape(tablecopy)[1] != 2:
			raise Warning('Table is not well defined')

		# Find nodes
		INC = get_INCTable(nbctrlpts)
		dod_total = []
		for i in range(nbDOF):
			dod = []
			for j in range(dimen):
				block_bound_dim = tablecopy[j, :, i]
				
				if block_bound_dim[0]: 
					dod.extend(np.where(INC[:, j] == 0)[0])

				if block_bound_dim[1]:
					dod.extend(np.where(INC[:, j] == nbctrlpts[j]-1)[0])

			# Rearrange
			dod = np.unique(dod); dod = np.array(dod, dtype=int)
			dod_total.append(dod)

		return dod_total
	
	def clear_Dirichlet(self):
		" Clears Dirichlet boundaries "
		self.thDirichletTable = np.zeros((self._dim, 2), dtype=bool)
		self.thDirichletBound  = None
		self.thdof  = []
		self.thdod  = []
		self.mchDirichletTable = np.zeros((self._dim, 2, self._dim), dtype=bool)
		self.mchDirichletBound  = None
		self.mchdof = [[] for i in range(self._dim)]
		self.mchdod = [[] for i in range(self._dim)]
		return
	

	# Heat problem

	def __update_thDirichletBound(self):
		nbctrlpts_total = np.product(self._nbctrlpts)
		dod = set(self.thdod)
		dof = set(np.arange(nbctrlpts_total, dtype=int)).difference(dod)
		self.thdod = np.sort(np.array(list(dod), dtype=int))
		self.thdof = np.sort(np.array(list(dof), dtype=int))
		return

	def add_DirichletConstTemperature(self, table=None, temperature=0.0):
		"This function is first tentative of adding constant boundary conditions "
		table = np.array(table, dtype=bool)
		if not np.any(table == True): raise Warning('At least one blocked face is needed')
		self.__activate_DirichletThermal()
		dod_total = self.__get_boundaryNodes(table, self._nbctrlpts, dimen=self._dim)[0]
		self.thDirichletTable += table
		
		if np.isscalar(temperature): 
			self.thDirichletBound[dod_total] = temperature
		else: raise Warning('Not possible')
		tmp = np.append(self.thdod, dod_total)
		self.thdod = np.array(tmp, dtype=int)
		self.__update_thDirichletBound()
		return 

	def getThermalBoundaryConditionInfo(self): 
		if self.thDirichletTable is None: raise Warning('Please define first total Dirichlet boundaries')
		self.__update_thDirichletBound()
		return  self.thdod, self.thDirichletBound[self.thdod], self.thdof
	
	# Mechanical problem

	def __update_mchDirichletBound(self):
		nbctrlpts_total = np.product(self._nbctrlpts)
		self.mchdof    = [[] for i in range(self._dim)]
		for i, dod in enumerate(self.mchdod):
			dod = set(dod)
			dof = set(np.arange(nbctrlpts_total, dtype=int)).difference(dod)
			self.mchdod[i] = np.sort(np.array(list(dod), dtype=int))
			self.mchdof[i] = np.sort(np.array(list(dof), dtype=int))
		return

	def add_DirichletDisplacement(self, table=None, displacement=0.0):
		"This function is first tentative of adding constant boundary conditions "
		table = np.array(table, dtype=bool)
		if not np.any(table == True): raise Warning('At least one blocked face is needed')
		self.__activate_DirichletMechanical()
		dod_total = self.__get_boundaryNodes(table, self._nbctrlpts, dimen=self._dim)
		self.mchDirichletTable += table
		
		if np.isscalar(displacement): 
			for i, dod in enumerate(dod_total):
				self.mchDirichletBound[dod, i] = displacement*np.ones(len(dod))
		else: 
			if np.size(displacement, axis=0) != len(dod_total) and np.size(displacement, axis=1) != len(dod_total[0]): 
				raise Warning('Not possible')
			for i, dod in enumerate(dod_total):
				self.mchDirichletBound[dod, i] = displacement	

		for i, dod in enumerate(dod_total):
			tmp = np.append(self.mchdod[i], dod)
			self.mchdod[i] = np.array(tmp, dtype=int)
		self.__update_mchDirichletBound()
		return 