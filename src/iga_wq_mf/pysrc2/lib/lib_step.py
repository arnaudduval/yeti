from lib.__init__ import *

class step():

	def __init__(self, nbctrlpts=np.array([1, 1, 1])):
		self._dim                 = np.count_nonzero(nbctrlpts>1)
		if self._dim == 1: raise Warning('Not possible')
		self._nbctrlpts           = nbctrlpts

		self._thDirichletTable    = np.zeros((self._dim, 2), dtype=bool)
		self._thDirichletBound    = None
		self._thdof  = []
		self._thdod  = []

		self._mchDirichletTable   = np.zeros((self._dim, 2, self._dim), dtype=bool)
		self._mchDirichletBound   = None
		self._mchdof = []
		self._mchdod = []
		
		return
	
	def activate_DirichletThermal(self):
		nbctrlpts_total = np.product(self._nbctrlpts)
		if self._thDirichletBound is None:
			self._thDirichletBound = np.zeros(nbctrlpts_total)
		return
	
	def activate_DirichletMechanical(self):
		nbctrlpts_total = np.product(self._nbctrlpts)
		dimen = np.size(self._nbctrlpts)
		if self._mchDirichletBound is None:
			self._mchDirichletBound = np.zeros((nbctrlpts_total, dimen))
		return
	
	def clear_Dirichlet(self):
		" Clears Dirichlet boundaries "
		self._thDirichletTable = None
		self._thdof = np.array([], dtype=int)
		self._thdod = np.array([], dtype=int)
		self._mchDirichletTable = None
		self._mchdof = np.array([], dtype=int)
		self._mchdod = np.array([], dtype=int)
		return
		
	def __get_boundaryNodes(self, table, nbctrlpts, dimen=3): 
		" Gets the indices of the blocked and free control points from table"

		def get_INCTable(nnzByDimension):
			" Sets topology table, also known as INC: NURBS coordinates. "
			# Create INC: NURBS coordinates
			nnz_total = np.prod(nnzByDimension)
			table = np.zeros((nnz_total, 3), dtype= int)
			for i3 in range(nnzByDimension[2]): 
				for i2 in range(nnzByDimension[1]): 
					for i1 in range(nnzByDimension[0]):
						genPos = i1 + i2*nnzByDimension[0] + i3*nnzByDimension[0]*nnzByDimension[1]
						table[genPos, :] = [i1, i2, i3]
			return table

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
	
	# Heat problem

	def update_thDirichletBound(self):
		nbctrlpts_total = np.product(self._nbctrlpts)
		dod = set(self._thdod)
		dof = set(np.arange(nbctrlpts_total, dtype=int)).difference(dod)
		self._thdod = np.array(list(dod), dtype=int)
		self._thdof = np.array(list(dof), dtype=int)
		return

	def add_DirichletTemperature(self, table=None, temperature=0.0):
		"This function is first tentative of adding constant boundary conditions "
		table = np.array(table, dtype=bool)
		if not np.any(table == True): raise Warning('At least one blocked face is needed')
		self.activate_DirichletThermal()
		dod = self.__get_boundaryNodes(table, self._nbctrlpts, dimen=self._dim)[0]
		self._thDirichletTable += table
		
		if np.isscalar(temperature): 
			self._thDirichletBound[dod] += temperature*np.ones(len(dod))
		else: 
			if len(temperature) != len(dod): raise Warning('Not possible')
			self._thDirichletBound[dod] += temperature		
		self._thdod = np.append(self._thdod, dod)
		self._thdod = np.array(self._thdod, dtype=int)
		self.update_thDirichletBound()
		return 

	def getThermalBoundaryConditionInfo(self): 
		if self._thDirichletTable is None: raise Warning('Please define first total Dirichlet boundaries')
		self.update_thDirichletBound()
		return  self._thdod, self._thDirichletBound[self._thdod], self._thdof
	
	# Mechanical problem
