from lib.__init__ import *

class step():

	def __init__(self, nbctrlpts):
		nbctrlpts_total           = np.product(nbctrlpts)
		self._thermalDirichlet    = None
		self._thermal_dof         = None
		self._thermal_dod         = None
		self._mechanicalDirichlet = None
		self._mechanical_dof      = None
		self._mechanical_dod      = None
		self._thermalNeumann      = None
		self._mechanicalNeumann   = None
		self._ThermalDirichletBound      = np.zeros(nbctrlpts_total)
		return
	
	def get_INCTable(self, nnzByDimension):
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
	
	def get_boundaryNodes(self, table, nbctrlpts, dimen=3): 
		" Gets the indices of the blocked and free control points from table"

		# The table of dirichlet boundaries must be at least 3D
		table = np.atleast_3d(table)

		# Get number of degree of freedom (DOF) per node
		nbDOF = np.shape(table)[2]
		if np.shape(table)[0] < dimen or np.shape(table)[1] != 2:
			raise Warning('Table is not well defined')

		# Find nodes
		nbctrlpts_total = np.product(nbctrlpts)
		INC = self.get_INCTable(nbctrlpts)
		dod_total = []; dof_total = []
		for i in range(nbDOF):
			dod = []
			for j in range(dimen):
				block_bound_dim = table[j, :, i]
				
				if block_bound_dim[0]: 
					dod.extend(np.where(INC[:, j] == 0)[0])

				if block_bound_dim[1]:
					dod.extend(np.where(INC[:, j] == nbctrlpts[j]-1)[0])

			# Rearrange
			dod = np.unique(dod)
			dof = set(np.arange(nbctrlpts_total, dtype= int)) - set(dod)
			dod_total.append(list(dod))
			dof_total.append(list(dof))

		return dof_total, dod_total

	def set_DirichletConditions(self, Dirichlet:dict):
		" Gets free and blocked control points "

		# Thermal 
		try: 
			TTable     = Dirichlet['thermal']
			Tdof, Tdod = self.get_boundaryNodes(table=TTable)
			Tdof, Tdod = Tdof[0], Tdod[0]
		except: 
			TTable, Tdof, Tdod = None, None, None
		self._thermalDirichlet = TTable
		self._thermal_dof = Tdof
		self._thermal_dod = Tdod

		# Mechanical 
		try: 
			MTable     = Dirichlet['mechanical']
			Mdof, Mdod = self.get_boundaryNodes(table=MTable)
		except: 
			MTable, Mdof, Mdod = None, None, None
		self._mechanicalDirichlet = MTable
		self._mechanical_dof = Mdof
		self._mechanical_dod = Mdod
		
		return
	
	def set_NeumannConditions(self, Neumann:dict):
		" Gets Neumann control points and quadrature points"
		
		# Thermal 
		try:    FTable = Neumann['thermal']
		except: FTable = None
		self._thermalNeumann = FTable

		# Mechanical 
		try:    FTable = Neumann['mechanical']
		except: FTable = None
		self._mechanicalNeumann = FTable

		return 

	def clear_Dirichlet(self):
		" Clears Dirichlet boundaries "
		self._thermalDirichlet = None
		self._thermal_dof = None
		self._thermal_dod = None
		self._mechanicalDirichlet = None
		self._mechanical_dof = None
		self._mechanical_dod = None
		return

	def clear_Neumann(self):
		" Clears Neumann boundaries "
		self._thermalNeumann = None
		self._mechanicalNeumann = None
		return
	
	def add_thermal_IBC(self, table=None, temperature=0.0):
		"This function is first tentative of adding constant boundary conditions "
		if self._thermalDirichlet is None: raise Warning('Please define first total Dirichlet boundaries')
		if np.any(self._thermalDirichlet-table<0): raise Warning('Table is not well defined')
		if not np.any(table>0): raise Warning('At least one blocked face is needed')
		dod = self.get_boundaryNodes(table=table)[1][0]
		self._ThermalDirichletBound[dod] += temperature*np.ones(len(dod))
		return 

	def get_thermal_IBC(self): 
		if self._thermalDirichlet is None: raise Warning('Please define first total Dirichlet boundaries')
		return self._ThermalDirichletBound[self._thermal_dod]
	