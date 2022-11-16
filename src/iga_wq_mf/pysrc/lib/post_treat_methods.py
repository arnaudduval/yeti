"""
.. This module contains functions used to run a simulation and plot the results
.. Joaquin Cornejo
"""

from .__init__ import *

# My libraries
from .create_geomdl import geomdlModel
from .fortran_mf_iga import fortran_mf_iga
from .fortran_mf_wq import fortran_mf_wq

def plot_iterative_solver(filename, inputs:dict, extension='.pdf', threshold=1.e-12, pLegend=True):
	
	savename = filename.split('.')[0] + extension
	residue = inputs["resPCG"]
	method_list = inputs["Methods"]

	new_method_list = []
	for pcgmethod in method_list:
		if pcgmethod   == "WP" : new_method_list.append('w.o. preconditioner')
		elif pcgmethod == "C"  : new_method_list.append('Classic FD method')
		elif pcgmethod == "TDS": new_method_list.append('Literature + scaling') 
		elif pcgmethod == "TDC": new_method_list.append('Literature') 
		elif pcgmethod == "JMC": new_method_list.append('This work') 
		elif pcgmethod == "JMS": new_method_list.append('This work + scaling')
	
	markers = ['o', 'v', 's', 'X', '+', 'p']

	# Set figure parameters
	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))
	for i, pcgmethod in enumerate(new_method_list):
		residue_method = np.asarray(residue[i])
		residue_method = residue_method[residue_method>threshold]
		ax.semilogy(np.arange(len(residue_method)), residue_method, '-',
					label=pcgmethod, marker=markers[i])

	ax.set_ylim(top=10.0, bottom=threshold)
	ax.set_xlabel('Number of iterations of BiCGSTAB solver')
	ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')

	if pLegend: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	ax.legend(loc=0)
	fig.tight_layout()
	fig.savefig(savename)

	return

class ThermalSimulation():

	def __init__(self, inputs:dict, folder=None):

		# Get data to run simulation
		self._degree = inputs.get('degree', 2)
		self._cuts = inputs.get('cuts', 2)
		self._geoName = inputs.get('name', 'TR')
		self._isIGA = inputs.get('isIGA', False)
		self._funPowDen = inputs.get('funPowDen', None)
		self._funTemp = inputs.get('funTemp', None)   
		self._method_list = inputs.get('IterMethods', ['WP'])
		self._isOnlyIter = inputs.get('isOnlyIter', True)
		self._thermalModel = None
		self._nbIter = 100
		self._threshold = 1e-15

		# Get filename
		filename = self.get_filename()
		self._filename = folder + filename

		return

	def setup_PCGparameters(self, nbIter=100, threshold=1e-15):
		self._nbIter = nbIter
		self._threshold = threshold
		return
	
	def get_filename(self):
		" Set filename using input information "
		# Get text file name
		filename = (self._geoName 
					+ '_p_' + str(self._degree) 
					+ '_nbel_' + str(2**self._cuts)
		)
		if self._isIGA: filename += '_IGAG'
		else: filename += '_IGAWQ'
		filename += '.txt'

		return filename

	def create_geometryModel(self):
		" Create model using YETI refinement "

		geometry = {'degree': self._degree*np.ones(3, dtype=int)}
		modelGeo = geomdlModel(self._geoName, **geometry)
		modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=self._cuts*np.ones(3, dtype=int))

		return modelIGA

	def create_thermalModel(self, material=None, Dirichlet=None):
		" Creates thermal model in IGA-Galerkin or IGA-WQ approach "

		modelIGA = self.create_geometryModel()
		if self._isIGA: thermalModel = fortran_mf_iga(modelIGA, material=material, Dirichlet=Dirichlet)
		else:           thermalModel = fortran_mf_wq(modelIGA, material=material, Dirichlet=Dirichlet)

		return thermalModel

	def compute_force_SM(self, model:fortran_mf_wq, funPowDen=None, funTemp=None):
		""" Compute force vector b = Fn - Knd.ud where F is source vector and K is conductivity matrix. 
			This equation is used in substitution method (SM)
		"""
		
		dof = model._thermal_dof
		dod = model._thermal_dod

		if funTemp is not None:  
			ud = model.interpolate_ControlPoints(funfield=funTemp)[dod]
			u = np.zeros(model._nb_ctrlpts_total); u[dod] = ud
			Fn = model.eval_source_vector(funPowDen, indi=dof) 
			Knd_ud = model.eval_Ku(u)[0]
			Fn -= Knd_ud[dof]
		else:
			ud = None
			Fn = model.eval_source_vector(funPowDen, indi=dof)
			
		return Fn, ud
	
	def run_iterative_solver(self, model:fortran_mf_wq, b=None, nbIterPCG=100, threshold=1e-15, methodPCG='WP'):
		" Solve steady heat problems using iterative solver "
		
		start = time.process_time()
		if model._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
		un, residue = model.MFsteadyHeat(b, nbIterPCG=nbIterPCG, threshold=threshold, methodPCG=methodPCG)
		stop = time.process_time()
		time_t = stop - start 

		return un, residue, time_t

	def run_direct_solver(self, model:fortran_mf_wq, b=None):
		" Solve steady heat problems using direct solver "

		start = time.process_time()
		if model._thermalDirichlet is None: raise Warning('Ill conditionned. It needs Dirichlet conditions')
		A = model.eval_conductivity_matrix(indi=model._thermal_dof, indj=model._thermal_dof)
		stop = time.process_time()
		time_assembly = stop - start 

		# Solve system
		start = time.process_time()
		un = sp.linalg.spsolve(A, b)
		stop = time.process_time()
		time_solver = stop - start 

		return un, time_assembly, time_solver

	def run_simulation(self, material=None, Dirichlet=None):
		" Runs simulation using given input information "

		doDirect = True
		if self._isOnlyIter: doDirect = False
		timeAssembly, timeDirect = 0.0, 0.0
		timeNoIter, timeIter  = 0.0, 0.0
		resPCG = np.zeros(self._nbIter+1)

		self._thermalModel = self.create_thermalModel(material=material, Dirichlet=Dirichlet)
		b = self.compute_force_SM(self._thermalModel, self._funPowDen, self._funTemp)[0]

		# Run direct solver 
		if doDirect: timeAssembly, timeDirect = self.run_direct_solver(self._thermalModel, b)[1:]

		# Run iterative methods
		timeNoIter = []
		for itmethod in self._method_list:
			time_temp = self.run_iterative_solver(self._thermalModel, b=b, methodPCG=itmethod, 
													nbIterPCG=0, threshold=self._threshold)[-1]
			timeNoIter.append(time_temp)

		timeIter, resPCG = [], []
		for itmethod in self._method_list:
			residue_t, time_temp = self.run_iterative_solver(self._thermalModel, b=b, methodPCG=itmethod, 
													nbIterPCG=self._nbIter, threshold=self._threshold)[1:]
			timeIter.append(time_temp)
			resPCG.append(residue_t)
				
		# Write file
		output = {'Methods': self._method_list, 'TimeAssembly': timeAssembly, 'TimeDirect': timeDirect, 
				'TimeNoIter':timeNoIter, 'TimeIter': timeIter, 'resPCG': resPCG}
		self.write_text_file(output)

		return 

	def write_text_file(self, inputs:dict): 
		" Writes and exports simulation data in .txt file "
		
		iterMethods  = inputs['Methods']
		timeAssembly = inputs['TimeAssembly']
		timeDirect   = inputs['TimeDirect']
		timeIter     = inputs['TimeIter']
		timeNoIter   = inputs['TimeNoIter']
		resPCG       = inputs['resPCG']
		
		with open(self._filename, 'w') as f:
			f.write('** RESULTS **\n')
			f.write('** Direct solver\n')
			f.write('*Time assembly\n')
			f.write('{:E}\n'.format(timeAssembly))
			f.write('*Time direct\n')
			f.write('{:E}\n'.format(timeDirect))
			f.write('** Iterative solver ' + ','.join([item for item in iterMethods]) + '\n')
			f.write('** Number of iterations ' + '{:d}\n'.format(self._nbIter))

			for i, method in enumerate(iterMethods):
				f.write('**' + method + '\n')
				f.write('*Time prepare ' + method +'\n')
				f.write('{:E}\n'.format(timeNoIter[i]))
				f.write('*Time iter ' + method +'\n')
				f.write('{:E}\n'.format(timeIter[i]))
				f.write('*Residue ' + method + '\n')
				f.writelines(['{:E}'.format(res) + '\n'
								for res in resPCG[i]]) 
		return

class SimulationData(): 

	def __init__(self, filename=None): 

		self._filename = filename
		
		# Get important data from title
		full_path = os.path.realpath(self._filename)
		name = os.path.basename(full_path).split('.')[0]
		self._interpret_title(name)

		self._dataSimulation = self.get_infos_simulation()
	
		return

	def _interpret_title(self, name):
		" Gets important information from title "

		# Split string
		ls = name.split('_')

		# Get data
		self._geoCase = ls[0]
		self._degree = int(ls[2])
		self._cuts = int(np.log2(int(ls[4])))

		if ls[5] == 'IGAG': self._isIGA = True
		else: self._isIGA = False

		return

	def get_infos_simulation(self):
		" Gets the information of the simulation "

		lines = self._get_cleanLines()
		iterMethods, nbIter = self._read_iter_methods(lines)

		time_noiter, time_itersolver, resPCG = [], [], []
		for itmethod in iterMethods:
			tp = self._read_time_preparation(lines, itmethod)
			ti = self._read_time_iterations(lines, itmethod)
			res = self._read_residue(lines, itmethod, nbIter)

			time_noiter.append(tp); time_itersolver.append(ti)
			resPCG.append(res)

		output = {'Methods': iterMethods, 'TimeNoIter':time_noiter, 
					'TimeIter': time_itersolver, 'resPCG': resPCG}
		return output

	def _get_cleanLines(self):
		" Gets the lines from the file that could have important information "

		inputFile  = open(self._filename, 'r')
		lines      = inputFile.readlines()
		linesClean = []
		theLine    = ''
		# Joining lines ending with ',' in a new list of lines
		# This allow to merge definitions written on more than a single line
		for i in range(0, len(lines)):
			words = lines[i].rstrip().split(',')
			# Removing trailing spaces
			lastWord = words[-1].split()
			if(len(lastWord) == 0):
				theLine = theLine + lines[i]
				# Removing '\n' character
				theLine = theLine.rstrip()
			else:
				theLine = theLine + lines[i]
				linesClean.append(theLine.rstrip())
				theLine = ''
		inputFile.close()

		return linesClean

	def _get_num_line(self, lines, str2find):
		i=0
		for line in lines:
			if line.startswith(str2find): break
			i += 1
		if i==len(lines):
			print("Error: keyword " + str2find + " is missing in data file.")
		return i

	def _read_time_assembly(self, lines):
		i = self._get_num_line(lines,'*Time assembly')
		time = float(lines[i+1])
		return time

	def _read_time_solve(self, lines): 
		i = self._get_num_line(lines,'*Time direct')
		time = float(lines[i+1])
		return time

	def _read_iter_methods(self, lines):
		i = self._get_num_line(lines,'** Iterative solver')
		ls = lines[i].split(' ')[-1]
		methods = ls.split(',')
		ls = lines[i+1].split(' ')[-1]
		nbIter = int(ls)
		return methods, nbIter

	def _read_time_preparation(self, lines, method):
		i = self._get_num_line(lines,'*Time prepare ' + method)
		time = float(lines[i+1])
		return time

	def _read_time_iterations(self, lines, method):
		i = self._get_num_line(lines,'*Time iter ' + method)
		time = float(lines[i+1])
		return time
	
	def _read_residue(self, lines, method, length=100):
		i = self._get_num_line(lines,'*Residue ' + method) + 1
		lines_data = lines[i:i+length+1]

		residue = []
		for line in lines_data:
			ls = line.split(',')
			residue.append(float(ls[0]))

		return residue
	