from .__init__ import *
from .lib_material import thermomat
from .lib_geomdl import Geomdl
from .lib_part import part
from .lib_boundary import boundaryCondition
from .lib_job import heatproblem

class encoder():

	def __init__(self, simuArgs:dict):

		self._name 		= simuArgs.get('name', '').lower()
		self._degree 	= simuArgs.get('degree', 2)
		self._nbcuts 	= simuArgs.get('nb_refinementByDirection', 2)
		self._nbctrlpts = int(self._degree + 2**self._nbcuts) 
		
		self._isGaussQuad = simuArgs.get('isGauss', False)
		self._funPowDen   = simuArgs.get('funPowerDensity', None)
		self._funTemp     = simuArgs.get('funTemperature', None)   
		self._iterMethods = simuArgs.get('IterMethods', ['WP'])
		
		filename = self.__make_filename()
		self._filename = simuArgs.get('folder', './') + filename 
		self._part = None
		
		return
	
	def __make_filename(self):
		" Make up filename using input information "
		# Get text file name
		filename = (self._name 
					+ '_p_' + str(self._degree) 
					+ '_nbel_' + str(2**self._nbcuts)
		)
		if self._isGaussQuad: filename += '_IGAG'
		else: filename += '_IGAWQ'
		filename += '.txt'
		return filename
	
	def __create_model(self):
		# Create geometry 
		geoArgs   = {'name':self._name, 
					'nb_refinementByDirection': self._nbcuts*np.ones(3, dtype=int),
					'degree': self._degree*np.ones(3, dtype=int)}
		modelgeo = Geomdl(geoArgs)
		modelIGA = modelgeo.getIGAParametrization()
		if self._isGaussQuad: quadArgs = {'quadrule': 'iga'}
		else:                 quadArgs = {'quadrule': 'wq'}
		self._part = part(modelIGA, quadArgs)
		return 

	def write_resultsFile(self, inputs:dict): 
		" Writes and exports simulation data in .txt file "
		
		nbIterPCG   = inputs['nbIterPCG']
		iterMethods = inputs['iterMethods']
		timeNoIter  = inputs['timeNoIter']
		timeIter    = inputs['timeIter']
		residuePCG  = inputs['resPCG']
		
		with open(self._filename, 'w') as f:
			f.write('** RESULTS **\n')
			f.write('** Iterative solver ' + ','.join([item for item in iterMethods]) + '\n')
			f.write('** Number of iterations ' + '{:d}\n'.format(nbIterPCG))

			for i, method in enumerate(iterMethods):
				f.write('**' + method + '\n')
				f.write('*Time prepare ' + method +'\n')
				f.write('{:E}\n'.format(timeNoIter[i]))
				f.write('*Time iter ' + method +'\n')
				f.write('{:E}\n'.format(timeIter[i]))
				f.write('*Residue ' + method + '\n')
				f.writelines(['{:E}'.format(res) + '\n'
								for res in residuePCG[i]]) 
		return
	
	def __eval_heatForce(self, problem:heatproblem, funPowDen=None, funTemp=None):
		""" Compute force vector b = Fn - Knd.ud where F is source vector and K is conductivity matrix. 
			This equation is used in substitution method (SM)
		"""
		
		dof = problem.boundary.thdof
		dod = problem.boundary.thdod
		nbctrlpts_total = problem.part.nbctrlpts_total

		if funTemp is not None:  
			ud = problem.L2projectionCtrlpts(funfield=funTemp)[dod]
			u  = np.zeros(nbctrlpts_total); u[dod] = ud
			Fn = problem.eval_volForce(funPowDen, indi=dof) 
			Knd_ud = problem.eval_mfConductivity(u, table=np.zeros((3, 2), dtype=bool), args=self._part.qpPhy)
			Fn    -= Knd_ud[dof]
		else:
			Fn = problem.eval_volForce(funPowDen, indi=dof)
		
		return Fn
	
	def __run_iterativeSolver(self, problem:heatproblem, b):
		" Solve steady heat problems using iterative solver "
		start = time.process_time()
		un, residue = problem.solveSteadyHeatProblemFT(b)
		stop = time.process_time()
		time_t = stop - start 
		return un, residue, time_t

	def simulate(self, material:thermomat, boundary:boundaryCondition, overwrite=True):
		" Runs simulation using given input information "

		if self._part is None: self.__create_model()
		problem = heatproblem(material, self._part, boundary)
		Fn      = self.__eval_heatForce(problem, self._funPowDen, self._funTemp)

		# Run iterative methods
		timeNoIter = []; problem.addSolverConstraints({'nbIterationsPCG':0})
		for im in self._iterMethods:
			problem._methodPCG = im
			time_temp = self.__run_iterativeSolver(problem, b=Fn)[-1]
			timeNoIter.append(time_temp)

		timeIter, resPCG = [], []; problem.addSolverConstraints({})
		for im in self._iterMethods:
			problem._methodPCG = im
			un, residue_t, time_temp = self.__run_iterativeSolver(problem, b=Fn)
			timeIter.append(time_temp)
			resPCG.append(residue_t)
				
		# Write file
		output = {'nbIterPCG': problem._nbIterPCG, 'iterMethods': self._iterMethods, 'timeNoIter':timeNoIter, 'timeIter': timeIter, 'resPCG': resPCG}
		if overwrite: self.write_resultsFile(output)

		return un

class decoder(): 

	def __init__(self, filename=None): 

		if filename is None: raise Warning('Not possible. Insert a valid path')
		self._filename = filename
		self._dataSimulation = {}
		full_path = os.path.realpath(filename)
		title     = os.path.basename(full_path).split('.')[0]
		self.__decrypt_title(title)
		self._dataSimulation = self.__getInfo()
	
		return

	def __decrypt_title(self, name):
		" Gets important information from title "

		# Split string
		ls = name.split('_')

		# Get data
		self._geoname = ls[0]
		self._degree  = int(ls[2])
		self._cuts    = int(np.log2(int(ls[4])))

		if ls[5] == 'IGAG': self._isGaussQuad = True
		else: self._isGaussQuad = False

		return

	def __getInfo(self):
		" Gets the information of the simulation "

		lines = self.__get_cleanLines(self._filename)
		iterMethods, nbiter = self.__read_iterMethods(lines)

		time_noiter, time_itersolver, resPCG = [], [], []
		for itmethod in iterMethods:
			tp = self.__read_timePreparation(lines, itmethod)
			ti = self.__read_timeIterations(lines, itmethod)
			res = self.__read_residue(lines, itmethod, nbiter)
			time_noiter.append(tp); time_itersolver.append(ti); resPCG.append(res)

		output = {'methods': iterMethods, 'timeNoIter':time_noiter, 
				'timeIter': time_itersolver, 'resPCG': resPCG}
		return output

	def __get_cleanLines(self, filename):
		" Gets the lines from the file that could have important information "

		inputFile  = open(filename, 'r')
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

	def __get_numLine(self, lines, str2find):
		i=0
		for line in lines:
			if line.startswith(str2find): break
			i += 1
		if i==len(lines):
			print('Error: keyword ' + str2find + ' is missing in data file.')
		return i

	def __read_iterMethods(self, lines):
		i = self.__get_numLine(lines,'** Iterative solver')
		ls = lines[i].split(' ')[-1]
		methods = ls.split(',')
		ls = lines[i+1].split(' ')[-1]
		nbIter = int(ls)
		return methods, nbIter

	def __read_timePreparation(self, lines, method):
		i = self.__get_numLine(lines,'*Time prepare ' + method)
		time = float(lines[i+1])
		return time

	def __read_timeIterations(self, lines, method):
		i = self.__get_numLine(lines,'*Time iter ' + method)
		time = float(lines[i+1])
		return time
	
	def __read_residue(self, lines, method, length=100):
		i = self.__get_numLine(lines,'*Residue ' + method) + 1
		lines_data = lines[i:i+length+1]

		residue = []
		for line in lines_data:
			ls = line.split(',')
			residue.append(float(ls[0]))

		return residue
	
	def plot_results(self, extension='.pdf', threshold=1.e-12, plotLegend=True):

		savename = self._filename.split('.')[0] + extension
		residue  = self._dataSimulation.get('resPCG')
		method_list = self._dataSimulation.get('methods')

		new_method_list = []
		for pcgmethod in method_list:
			if pcgmethod   == 'WP' : new_method_list.append('w.o. preconditioner')
			elif pcgmethod == 'C'  : new_method_list.append('Classic FD method')
			elif pcgmethod == 'TDS': new_method_list.append('Literature + scaling') 
			elif pcgmethod == 'TDC': new_method_list.append('Literature') 
			elif pcgmethod == 'JMC': new_method_list.append('This work') 
			elif pcgmethod == 'JMS': new_method_list.append('This work + scaling')
		
		markers = ['o', 'v', 's', 'X', '+', 'p']

		# Set figure parameters
		fig, ax = plt.subplots(nrows=1, ncols=1)
		for i, pcgmethod in enumerate(new_method_list):
			residue_method = np.asarray(residue[i])
			residue_method = residue_method[residue_method>threshold]
			ax.semilogy(np.arange(len(residue_method)), residue_method, '-',
						label=pcgmethod, marker=markers[i])

		ax.set_ylim(top=10.0, bottom=threshold)
		ax.set_xlabel('Number of iterations of BiCGSTAB solver')
		ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')

		# if plotLegend: ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		ax.legend(loc='best')
		fig.tight_layout()
		fig.savefig(savename)

		return