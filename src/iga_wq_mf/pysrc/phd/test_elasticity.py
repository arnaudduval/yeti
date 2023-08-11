from pysrc.lib.__init__ import *
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_job import mechaproblem
from pysrc.lib.lib_simulation import decoder

def forceSurfFun(P:list):
	ref  = np.array([0.0, 4e1])
	prop = np.zeros((2, np.size(P, axis=1)))
	for i in range(2): prop[i, :] = ref[i] 
	return prop

class simulate():

	def __init__(self, simuArgs:dict):

		self._name 		= simuArgs.get('name', '').lower()
		self._degree 	= simuArgs.get('degree', 2)
		self._nbcuts 	= simuArgs.get('nb_refinementByDirection', 2)
		self._nbctrlpts = int(self._degree + 2**self._nbcuts) 
		
		self._isGaussQuad = simuArgs.get('isGauss', False)
		self._funSurf     = simuArgs.get('funSurf', None)
		self._iterMethods = simuArgs.get('IterMethods', ['WP'])
		
		self._filename = simuArgs.get('folder', './') + self.__make_filename() 
		
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
	
	def __run_iterativeSolver(self, problem:mechaproblem, Fext):
		" Solve steady heat problems using iterative solver "
		start = time.process_time()
		sol, residue, _ = problem.solveElasticityProblemFT(Fext)
		stop = time.process_time()
		time_t = stop - start 
		return sol, residue, time_t

	def simulate(self, material:mechamat, boundary:boundaryCondition, overwrite=True):
		" Runs simulation using given input information "

		geoArgs  = {'name':self._name, 
					'nb_refinementByDirection': self._nbcuts*np.ones(3, dtype=int),
					'degree': self._degree*np.ones(3, dtype=int)}
		modelgeo = Geomdl(geoArgs)
		modelIGA = modelgeo.getIGAParametrization()
		if self._isGaussQuad: quadArgs = {'quadrule': 'iga'}
		else:                 quadArgs = {'quadrule': 'wq'}
		modelPhy = part(modelIGA, quadArgs)

		problem = mechaproblem(material, modelPhy, boundary)
		Fext    = problem.compute_surfForce(self._funSurf, nbFacePosition=1)[0]

		# Run iterative methods
		timeNoIter = []; problem.addSolverConstraints({'nbIterationsPCG':0})
		for im in self._iterMethods:
			problem._methodPCG = im
			time_temp = self.__run_iterativeSolver(problem, Fext=Fext)[-1]
			timeNoIter.append(time_temp)

		timeIter, resPCG = [], []; problem.addSolverConstraints({'nbIterationsPCG':100, 'PCGThreshold':1e-12})
		for im in self._iterMethods:
			problem._methodPCG = im
			un, residue_t, time_temp = self.__run_iterativeSolver(problem, Fext=Fext)
			timeIter.append(time_temp)
			resPCG.append(residue_t)
				
		# Write file
		output = {'nbIterPCG': problem._nbIterPCG, 'iterMethods': self._iterMethods, 'timeNoIter':timeNoIter, 'timeIter': timeIter, 'resPCG': resPCG}
		if overwrite: self.write_resultsFile(output)

		return un

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

dataExist   = False
degree_list = np.arange(6, 7)
cuts_list   = np.arange(6, 7)
name_list   = ['qa', 'sq']
IterMethods = ['WP', 'C', 'JMC', 'TDC']
matArgs     = {'elastic_modulus':1e3, 'elastic_limit':1e10, 'poisson_ratio':0.3}

for cuts in cuts_list:
	for degree in degree_list:
		for name in name_list: 

			inputs = {'degree': degree, 'nb_refinementByDirection': cuts, 'name': name, 'isGauss': False, 
					'funSurf': forceSurfFun, 'IterMethods': IterMethods, 'folder': folder}
			simulation = simulate(inputs)  
			
			if not dataExist:
				mat = mechamat(matArgs)
				nbctrlpts = np.zeros(3, dtype=int); nbctrlpts[:2] = simulation._nbctrlpts
				boundary = boundaryCondition(nbctrlpts)
				table = np.zeros((2, 2, 2), dtype=int)
				table[0, 0, 0] = 1
				table[1, 0, 1] = 1
				boundary.add_DirichletDisplacement(table=table)
				simulation.simulate(material=mat, boundary=boundary, overwrite=True)

			else :
				simuOutput = decoder(simulation._filename)
				simuOutput.plot_results(extension='.pdf', plotLegend=True)
				