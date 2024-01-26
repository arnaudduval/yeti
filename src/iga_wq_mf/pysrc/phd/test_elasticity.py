from pysrc.lib.__init__ import *
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_job import mechaproblem
from pysrc.phd.lib_simulation import decoder, simulate

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
TRACTION, RINT, REXT = 1.0, 1.0, 2
YOUNG, POISSON = 1e3, 0.1
GEONAME = 'QA'
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':2e10, 'poisson_ratio': POISSON, 
			'plasticLaw': {'Isoname':'linear', 'Eiso':YOUNG/10}}

DEGREE_LIST = np.arange(5, 6)
CUTS_LIST   = np.arange(6, 7)
ITERMETHODS = ['WP', 'C', 'JMC']
dataExist   = True

def forceSurf_infPlate(P:list):
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = RINT**2/r_square # Already squared
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = TRACTION/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = TRACTION/2*3*np.sin(3*theta)*(b**2 - b)
	return F

class mechasimulate(simulate):
	def __init__(self, simuArgs:dict):
		super().__init__(simuArgs)
		return

	def simulate(self, material:mechamat, boundary:boundaryCondition, overwrite=True):
		" Runs simulation using given input information "

		blockPrint()
		geoArgs  = {'name':self._name, 
					'nb_refinementByDirection': self._nbcuts*np.ones(3, dtype=int),
					'degree': self._degree*np.ones(3, dtype=int),
					'extra':{'Rin':RINT, 'Rex':REXT}
					}
		modelgeo = Geomdl(geoArgs)
		modelIGA = modelgeo.getIGAParametrization()
		if self._isGaussQuad: quadArgs = {'quadrule': 'iga'}
		else:                 quadArgs = {'quadrule': 'wq'}
		modelPhy = part(modelIGA, quadArgs)

		problem = mechaproblem(material, modelPhy, boundary)
		Fext    = problem.compute_surfForce(self._funforce, nbFacePosition=1)[0]

		# Run iterative methods	
		timeNoIter = []; problem.addSolverConstraints({'nIterKrylov':0})
		for im in self._iterMethods:
			problem._linPreCond = im
			time_temp = self.run_iterativeSolver(problem, Fext=Fext)[-1]
			timeNoIter.append(time_temp)
		enablePrint()
		timeIter, resPCG = [], []; problem.addSolverConstraints({'nIterKrylov':100, 'thresholdKrylov':1e-12})
		print(self._degree, self._nbcuts)
		for im in self._iterMethods:
			problem._linPreCond = im
			un, residue_t, time_temp = self.run_iterativeSolver(problem, Fext=Fext)
			timeIter.append(time_temp)
			resPCG.append(residue_t)
			print(im, time_temp, len(residue_t[residue_t>0.0]))
		print('--')
				
		# Write file
		output = {'nbIterPCG': problem._itersLin, 'iterMethods': self._iterMethods, 'timeNoIter':timeNoIter, 'timeIter': timeIter, 'resPCG': resPCG}
		if overwrite: self.write_resultsFile(output)

		return un

for cuts in CUTS_LIST:
	for degree in DEGREE_LIST:

		inputs = {'degree': degree, 'nb_refinementByDirection': cuts, 'name': GEONAME, 'isGauss': False, 
				'funforce': forceSurf_infPlate, 'IterMethods': ITERMETHODS, 'folder': folder}
		simulation = mechasimulate(inputs)  
		
		if not dataExist:
			mat = mechamat(MATARGS)
			nbctrlpts = np.ones(3, dtype=int); nbctrlpts[:2] = simulation._nbctrlpts
			boundary = boundaryCondition(nbctrlpts)
			table = np.zeros((2, 2, 2), dtype=int)
			table[1, 1, 0] = 1; table[1, 1, 1] = 1
			boundary.add_DirichletDisplacement(table=table)
			simulation.simulate(material=mat, boundary=boundary, overwrite=True)

		else :
			simuOutput = decoder(simulation._filename)
			simuOutput.plot_results(extension='_EL.pdf', plotLegend=True)
			