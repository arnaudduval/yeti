from pysrc.lib.__init__ import *
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_job import mechaproblem
from pysrc.phd.lib_simulation import decoder, simulate

def forceSurf_infPlate(P:list):
	Tx, a = 1.0, 1.0
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	b = a**2/r_square
	theta = np.arcsin(y/np.sqrt(r_square))

	F = np.zeros((2, nnz))
	F[0, :] = Tx/2*(2*np.cos(theta) - b*(2*np.cos(theta) + 3*np.cos(3*theta)) + 3*b**2*np.cos(3*theta))
	F[1, :] = Tx/2*3*np.sin(3*theta)*(b**2 - b)
	return F

def exactDisplacement_infPlate(P:list):
	Tx, a = 1.0, 1.0
	E, nu = 1e3, 0.3
	x = P[0, :]; y = P[1, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	theta = np.arcsin(y/np.sqrt(r_square))
	b = a**2/r_square # Already squared
	c = Tx*(1.0 + nu)*np.sqrt(r_square)/(2*E)

	disp = np.zeros((2, nnz))
	disp[0, :] = c*(2*(1-nu)*np.cos(theta) + b*(4*(1-nu)*np.cos(theta) + np.cos(3*theta)) - b**2*np.cos(3*theta))
	disp[1, :] = c*(-2*nu*np.sin(theta) + b*(2*(-1 + 2*nu)*np.sin(theta) + np.sin(3*theta)) - b**2*np.sin(3*theta))
	
	return disp

class mechasimulate(simulate):
	def __init__(self, simuArgs:dict):
		super().__init__(simuArgs)
		return

	def simulate(self, material:mechamat, boundary:boundaryCondition, overwrite=True):
		" Runs simulation using given input information "

		blockPrint()
		geoArgs  = {'name':self._name, 
					'nb_refinementByDirection': self._nbcuts*np.ones(3, dtype=int),
					'degree': self._degree*np.ones(3, dtype=int)}
		modelgeo = Geomdl(geoArgs)
		modelIGA = modelgeo.getIGAParametrization()
		if self._isGaussQuad: quadArgs = {'quadrule': 'iga'}
		else:                 quadArgs = {'quadrule': 'wq'}
		modelPhy = part(modelIGA, quadArgs)

		problem = mechaproblem(material, modelPhy, boundary)
		Fext    = problem.compute_surfForce(self._funforce, nbFacePosition=1)[0]

		# Run iterative methods
		timeNoIter = []; problem.addSolverConstraints({'nbIterationsPCG':0})
		for im in self._iterMethods:
			problem._methodPCG = im
			time_temp = self.run_iterativeSolver(problem, Fext=Fext)[-1]
			timeNoIter.append(time_temp)
		enablePrint()
		timeIter, resPCG = [], []; problem.addSolverConstraints({'nbIterationsPCG':100, 'PCGThreshold':1e-12})
		print(self._degree, self._nbcuts)
		for im in self._iterMethods:
			problem._methodPCG = im
			un, residue_t, time_temp = self.run_iterativeSolver(problem, Fext=Fext)
			timeIter.append(time_temp)
			resPCG.append(residue_t)
			print(im, time_temp, len(residue_t[residue_t>0.0]))
		print('--')
				
		# Write file
		output = {'nbIterPCG': problem._nbIterPCG, 'iterMethods': self._iterMethods, 'timeNoIter':timeNoIter, 'timeIter': timeIter, 'resPCG': resPCG}
		if overwrite: self.write_resultsFile(output)

		return un

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

dataExist   = False 
degree_list = np.arange(4, 7)
cuts_list   = np.arange(5, 8)
name_list   = ['qa']
IterMethods = ['JMC']
matArgs     = {'elastic_modulus':1e3, 'elastic_limit':1e10, 'poisson_ratio':0.3}

for cuts in cuts_list:
	for degree in degree_list:
		for name in name_list: 

			inputs = {'degree': degree, 'nb_refinementByDirection': cuts, 'name': name, 'isGauss': False, 
					'funforce': forceSurf_infPlate, 'IterMethods': IterMethods, 'folder': folder}
			simulation = mechasimulate(inputs)  
			
			if not dataExist:
				mat = mechamat(matArgs)
				nbctrlpts = np.ones(3, dtype=int); nbctrlpts[:2] = simulation._nbctrlpts
				boundary = boundaryCondition(nbctrlpts)
				table = np.zeros((2, 2, 2), dtype=int)
				table[0, 0, 0] = 1; table[1, 0, 1] = 1
				boundary.add_DirichletDisplacement(table=table)
				simulation.simulate(material=mat, boundary=boundary, overwrite=True)

			else :
				simuOutput = decoder(simulation._filename)
				simuOutput.plot_results(extension='_EL.pdf', plotLegend=False)
				