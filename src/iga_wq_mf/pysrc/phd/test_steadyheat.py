from pysrc.lib.__init__ import *
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_job import heatproblem
from pysrc.phd.lib_simulation import decoder, simulate

def conductivityProperty(P:list):
	dimen = np.size(P, axis=0)
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((dimen, dimen, np.size(P, axis=1)))
	for i in range(dimen): 
		for j in range(dimen):
			Kprop[i, j, :] = Kref[i, j] 
	# x = P[0, :]; y = P[1, :]; z = P[2, :]
	# Kprop[0, 0, :] += 0.75*np.cos(np.pi*y)
	# Kprop[1, 1, :] += 2*np.exp(-(z-0.5)**2)
	# if dimen > 2: Kprop[2, 2, :] += 2.5*np.cos(np.pi*x)**2
	return Kprop 

def powerDensity_cube(P: list):
	""" u = sin(pi*x)*sin(pi*y)*sin(pi*z)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]

	# Anisotropy
	f = (6*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 
	- (np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*z)*np.sin(np.pi*y))/5 
	- (np.pi**2*np.cos(np.pi*y)*np.cos(np.pi*z)*np.sin(np.pi*x))/2 
	- np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(np.pi*z)
	)

	return f

def powerDensity_prism(P: list):
	""" u = (-5*x+6*y+45)*(5*x+6*y-45)*x*(x-6)*sin(pi*z)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]

	# Anisotropy
	f = (4*x*np.sin(np.pi*z)*(5*x + 6*y - 45) 
	- 94*x*np.sin(np.pi*z)*(x - 6) 
	- 16*x*np.sin(np.pi*z)*(6*y - 5*x + 45) 
	- 2*np.sin(np.pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
	- 16*np.sin(np.pi*z)*(x - 6)*(6*y - 5*x + 45) 
	+ 4*np.sin(np.pi*z)*(x - 6)*(5*x + 6*y - 45) 
	- (np.pi*np.cos(np.pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/5 
	- 4*x*np.pi*np.cos(np.pi*z)*(x - 6)*(6*y - 5*x + 45) 
	- 2*x*np.pi*np.cos(np.pi*z)*(x - 6)*(5*x + 6*y - 45) 
	- (x*np.pi*np.cos(np.pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/5 
	+ 3*x*np.pi**2*np.sin(np.pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
	)/1000
	
	return f

def powerDensity_thickRing(P: list):
	""" u = sin(5*pi*x)*sin(5*pi*y)*sin(5*pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]
	z = P[2, :] 

	# Anisotropy
	f = (150*np.pi**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
	- 16*y**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z) 
	- 6*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 6*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 8*x*y*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z) 
	- 8*x**2*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z) 
	- 25*np.pi**2*np.cos(5*np.pi*x)*np.cos(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- 5*np.pi**2*np.cos(5*np.pi*x)*np.cos(5*np.pi*z)*np.sin(5*np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- (25*np.pi**2*np.cos(5*np.pi*y)*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4))/2 
	- 20*x*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 10*x*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 2*x*np.pi*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*(x**2 + y**2 - 1) 
	- 20*x*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 10*x*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 2*x*np.pi*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*(x**2 + y**2 - 4) 
	- 10*y*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 40*y*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 1) 
	- 5*y*np.pi*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*(x**2 + y**2 - 1) 
	- 10*y*np.pi*np.cos(5*np.pi*x)*np.sin(5*np.pi*y)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 40*y*np.pi*np.cos(5*np.pi*y)*np.sin(5*np.pi*x)*np.sin(5*np.pi*z)*(x**2 + y**2 - 4) 
	- 5*y*np.pi*np.cos(5*np.pi*z)*np.sin(5*np.pi*x)*np.sin(5*np.pi*y)*(x**2 + y**2 - 4)
	)

	return f

def powerDensity_quartCircle(P: list):
	""" u = sin(pi*x)*sin(pi*y)*(x**2+y**2-1)*(x**2+y**2-16)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]

	f = (3*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 16) 
	- 16*y**2*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 16) 
	- 8*x*y*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 16) 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 16) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 16) 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 16) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 16) 
	- 8*x**2*np.sin(np.pi*x)*np.sin(np.pi*y)
	)

	return f

class heatsimulate(simulate):

	def __init__(self, simuArgs:dict):
		super().__init__(simuArgs)
		return

	def simulate(self, material:heatmat, boundary:boundaryCondition, overwrite=True):
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

		problem = heatproblem(material, modelPhy, boundary)
		Fext    = problem.compute_volForce(self._funforce)

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

for cuts in cuts_list:
	for degree in degree_list:
		for name in name_list: 

			if name   == 'cb' : funpow = powerDensity_cube 
			elif name == 'vb' : funpow = powerDensity_prism 
			elif name == 'tr' : funpow = powerDensity_thickRing 
			elif name == 'qa' : funpow = powerDensity_quartCircle

			inputs = {'degree': degree, 'nb_refinementByDirection': cuts, 'name': name, 'isGauss': False, 
					'funforce': funpow, 'IterMethods': IterMethods, 'folder': folder}
			simulation = heatsimulate(inputs)  
			
			if not dataExist:
				mat = heatmat()
				mat.addConductivity(conductivityProperty, isIsotropic=False)
				nbctrlpts = simulation._nbctrlpts*np.ones(3, dtype=int)
				if name == 'qa': nbctrlpts[-1] = 1
				boundary = boundaryCondition(nbctrlpts)
				boundary.add_DirichletConstTemperature(table=np.ones((3, 2), dtype=bool))
				simulation.simulate(material=mat, boundary=boundary, overwrite=True)

			else :
				simuOutput = decoder(simulation._filename)
				simuOutput.plot_results(extension='_TH.pdf', plotLegend=True)
				