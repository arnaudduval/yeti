from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import heatproblem

def conductivityProperty(args:dict):
	P = args.get('position')
	Kref  = np.array([[1, 0.5],[0.5, 2]])
	Kprop = np.zeros((2, 2, np.size(P, axis=1)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j] 
	return Kprop 

def powerDensity_quartCircle(args:dict):
	""" u = sin(pi*x)*sin(pi*y)*sin(pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	P = args.get('position')
	x = P[0, :]
	y = P[1, :]

	f = (3*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- 16*y**2*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 8*x*y*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 4) 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 4) 
	- 8*x**2*np.sin(np.pi*x)*np.sin(np.pi*y)
	)

	return f

def exactTemperature_quartCircle(P: list):
	""" u = sin(pi*x)*sin(pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]

	u = np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 -1)*(x**2 + y**2 - 4)
	return u

def exactDiffTemperature_quartCircle(P: list):
	""" u = sin(pi*x)*sin(pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
		uders = grad(u)
	"""
	x = P[0, :]
	y = P[1, :]

	uders = np.zeros((1, 2, np.size(P, axis=1)))

	uders[0, 0, :] = (2*x*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
				+ 2*x*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
				+ np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
				)

	uders[0, 1, :] = (2*y*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
				+ 2*y*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4)
				+ np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
				)

	return uders

# Set global variables
degree, cuts = 2, 4
geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'Rin':1.0, 'Rex':2.0}
			}

material = heatmat()
material.addConductivity(conductivityProperty, isIsotropic=False)				

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs={'quadrule': 'iga', 'type': 'leg'})

# Set Dirichlet boundaries
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))

# Solve elastic problem
start=time.time()
problem = heatproblem(material, modelPhy, boundary)
Fext = problem.compute_volForce(powerDensity_quartCircle)
temperature = np.zeros(problem.part.nbctrlpts_total)
problem.solveFourierSteadyProblem(temperature, Fext=Fext)
