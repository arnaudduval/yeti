from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem

def conductivityProperty(temperature):
	cst = 10.0
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j]*cst*(1.0 + 2.0/(1.0 + np.exp(-5.0*(temperature-1.0))))
	return Kprop 

def capacityProperty(temperature):
	cst = 1.0
	Cprop = cst*(1 + np.exp(-2.0*abs(temperature)))
	return Cprop

def exactTemperature_quartCircle(args:dict):
	""" u = sin(pi*x)*sin(pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	position = args['Position']; timespan = args['Time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); u = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		u[:, i] = (np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 -1)*(x**2 + y**2 - 4))*timespan[i]
	return np.ravel(u, order='F')

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 3, 4

# Create model 
geoArgs = {'name': 'sq', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
quadArgs  = {'quadrule': 'iga', 'type': 'leg'}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)

# Create time span
nbel = 8
crv = createUniformCurve(degree, nbel, 1.0)
timespan = part1D(crv, {'quadArgs':{'quadrule': 'iga', 'type': 'leg'}})

# Add material 
material = heatmat()
material.addConductivity(conductivityProperty, isIsotropic=False) 
material.addCapacity(capacityProperty, isIsotropic=False) 

# Block boundaries
dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timespan.nbctrlpts])
boundary = boundaryCondition(stnbctrlpts)
boundary.add_DirichletConstTemperature(table=dirichlet_table)

# ---------------------
# Transient model
# ---------------------
problem = stheatproblem(material, modelPhy, timespan, boundary)

# External heat force
Fext = problem.compute_volForce(exactTemperature_quartCircle, 
								{'Position':problem.part.qpPhy, 'Time':problem.time.qpPhy})
temperature = np.zeros(np.prod(stnbctrlpts))
problem.solveFourierSTHeatProblem(temperature, Fext)