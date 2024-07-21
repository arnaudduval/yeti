from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import stheatproblem

def conductivityProperty(temperature):
	cst = 10.0
	Kref  = np.array([[1., 0.0, 0.0],[0.0, 1.0, 0.0], [0., 0., 1.0]])
	Kprop = np.zeros((3, 3, len(temperature)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j]*cst*(1.0 + 2.0/(1.0 + np.exp(-5.0*(temperature-1.0))))
	return Kprop 

def capacityProperty(temperature):
	cst = 1.0
	Cprop = cst*(1 + np.exp(-2.0*abs(temperature)))
	# Cprop = np.ones(shape=np.shape(temperature))
	return Cprop

def exactTemperature(args:dict):
	"""
		f = -div(lambda * grad(u))
	"""
	position = args['Position']; timespan = args['Time']
	x = position[0, :]; y = position[1, :]; z = position[2, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); u = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		u[:, i] = (np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z))*timespan[i]
	return np.ravel(u, order='F')

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 2, 4

# Create model 
geoArgs = {'name': 'cb', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
quadArgs  = {'quadrule': 'iga', 'type': 'leg'}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)

# Create time span
nbel = 4
crv = createUniformOpenCurve(degree, nbel, 1.0)
timespan = part1D(crv, {'quadArgs':{'quadrule': 'iga', 'type': 'leg'}})

# Add material 
material = heatmat()
material.addConductivity(conductivityProperty, isIsotropic=False) 
material.addCapacity(capacityProperty, isIsotropic=False) 

# Block boundaries
dirichlet_table = np.ones((4, 2)); dirichlet_table[-1, 1] = 0
stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timespan.nbctrlpts])
boundary = boundaryCondition(stnbctrlpts)
boundary.add_DirichletConstTemperature(table=dirichlet_table)

# Transient model
problem = stheatproblem(material, modelPhy, timespan, boundary)

# External heat force
Fext = problem.compute_volForce(exactTemperature, 
								{'Position':problem.part.qpPhy, 'Time':problem.time.qpPhy})
u_guess = np.zeros(np.prod(stnbctrlpts))
u_sol, resPCG = problem.solveFourierSTHeatProblem(u_guess, Fext)
u_sol = np.reshape(u_sol, (problem.part.nbctrlpts_total, problem.time.nbctrlpts), order='F')
modelPhy.postProcessingPrimal(fields={'Ulast': u_sol[:, -1], 'Ustart': u_sol[:, 0]}, folder=folder)