from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import stheatproblem

def conductivityProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j]
	return Kprop 

def capacityProperty(args):
	temperature = args['temperature']
	Cprop = np.ones(shape=np.shape(temperature))
	return Cprop

def exactTemperature(qpPhy):
	x = qpPhy[0, :]; y = qpPhy[1, :]; t = qpPhy[2, :]
	u = (np.sin(np.pi*x)*np.sin(np.pi*y))*t
	return u

def powerDensity(args:dict):
	position = args['position']; timespan = args['time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = timespan[i]
		f[:, i] = (np.sin(np.pi*x)*np.sin(np.pi*y))*(1. + 2*np.pi**2*t)
	return np.ravel(f, order='F')

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2spacetime/'
if not os.path.isdir(folder): os.mkdir(folder)

def simulate(degree, cuts, quadArgs, problemArgs={}):
	# Create model 
	geoArgs = {'name': 'sq', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Create time span
	crv = createUniformOpenCurve(degree, 2**cuts, 1.0)
	timespan = part1D(crv, {'quadArgs':quadArgs})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
	stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timespan.nbctrlpts_total])
	boundary = boundaryCondition(stnbctrlpts)
	boundary.add_DirichletConstTemperature(table=dirichlet_table)

	problem = stheatproblem(material, modelPhy, timespan, boundary)
	# External heat force
	Fext = problem.compute_volForce(powerDensity, {'position':problem.part.qpPhy, 'time':problem.time.qpPhy})
	u_guess = np.zeros(np.prod(stnbctrlpts)); u_guess[boundary.thdod] = 0.0
	problem.addSolverConstraints(solverArgs=problemArgs)
	output = problem.solveFourierSTHeatProblem(u_guess, Fext, isfull=False, isadaptive=True)
	return problem, output

degree, cuts = 6, 5
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
start = time.time()
# for j, pcgmethod in enumerate(['C', 'JMC', 'TDC']):
for j, pcgmethod in enumerate(['JMC']):
	blockPrint()
	problemArgs = {'linearsolver':'GMRES', 'preconditioner':pcgmethod, 'thres_linear':1e-12}
	start = time.time()
	problem, output = simulate(degree, cuts, quadArgs, problemArgs=problemArgs)
	stop = time.time()
	resPCG = output['KrylovRes']
	enablePrint()
	print('time:%.2e'%(stop-start))

	if pcgmethod == 'WP': pcgname = 'w.o. preconditioner'
	elif pcgmethod == 'C' : pcgname = 'Classic FD method'
	elif pcgmethod == 'TDC': pcgname = 'Literature'
	elif pcgmethod == 'JMC': pcgname = 'This work'
	ax.semilogy(resPCG[0], marker=MARKERLIST[j], label=pcgname)
finish = time.time()
print('Time: %.4e' %(finish-start))