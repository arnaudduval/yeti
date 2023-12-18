from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem

def conductivityProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j]*(1+0.5*np.exp(-0.25*np.abs(temperature))*(np.sin(temperature)))
			# Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
			# Kprop[i, j, :] = Kref[i, j]
	return Kprop 

def capacityProperty(args):
	temperature = args['temperature']
	Cprop = 1 + np.exp(-0.1*temperature**2)+0.25*np.sin(10*temperature)
	# Cprop = (1.0 + np.exp(-np.abs(temperature)))
	# Cprop = np.ones(shape=np.shape(temperature))
	return Cprop

def conductivityDersProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j]*np.exp(-0.25*np.abs(temperature))*(0.5*np.cos(temperature)
													-0.125*np.sign(temperature)*np.sin(temperature))
			# Kprop[i, j, :] = -Kref[i, j]*2.0*np.sign(temperature)*np.exp(-np.abs(temperature))
	return Kprop 

def capacityDersProperty(args):
	temperature = args['temperature']
	Cprop = 2.5*np.cos(10*temperature)-0.2*np.exp(-0.1*temperature**2)*temperature
	# Cprop = -np.sign(temperature)*np.exp(-np.abs(temperature))
	return Cprop

def powerDensity(args:dict):
	position = args['Position']; timespan = args['Time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		z = timespan[i]
		f[:, i] = 0.001*(4*x*np.sin(np.pi*z)*(5*x + 6*y - 45) 
				- 94*x*np.sin(np.pi*z)*(x - 6) 
				- 16*x*np.sin(np.pi*z)*(6*y - 5*x + 45) 
				- 2*np.sin(np.pi*z)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
				- 16*np.sin(np.pi*z)*(x - 6)*(6*y - 5*x + 45) 
				+ 4*np.sin(np.pi*z)*(x - 6)*(5*x + 6*y - 45) 
				+ x*np.pi*np.cos(np.pi*z)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
				)
	return np.ravel(f, order='F')

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 3, 4

# Create model 
geoArgs = {'name': 'tp', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
quadArgs = {'quadrule': 'iga', 'type': 'leg'}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)

# Create time span
nbel = int(2**cuts)
crv = createUniformCurve(degree, nbel, 1.)
timespan = part1D(crv, {'quadArgs':{'quadrule': 'iga', 'type': 'leg'}})

# Add material 
material = heatmat()
material.addConductivity(conductivityProperty, isIsotropic=False) 
material.addCapacity(capacityProperty, isIsotropic=False) 
material.addConductivityDers(conductivityDersProperty, isIsotropic=False) 
material.addCapacityDers(capacityDersProperty, isIsotropic=False) 

# Block boundaries
dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
dirichlet_table[0, 1] = 0; 	dirichlet_table[0, 0] = 0
stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timespan.nbctrlpts])
boundary = boundaryCondition(stnbctrlpts)
boundary.add_DirichletConstTemperature(table=dirichlet_table)

# ---------------------
# Transient model
# ---------------------
problem = stheatproblem(material, modelPhy, timespan, boundary)

# External heat force
Fext = problem.compute_volForce(powerDensity, 
								{'Position':problem.part.qpPhy, 
								'Time':problem.time.qpPhy})
u_guess = np.zeros(np.prod(stnbctrlpts)); u_guess[boundary.thdod] = 0.0

##################################################################
# problem._Krylov = 'GMRES'
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# for j, pcgmethod in enumerate(['C', 'JMC', 'TDC']):
# 	problem._KrylovPreconditioner = pcgmethod
# 	u_sol, resPCG = problem.solveFourierSTHeatProblem(u_guess, Fext, isfull=False, isadaptive=False)
# 	if pcgmethod == 'C'    : pcgname = 'Classic FD method'
# 	elif pcgmethod == 'TDC': pcgname = 'Literature'
# 	elif pcgmethod == 'JMC': pcgname = 'This work'
# 	ax.semilogy(resPCG[0], marker=MARKERLIST[j], label=pcgname)

# ax.set_xlim(right=100, left=0)
# ax.set_ylim(top=10.0, bottom=1e-12)
# ax.set_xlabel('Number of iterations of ' + problem._Krylov + ' solver')
# ax.set_ylabel('Relative residue')
# ax.legend()
# fig.tight_layout()
# fig.savefig(folder+problem._Krylov+'residueLinear'+'.pdf')

##################################################################
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# problem._KrylovPreconditioner = 'JMC'
# u_sol, resPCG = problem.solveFourierSTHeatProblem(u_guess, Fext, 
# 												isfull=True, 
# 												isadaptive=False)
# print(u_sol.max(), u_sol.min())

# # Create continous resPCG
# resPCGclean = np.array([])
# for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])

# ax.semilogy(resPCGclean)
# ax.set_xlim(right=400, left=0)
# ax.set_ylim(top=10.0, bottom=1e-12)
# ax.set_xlabel('Number of iterations of ' + problem._Krylov + ' solver')
# ax.set_ylabel('Relative residue')
# fig.tight_layout()
# fig.savefig(folder+problem._Krylov+'NL_FS2'+'.pdf')

# u_sol = np.reshape(u_sol, (problem.part.nbctrlpts_total, problem.time.nbctrlpts), order='F')
# modelPhy.exportResultsCP(fields={'Ulast': u_sol[:, -1], 'Ustart': u_sol[:, 0]}, folder=folder)

##################################################################
# u_guess = np.random.uniform(-1., 1., np.prod(stnbctrlpts)); u_guess[boundary.thdod] = 0.0
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
problem._Krylov = 'GMRES'
problem._KrylovPreconditioner = 'TDC'
problem._nIterNewton = 7

for isfull in [True, False]:
	u_sol, resPCG = problem.solveFourierSTHeatProblem(u_guess, Fext, 
													isfull=isfull, 
													isadaptive=False)
	# Create continous resPCG
	if isfull: name='Inconsistent\npreconditioner'
	else: name='Consistent\npreconditioner'
	resPCGclean = np.array([])
	for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])
	ax.semilogy(resPCGclean, label=name)

# ax.set_xlim(right=50, left=0)
ax.set_ylim(top=10.0, bottom=1e-12)
ax.set_xlabel('Number of iterations of ' + problem._Krylov + ' solver')
ax.set_ylabel('Relative residue')
ax.legend()
fig.tight_layout()
fig.savefig(folder+problem._Krylov+'residueNL'+'.pdf')
