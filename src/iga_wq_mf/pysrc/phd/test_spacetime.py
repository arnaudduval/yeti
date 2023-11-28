from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem

def conductivityProperty(temperature):
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
	return Kprop 

def capacityProperty(temperature):
	Cprop = (1 + np.exp(-np.abs(temperature)))
	return Cprop

def powerDensity(args:dict):
	position = args['Position']; timespan = args['Time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		f[:, i] = 10*((x - 0.5)**2 + (y - 0.5)**2)*timespan[i]
	return np.ravel(f, order='F')

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 2, 5

# Create model 
geoArgs = {'name': 'qa', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
quadArgs = {'quadrule': 'iga', 'type': 'leg'}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)

# Create time span
nbel = 32
crv = createUniformCurve(degree, nbel, 0.5)
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
Fext = problem.compute_volForce(powerDensity, 
								{'Position':problem.part.qpPhy, 
								'Time':problem.time.qpPhy})

u_guess = np.zeros(np.prod(stnbctrlpts)); u_guess[boundary.thdod] = 0.0
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
for pcgmethod in ['C', 'JMC', 'TDC']:
	problem._methodPCG = pcgmethod
	u_sol, resPCG = problem.solveFourierSTHeatProblem(u_guess, Fext)
	if pcgmethod   == 'WP' : pcgname = 'w.o. preconditioner'
	elif pcgmethod == 'C'  : pcgname = 'Classic FD method'
	elif pcgmethod == 'TDC': pcgname = 'Literature'
	elif pcgmethod == 'JMC': pcgname = 'This work'
		
	for i in range(len(resPCG)):
		ax.semilogy(resPCG[i], marker=MARKERLIST[i], label=pcgname)
	ax.set_xlim(right=100, left=0)
	ax.set_ylim(top=10.0, bottom=1e-12)
	ax.set_xlabel('Number of iterations of BiCGSTAB solver')
	ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_2}{||b||_2}$')
	ax.legend()
	fig.tight_layout()
	fig.savefig(folder+'PCGresidue2'+'.pdf')

# u_sol = np.reshape(u_sol, (problem.part.nbctrlpts_total, problem.time.nbctrlpts), order='F')
# modelPhy.exportResultsCP(fields={'Ulast': u_sol[:, -1], 'Ustart': u_sol[:, 0]}, folder=folder)