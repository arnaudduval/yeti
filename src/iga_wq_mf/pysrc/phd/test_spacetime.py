from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem
import pickle

def conductivityProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			# Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
			Kprop[i, j, :] = Kref[i, j]*(1+0.5*np.exp(-0.25*np.abs(temperature))*(np.sin(temperature)))
			# Kprop[i, j, :] = Kref[i, j]
	return Kprop 

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

def capacityProperty(args):
	temperature = args['temperature']
	# Cprop = (1.0 + np.exp(-np.abs(temperature)))
	Cprop = 1 + np.exp(-0.1*temperature**2)+0.25*np.sin(10*temperature)
	# Cprop = np.ones(shape=np.shape(temperature))
	return Cprop

def capacityDersProperty(args):
	temperature = args['temperature']
	# Cprop = -np.sign(temperature)*np.exp(-np.abs(temperature))
	Cprop = 2.5*np.cos(10*temperature)-0.2*np.exp(-0.1*temperature**2)*temperature
	return Cprop

def exactTemperature(qpPhy):
	x = qpPhy[0, :]; y = qpPhy[1, :]; t = qpPhy[2, :]
	u = 0.001*(-5*x + 6*y + 45)*(5*x + 6*y - 45)*x*(x-6)*np.sin(np.pi*t)
	return u

def powerDensity(args:dict):
	position = args['Position']; timespan = args['Time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); f = np.zeros((nc_sp, nc_tm))
	c = 0.001
	for i in range(nc_tm):
		t = timespan[i]

		u = c*x*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)

		f[:, i] = (
			4*np.sign(u)*np.exp(-np.abs(u))*(
				6*c*x*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45) 
				+ 6*c*x*np.sin(np.pi*t)*(x - 6)*(5*x + 6*y - 45)
				)**2
			#
			- (2*np.exp(-np.abs(u)) + 1)*(
				10*c*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45) 
				- 10*c*np.sin(np.pi*t)*(x - 6)*(5*x + 6*y - 45) 
				+ 2*c*np.sin(np.pi*t)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
				- 50*c*x*np.sin(np.pi*t)*(x - 6) + 10*c*x*np.sin(np.pi*t)*(6*y - 5*x + 45) 
				- 10*c*x*np.sin(np.pi*t)*(5*x + 6*y - 45)
				) 
			#
			- 2*(np.exp(-np.abs(u)) + 1/2)*(
				6*c*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45) 
				+ 6*c*np.sin(np.pi*t)*(x - 6)*(5*x + 6*y - 45) 
				+ 6*c*x*np.sin(np.pi*t)*(6*y - 5*x + 45) 
				+ 6*c*x*np.sin(np.pi*t)*(5*x + 6*y - 45)
				) 
			#	
			+ 2*np.sign(u)*np.exp(-np.abs(u))*(
				c*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
				+ 5*c*x*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45) 
				- 5*c*x*np.sin(np.pi*t)*(x - 6)*(5*x + 6*y - 45) 
				+ c*x*np.sin(np.pi*t)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
				)**2 
			#
			+ 2*np.sign(u)*np.exp(-np.abs(u))*(
				6*c*x*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45) 
				+ 6*c*x*np.sin(np.pi*t)*(x - 6)*(5*x + 6*y - 45)
				)*(
				c*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
				+ 5*c*x*np.sin(np.pi*t)*(x - 6)*(6*y - 5*x + 45) 
				- 5*c*x*np.sin(np.pi*t)*(x - 6)*(5*x + 6*y - 45) 
				+ c*x*np.sin(np.pi*t)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
				) 
			#
			- 72*c*x*np.sin(np.pi*t)*(4*np.exp(-np.abs(u)) + 2)*(x - 6) 
			#
			+ c*x*np.pi*np.cos(np.pi*t)*(np.exp(-np.abs(u)) + 1)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
		)

	return np.ravel(f, order='F')

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

def simulate(degree, cuts, quadArgs, problemArgs={}, isRandom=False):
	# Create model 
	geoArgs = {'name': 'tp', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	# modelPhy = part(modelIGA, quadArgs=quadArgs)
	modelPhy = part(modelIGA, quadArgs={'quadrule': 'iga', 'type': 'leg'})

	# Create time span
	crv = createUniformCurve(degree, 2**cuts, 1.)
	timespan = part1D(crv, {'quadArgs': quadArgs})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 
	material.addConductivityDers(conductivityDersProperty, isIsotropic=False) 
	material.addCapacityDers(capacityDersProperty, isIsotropic=False) 

	# Block boundaries
	dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
	stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timespan.nbctrlpts])
	boundary = boundaryCondition(stnbctrlpts)
	boundary.add_DirichletConstTemperature(table=dirichlet_table)

	problem = stheatproblem(material, modelPhy, timespan, boundary)
	# External heat force
	Fext = problem.compute_volForce(powerDensity, 
									{'Position':problem.part.qpPhy, 
									'Time':problem.time.qpPhy})
	if isRandom: u_guess = np.random.uniform(-2, 5, np.prod(stnbctrlpts))
	else: u_guess = np.zeros(np.prod(stnbctrlpts))
	u_guess[boundary.thdod] = 0.0
	if len(problemArgs) == 0: problemArgs={'isfull':True, 'isadaptive':False, 'NewtonIter':10}
	problem._nIterNewton = problemArgs['NewtonIter']
	problem._Krylov = problemArgs.get('Krylov', 'BICG')
	u_sol, resPCG, resNewton = problem.solveFourierSTHeatProblem(u_guess, Fext, isfull=problemArgs['isfull'], 
															isadaptive=problemArgs['isadaptive'])
	return problem, u_sol, resPCG, resNewton

normalPlot  = {'marker': 'o', 'linestyle': '-', 'markersize': 10}
onlyMarker1 = {'marker': '.', 'linestyle': ':', 'markersize': 6}
onlyMarker2 = {'marker': 'x', 'linestyle': 'None', 'markersize': 6}

with open(folder + 'refpart.pkl', 'rb') as inp:
	part_ref = pickle.load(inp)
with open(folder + 'reftime.pkl', 'rb') as inp:
	time_ref = pickle.load(inp)
u_ref = np.load(folder + 'refu.npy')

# # ===========================================
# degree_list = np.array([1, 2, 3, 4])
# cuts_list   = np.arange(1, 6)
# fig, ax = plt.subplots(figsize=(8, 6))
# for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 1, 2], [normalPlot, onlyMarker1, onlyMarker2]):
# 	quadArgs = {'quadrule': quadrule, 'type': quadtype}
# 	error_list = np.ones(len(cuts_list))

# 	for i, degree in enumerate(degree_list):
# 		color = COLORLIST[i]
# 		for j, cuts in enumerate(cuts_list):
# 			problem, displacement, _, _ = simulate(degree, cuts, quadArgs)
# 			# error_list[j] = problem.normOfError(displacement, normArgs={'type':'L2', 
# 			# 												'part_ref':part_ref, 'time_ref':time_ref, 'u_ref':u_ref}, 
# 			# 												isRelative=False)
# 			error_list[j] = problem.normOfError(displacement, normArgs={'type':'L2', 
# 															'exactFunction':exactTemperature}, 
# 															isRelative=False)
			
# 		if quadrule == 'iga': 
# 			ax.loglog(2**cuts_list, error_list, label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
# 						markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			
# 		else: 
# 			ax.loglog(2**cuts_list, error_list, color=color, marker=plotpars['marker'], markerfacecolor='w',
# 					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
		
# 		ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L_2(\Pi)}$')
# 		ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
# 		ax.set_ylim(top=1e2, bottom=1e-9)
# 		ax.legend()
# 		fig.tight_layout()
# 		fig.savefig(folder + 'NLConvergence_meshsize_timeleg' + '.pdf')

# # ===========================================
# fig, ax = plt.subplots(figsize=(8, 6))
# degree, cuts = 4, 3
# quadArgs = {'quadrule': 'iga', 'type': 'leg'}
# Niterlist = range(1, 9)
# legendname = ['Inexact regular', 'Inexact modified', 'Exact regular', 'Exact modified']
# caseplot = 3

# for i, isadaptive in enumerate([True, False]):
# 	for j, isfull in enumerate([True, False]):
# 		l = j + i*2
# 		krylov_list = np.ones(len(Niterlist))
# 		newton_list = np.ones(len(Niterlist))
# 		error_list  = np.ones(len(Niterlist))
# 		for k, niter in enumerate(Niterlist):
# 			problemArgs = {'isfull':isfull, 'isadaptive':isadaptive, 'NewtonIter':niter}
# 			problem, displacement, resPCG, resNewton = simulate(degree, cuts, quadArgs, problemArgs=problemArgs)
# 			newton_list[k] = resNewton[-1]
# 			# error_list[k] = problem.normOfError(displacement, normArgs={'type':'L2', 
# 			# 									'part_ref':part_ref, 'time_ref':time_ref, 'u_ref':u_ref}, 
# 			# 									isRelative=False)
# 			error_list[k] = problem.normOfError(displacement, normArgs={'type':'L2', 
# 															'exactFunction':exactTemperature}, 
# 															isRelative=False)
# 			resPCGclean = np.array([])
# 			for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])
# 			krylov_list[k] = len(resPCGclean)

# 		if caseplot == 1:
# 			yy = error_list; xx = np.append([0], krylov_list[:-1])
# 			ylim = [20, 1e-4]; xlim = [200, 0]
# 			ylabel = r'$\displaystyle ||u - u^h||_{L_2(\Pi)}$'
# 			xlabel = 'Total number of iterations of ' + problem._Krylov + ' solver'
# 		elif caseplot == 2:
# 			yy = newton_list; xx = np.append([0], krylov_list[:-1])
# 			ylim = [1e2, 1e-4]; xlim = [200, 0]
# 			ylabel = 'Norm of Newton residue'
# 			xlabel = 'Total number of iterations of ' + problem._Krylov + ' solver'
# 		elif caseplot == 3:
# 			yy = resNewton; xx = np.arange(0, len(resNewton))
# 			ylim = [1e2, 1e-4]; xlim = [8, 0]
# 			ylabel = 'Norm of Newton residue'
# 			xlabel = 'Number of Newton iterations'

# 		ax.semilogy(xx, yy, marker=MARKERLIST[l], label=legendname[l])
# 		ax.set_xlim(right=xlim[0], left=xlim[1])
# 		ax.set_ylim(top=ylim[0], bottom=ylim[1])
# 		ax.set_xlabel(xlabel)
# 		ax.set_ylabel(ylabel)
# 		ax.legend()
# 		fig.tight_layout()
# 		fig.savefig(folder+'NLConvergence_iters'+str(degree)+str(caseplot)+'.pdf')

#################################################################
## OLD TEST
#################################################################
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# degree, cuts = 4, 3
# quadArgs = {'quadrule': 'iga', 'type': 'leg'}
# problemArgs = {'isfull':True, 'isadaptive':False, 'NewtonIter':10}
# problem, displacement, resPCG = simulate(degree, cuts, quadArgs, problemArgs=problemArgs)

# # Create continous resPCG
# resPCGclean = np.array([])
# for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])

# ax.semilogy(resPCGclean)
# ax.set_xlim(right=400, left=0)
# ax.set_ylim(top=10.0, bottom=1e-12)
# ax.set_xlabel('Number of iterations of ' + problem._Krylov + ' solver')
# ax.set_ylabel('Relative residue')
# fig.tight_layout()
# fig.savefig(folder+problem._Krylov+'NL_FS3'+'.pdf')

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
degree, cuts = 4, 3
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))

for isfull in [True, False]:
	problemArgs = {'isfull':isfull, 'isadaptive':False, 'NewtonIter':1, 'Krylov':'BICG'}
	problem, displacement, resPCG, resNewton = simulate(degree, cuts, quadArgs, problemArgs=problemArgs, isRandom=True)

	# Create continous resPCG
	if isfull: name=r'$\mathsf{P}$' + ' for ' +  r'$\mathsf{A}+\mathsf{B}$'
	else: name=r'$\mathsf{P}$' + ' for ' +  r'$\mathsf{A}$'
	resPCGclean = np.array([])
	for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])
	ax.semilogy(resPCGclean, label=name)

ax.set_xlim(right=50, left=0)
ax.set_ylim(top=10.0, bottom=1e-12)
ax.set_xlabel('Number of iterations of ' + problem._Krylov + ' solver')
ax.set_ylabel('Relative residue')
ax.legend()
fig.tight_layout()
fig.savefig(folder+problem._Krylov+'residueNL2'+'.pdf')
