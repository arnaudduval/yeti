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
			Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
			# Kprop[i, j, :] = Kref[i, j]
	return Kprop 

def capacityProperty(args):
	temperature = args['temperature']
	Cprop = (1.0 + np.exp(-np.abs(temperature)))
	# Cprop = np.ones(shape=np.shape(temperature))
	return Cprop

def conductivityDersProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = -Kref[i, j]*2.0*np.sign(temperature)*np.exp(-np.abs(temperature))
	return Kprop 

def capacityDersProperty(args):
	temperature = args['temperature']
	Cprop = -np.sign(temperature)*np.exp(-np.abs(temperature))
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

def simulate(degree, cuts, quadArgs, problemArgs={}):
	# Create model 
	geoArgs = {'name': 'tp', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

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
	dirichlet_table[0, 0] = 0; dirichlet_table[0, 1] = 0
	stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timespan.nbctrlpts])
	boundary = boundaryCondition(stnbctrlpts)
	boundary.add_DirichletConstTemperature(table=dirichlet_table)

	problem = stheatproblem(material, modelPhy, timespan, boundary)
	# External heat force
	Fext = problem.compute_volForce(powerDensity, 
									{'Position':problem.part.qpPhy, 
									'Time':problem.time.qpPhy})
	u_guess = np.zeros(np.prod(stnbctrlpts)); u_guess[boundary.thdod] = 0.0
	if len(problemArgs) == 0: problemArgs={'isfull':True, 'isadaptive':True, 'NewtonIter':10}
	problem._nIterNewton = problemArgs['NewtonIter']
	u_sol, resPCG, resNewton = problem.solveFourierSTHeatProblem(u_guess, Fext, isfull=problemArgs['isfull'], 
															isadaptive=problemArgs['isadaptive'])
	return problem, u_sol, resPCG, resNewton

normalPlot  = {'marker': 'o', 'linestyle': '-', 'markersize': 10}
onlyMarker1 = {'marker': '.', 'linestyle': ':', 'markersize': 6}

with open(folder + 'refpart.pkl', 'rb') as inp:
	part_ref = pickle.load(inp)
with open(folder + 'reftime.pkl', 'rb') as inp:
	time_ref = pickle.load(inp)
u_ref = np.load(folder + 'refu.npy')

# # ===========================================
# degree_list = np.array([1, 2, 3, 4])
# cuts_list   = np.arange(1, 6)
# fig, ax = plt.subplots(figsize=(8, 6))
# for quadrule, quadtype, plotpars in zip(['iga', 'wq'], ['leg', 1], [normalPlot, onlyMarker1]):
# 	quadArgs = {'quadrule': quadrule, 'type': quadtype}
# 	error_list = np.ones(len(cuts_list))

# 	for i, degree in enumerate(degree_list):
# 		color = COLORLIST[i]
# 		for j, cuts in enumerate(cuts_list):
# 			problem, displacement, _, _ = simulate(degree, cuts, quadArgs)
# 			error_list[j] = problem.normOfError(displacement, normArgs={'type':'L2', 
# 															'part_ref':part_ref, 'time_ref':time_ref, 'u_ref':u_ref}, 
# 															isRelative=False)
			
# 		if quadrule == 'iga': 
# 			ax.loglog(2**cuts_list, error_list, label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
# 						markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			
# 		else: 
# 			ax.loglog(2**cuts_list, error_list, color=color, marker=plotpars['marker'], markerfacecolor='w',
# 					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
		
# 		ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L_2(\Pi)}$')
# 		ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
# 		ax.set_ylim(top=1e1, bottom=1e-3)
# 		ax.legend()
# 		fig.tight_layout()
# 		fig.savefig(folder + 'NLConvergence_meshsize' + '.pdf')

# # ===========================================
fig, ax = plt.subplots(figsize=(8, 6))
degree, cuts = 4, 3
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
Niterlist = range(1, 9)
legendname = ['Inexact regular', 'Inexact modified', 'Exact regular', 'Exact modified']
caseplot = 1

for i, isadaptive in enumerate([True, False]):
	for j, isfull in enumerate([True, False]):
		l = j + i*2
		krylov_list = np.ones(len(Niterlist))
		newton_list = np.ones(len(Niterlist))
		error_list  = np.ones(len(Niterlist))
		for k, niter in enumerate(Niterlist):
			problemArgs = {'isfull':isfull, 'isadaptive':isadaptive, 'NewtonIter':niter}
			problem, displacement, resPCG, resNewton = simulate(degree, cuts, quadArgs, problemArgs=problemArgs)
			newton_list[k] = resNewton[-1]
			error_list[k] = problem.normOfError(displacement, normArgs={'type':'L2', 
												'part_ref':part_ref, 'time_ref':time_ref, 'u_ref':u_ref}, 
												isRelative=False)
			resPCGclean = np.array([])
			for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])
			krylov_list[k] = len(resPCGclean)

		if caseplot == 1:
			yy = error_list; xx = np.append([0], krylov_list[:-1])
			ylim = [20, 1e-2]; xlim = [200, 0]
			ylabel = r'$\displaystyle ||u - u^h||_{L_2(\Pi)}$'
			xlabel = 'Total number of iterations of ' + problem._Krylov + ' solver'
		elif caseplot == 2:
			yy = newton_list; xx = np.append([0], krylov_list[:-1])
			ylim = [1e2, 1e-4]; xlim = [200, 0]
			ylabel = 'Norm of Newton residue'
			xlabel = 'Total number of iterations of ' + problem._Krylov + ' solver'
		elif caseplot == 3:
			yy = resNewton; xx = np.arange(0, len(resNewton))
			ylim = [1e2, 1e-4]; xlim = [8, 0]
			ylabel = 'Norm of Newton residue'
			xlabel = 'Number of Newton iterations'

		ax.semilogy(xx, yy, marker=MARKERLIST[l], label=legendname[l])
		ax.set_xlim(right=xlim[0], left=xlim[1])
		ax.set_ylim(top=ylim[0], bottom=ylim[1])
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.legend()
		fig.tight_layout()
		fig.savefig(folder+'NLConvergence_iters'+str(degree)+str(caseplot)+'.pdf')

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
# u_guess = np.random.uniform(-1., 1., np.prod(stnbctrlpts)); u_guess[boundary.thdod] = 0.0
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# problem._Krylov = 'GMRES'
# problem._KrylovPreconditioner = 'TDC'
# problem._nIterNewton = 7

# for isfull in [True, False]:
# 	u_sol, resPCG = problem.solveFourierSTHeatProblem(u_guess, Fext, 
# 													isfull=isfull, 
# 													isadaptive=False)
# 	# Create continous resPCG
# 	if isfull: name='Inconsistent\npreconditioner'
# 	else: name='Consistent\npreconditioner'
# 	resPCGclean = np.array([])
# 	for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])
# 	ax.semilogy(resPCGclean, label=name)

# # ax.set_xlim(right=50, left=0)
# ax.set_ylim(top=10.0, bottom=1e-12)
# ax.set_xlabel('Number of iterations of ' + problem._Krylov + ' solver')
# ax.set_ylabel('Relative residue')
# ax.legend()
# fig.tight_layout()
# fig.savefig(folder+problem._Krylov+'residueNL'+'.pdf')
