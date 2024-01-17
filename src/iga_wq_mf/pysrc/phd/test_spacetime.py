from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem
import pickle

NLCASE = 3 # 0, 1 or 2
c = 0.005

def conductivityProperty(args, nlcase=NLCASE):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			if nlcase==0: Kprop[i, j, :] = Kref[i, j]
			if nlcase==1: Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
			if nlcase>=2: Kprop[i, j, :] = Kref[i, j]*(1+0.5*np.exp(-0.25*np.abs(temperature))*(np.sin(temperature)))
	return Kprop 

def conductivityDersProperty(args, nlcase=NLCASE):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			if nlcase==0: Kprop[i, j, :] = np.zeros(len(temperature))
			if nlcase==1: Kprop[i, j, :] = -Kref[i, j]*2.0*np.sign(temperature)*np.exp(-np.abs(temperature))
			if nlcase>=2: Kprop[i, j, :] = Kref[i, j]*np.exp(-0.25*np.abs(temperature))*(0.5*np.cos(temperature)
													-0.125*np.sign(temperature)*np.sin(temperature))
	return Kprop 

def capacityProperty(args, nlcase=NLCASE):
	temperature = args['temperature']
	if nlcase==0: Cprop = np.ones(len(temperature))
	if nlcase==1: Cprop = (1.0 + np.exp(-np.abs(temperature)))
	if nlcase>=2: Cprop = 1 + np.exp(-0.1*temperature**2)+0.25*np.sin(10*temperature)
	return Cprop

def capacityDersProperty(args, nlcase=NLCASE):
	temperature = args['temperature']
	if nlcase==0: Cprop = np.zeros(len(temperature))
	if nlcase==1: Cprop = -np.sign(temperature)*np.exp(-np.abs(temperature))
	if nlcase>=2: Cprop = 2.5*np.cos(10*temperature)-0.2*np.exp(-0.1*temperature**2)*temperature
	return Cprop

def exactTemperature(qpPhy, nlcase=NLCASE):
	x = qpPhy[0, :]; y = qpPhy[1, :]; t = qpPhy[2, :]
	if nlcase<=2: u = c*(-5*x + 6*y + 45)*(5*x + 6*y - 45)*x*(x-6)*np.sin(np.pi*t)
	if nlcase>2:  u = c*(-6*x + y + 10)*(6*x + y - 10)*x*(2*x-3)*np.sin(np.pi*t)
	return u

def powerDensity(args:dict, nlcase=NLCASE):
	position = args['Position']; timespan = args['Time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); f = np.zeros((nc_sp, nc_tm))
	if nlcase<=2:
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
	if nlcase>2: 
		for i in range(nc_tm):
			t = timespan[i]
			u = c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10)*(6*x + y - 10)
			f[:, i] = (
				((np.exp(-np.abs(u)/4)*np.sin(u))/2 + 1)*(
					72*c*x*np.sin(np.pi*t)*(2*x - 3) 
					- 12*c*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					+ 12*c*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10) 
					- 24*c*x*np.sin(np.pi*t)*(y - 6*x + 10) 
					+ 24*c*x*np.sin(np.pi*t)*(6*x + y - 10) 
					- 4*c*np.sin(np.pi*t)*(y - 6*x + 10)*(6*x + y - 10)
					) 
				#
				- 2*((np.exp(-np.abs(u)/4)*np.sin(u))/4 + 1/2)*(
					c*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					+ c*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10) 
					+ 2*c*x*np.sin(np.pi*t)*(y - 6*x + 10) 
					+ 2*c*x*np.sin(np.pi*t)*(6*x + y - 10)
					) 
				#
				- (np.exp(-np.abs(u)/4)*np.cos(u)*(c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					+ c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10)) 
					- (np.exp(-np.abs(u)/4)*np.sign(u)*np.sin(u)*(c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					+ c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10)))/4
					)*(
					c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					+ c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10)
					) 
				#
				- ((np.exp(-np.abs(u)/4)*np.cos(u)*(c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
						+ c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10))
						)/4 
					- (np.exp(-np.abs(u)/4)*np.sign(u)*np.sin(u)*(c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
						+ c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10))
						)/16
					)*(
					2*c*x*np.sin(np.pi*t)*(y - 6*x + 10)*(6*x + y - 10) 
					+ 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					- 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10) 
					+ c*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10)*(6*x + y - 10)
					) 
				#
				- ((np.exp(-np.abs(u)/4)*np.cos(u)*(2*c*x*np.sin(np.pi*t)*(y - 6*x + 10)*(6*x + y - 10) 
						+ 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
						- 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10) 
						+ c*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10)*(6*x + y - 10))
					)/4 
					- (np.exp(-np.abs(u)/4)*np.sign(u)*np.sin(u)*(2*c*x*np.sin(np.pi*t)*(y - 6*x + 10)*(6*x + y - 10) 
						+ 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
						- 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10) 
						+ c*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10)*(6*x + y - 10)))/16
					)*(
					c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					+ c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10)
					) 
				#
				- ((np.exp(-np.abs(u)/4)*np.cos(u)*(2*c*x*np.sin(np.pi*t)*(y - 6*x + 10)*(6*x + y - 10) 
					+ 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					- 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10) 
					+ c*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10)*(6*x + y - 10))
					)/2 
					- (np.exp(-np.abs(u)/4)*np.sign(u)*np.sin(u)*(2*c*x*np.sin(np.pi*t)*(y - 6*x + 10)*(6*x + y - 10) 
						+ 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
						- 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10) 
						+ c*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10)*(6*x + y - 10)))/8
					)*(
					2*c*x*np.sin(np.pi*t)*(y - 6*x + 10)*(6*x + y - 10) 
					+ 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10) 
					- 6*c*x*np.sin(np.pi*t)*(2*x - 3)*(6*x + y - 10) 
					+ c*np.sin(np.pi*t)*(2*x - 3)*(y - 6*x + 10)*(6*x + y - 10)
					) 
				#
				- 2*c*x*np.sin(np.pi*t)*(np.exp(-np.abs(u)/4)*np.sin(u) + 2)*(2*x - 3) 
				+ c*x*np.pi*np.cos(np.pi*t)*(2*x - 3)*(np.exp(-(c**2*x**2*np.sin(np.pi*t)**2*(2*x - 3)**2*(y - 6*x + 10)**2*(6*x + y - 10)**2)/10) 
				+ np.sin(10*u)/4 + 1)*(y - 6*x + 10)*(6*x + y - 10) 
			)
	return np.ravel(f, order='F')

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/'
if not os.path.isdir(folder): os.mkdir(folder)

def simulate(degree, cuts, quadArgs, problemArgs={}, uguess=None, isRandom=False, nlcase=NLCASE):
	# Create model 
	if nlcase<=2:
		geoArgs = {'name': 'tp', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
	if nlcase>2:
		geoArgs = {'name': 'tp', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int),
					'extra':{'XY':np.array([[0.0, -10], [1.5, -1], [1.5, 1], [0.0, 10]])}}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	# modelPhy = part(modelIGA, quadArgs={'quadrule': 'iga', 'type': 'leg'})

	# Create time span
	if nlcase<=2:
		crv = createUniformCurve(degree, 2**cuts, 1.)
	if nlcase>2:
		crv = createUniformCurve(degree, 2**cuts, 2.)

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
	if uguess is None and isRandom: uguess = np.random.uniform(-2, 5, np.prod(stnbctrlpts))
	if uguess is None and not isRandom: uguess = np.zeros(np.prod(stnbctrlpts))
	uguess[boundary.thdod] = 0.0
	if len(problemArgs) == 0: problemArgs={'isfull':True, 'isadaptive':True}
	isfull = problemArgs.get('isfull', True)
	isadaptive = problemArgs.get('isadaptive', True)
	problem._nIterNewton = problemArgs.get('NewtonIter', 10)
	problem._Krylov = problemArgs.get('Krylov', 'BICG')
	problem._KrylovPreconditioner = problemArgs.get('KrylovPrecond', 'JMC')
	problem._thresholdKrylov = problemArgs.get('KrylovThreshold', 1e-10)
	u_sol, resPCG, resNewton = problem.solveFourierSTHeatProblem(uguess, Fext, isfull=isfull, isadaptive=isadaptive)
	return problem, u_sol, resPCG, resNewton

# normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
# onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
# onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}

# degree_list = np.array([1, 2, 3, 4])
# cuts_list   = np.arange(1, 7)
# fig, ax = plt.subplots(figsize=(8, 6))
# # for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 1, 2], [normalPlot, onlyMarker1, onlyMarker2]):
# for quadrule, quadtype, plotpars in zip(['iga'], ['leg'], [normalPlot]):
# 	quadArgs = {'quadrule': quadrule, 'type': quadtype}
# 	error_list = np.ones(len(cuts_list))

# 	for i, degree in enumerate(degree_list):
# 		color = COLORLIST[i]
# 		for j, cuts in enumerate(cuts_list):
# 			nbels = 2**cuts_list
# 			problem, displacement, _, _ = simulate(degree, cuts, quadArgs)
# 			error_list[j] = problem.normOfError(displacement, normArgs={'type':'L2', 
# 															'exactFunction':exactTemperature}, 
# 															isRelative=False)
			
# 		if quadrule == 'iga': 
# 			ax.loglog(nbels, error_list, label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
# 						markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			
# 			slope = np.polyfit(np.log10(nbels[2:]),np.log10(error_list[2:]), 1)[0]
# 			slope = round(slope, 1)
# 			annotation.slope_marker((nbels[-2], error_list[-2]), slope, 
# 							poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)			
# 		else: 
# 			ax.loglog(nbels, error_list, color=color, marker=plotpars['marker'], markerfacecolor='w',
# 					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			
# 		fig.savefig(folder + 'SPTNonLinearConvergenceL2'+str(NLCASE)+'.pdf')

# # ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
# # 				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ 2")
# # ax.loglog([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
# # 		markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label="IGA-WQ 4")

# ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L_2(\Pi)}$')
# ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
# # ax.set_ylim(top=1e2, bottom=1e-9)
# ax.set_ylim(top=1e0, bottom=1e-11)
# ax.legend(loc='lower left')
# fig.tight_layout()
# fig.savefig(folder + 'SPTNonLinearConvergenceL2'+str(NLCASE)+'.pdf')

# # ===========================================
fig1, ax1 = plt.subplots(figsize=(8, 6))
fig2, ax2 = plt.subplots(figsize=(8, 6))
fig3, ax3 = plt.subplots(figsize=(8, 6))
figs = [fig1, fig2, fig3]
axs  = [ax1, ax2, ax3]

degree, cuts = 4, 4
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
Niterlist = range(1, 9)
legendname = ['Exact regular', 'Exact modified', 'Inexact regular', 'Inexact modified']
linestyle_list = ['-', '--', '-', '--']
marker_list = ['s', 's', 'o', 'o']

for i, isadaptive in enumerate([False, True]):
	for j, isfull in enumerate([True, False]):
		l = j + i*2
		krylov_list = np.ones(len(Niterlist))
		newton_list = np.ones(len(Niterlist))
		error_list  = np.ones(len(Niterlist))
		for k, niter in enumerate(Niterlist):
			problemArgs = {'isfull':isfull, 'isadaptive':isadaptive, 'NewtonIter':niter}
			problem, displacement, resPCG, resNewton = simulate(degree, cuts, quadArgs, problemArgs=problemArgs)
			newton_list[k] = resNewton[-1]
			error_list[k]  = problem.normOfError(displacement, normArgs={'type':'L2', 
															'exactFunction':exactTemperature}, 
															isRelative=False)
			resPCGclean = np.array([])
			for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])
			krylov_list[k] = len(resPCGclean)
		
		for caseplot, fig, ax in zip(range(1, 4), figs, axs):
			if caseplot == 1:
				yy = error_list; xx = np.append([0], krylov_list[:-1])
				# ylim = [20, 1e-3]; xlim = [150, 0]
				# ylim = [20, 1e-6]; xlim = [200, 0]
				ylim = [None, None]; xlim = [250, 0]
				ylabel = r'$\displaystyle ||u - u^h||_{L_2(\Pi)}$'
				xlabel = 'Total number of iterations of ' + problem._Krylov + ' solver'
			elif caseplot == 2:
				yy = newton_list; xx = np.append([0], krylov_list[:-1])
				# ylim = [1e2, 1e-6]; xlim = [150, 0]
				# ylim = [1e1, 1e-7]; xlim = [200, 0]
				ylim = [1e1, 1e-7]; xlim = [250, 0]
				ylabel = 'Norm of Newton residue'
				xlabel = 'Total number of iterations of ' + problem._Krylov + ' solver'
			elif caseplot == 3:
				yy = resNewton; xx = np.arange(0, len(resNewton))
				# ylim = [1e2, 1e-6]; xlim = [8, 0]
				# ylim = [1e1, 1e-7]; xlim = [8, 0]
				ylim = [1e1, 1e-7]; xlim = [8, 0]
				ylabel = 'Norm of Newton residue'
				xlabel = 'Number of Newton iterations'

			ax.semilogy(xx, yy, label=legendname[l], marker=marker_list[l], linestyle=linestyle_list[l])
			# ax.set_xlim(right=xlim[0], left=xlim[1])
			# ax.set_ylim(top=ylim[0], bottom=ylim[1])
			ax.set_xlabel(xlabel)
			ax.set_ylabel(ylabel)
			ax.legend()
			fig.tight_layout()
			fig.savefig(folder+'NLConvergence_iters'+str(degree)+str(caseplot)+str(NLCASE)+'.pdf')

#################################################################
## OLD TEST
##################################################################
# degree, cuts = 4, 4
# quadArgs = {'quadrule': 'iga', 'type': 'leg'}
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# for j, pcgmethod in enumerate(['WP', 'C', 'JMC', 'TDC']):
# 	problemArgs = {'isfull':False, 'isadaptive':False, 'NewtonIter':1, 
# 				'Krylov':'BICG', 'KrylovPrecond':pcgmethod, 'KrylovThreshold':1e-12}
# 	problem, displacement, resPCG, resNewton = simulate(degree, cuts, quadArgs, problemArgs=problemArgs)
# 	if pcgmethod == 'WP': pcgname = 'w.o. preconditioner' 
# 	elif pcgmethod == 'C' : pcgname = 'Classic FD method'
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

## ##################################################################
# degree, cuts = 4, 4
# quadArgs = {'quadrule': 'iga', 'type': 'leg'}
# stnbctrlpts = [int(2**cuts+degree) for i in range(3)]
# uguess =  np.random.uniform(-2, 5, np.prod(stnbctrlpts))
# for krylovmet in ['BICG', 'GMRES']:
# 	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7, 4))
# 	for j, isfull in enumerate([True, False]):
# 		problemArgs = {'isfull':isfull, 'isadaptive':False, 'NewtonIter':1, 'Krylov':krylovmet, 'KrylovThreshold':1e-12}
# 		problem, displacement, resPCG, resNewton = simulate(degree, cuts, quadArgs, problemArgs=problemArgs, uguess=uguess, isRandom=True)

# 		# Create continous resPCG
# 		if isfull: name=r'$\mathsf{P}$' + ' for ' +  r'$\mathsf{A}+\mathsf{B}$'
# 		else: name=r'$\mathsf{P}$' + ' for ' +  r'$\mathsf{A}$'
# 		resPCGclean = np.array([])
# 		for pcglist in resPCG: resPCGclean = np.append(resPCGclean, pcglist[np.nonzero(pcglist)])
# 		ax.semilogy(resPCGclean, marker='x', linestyle=':', label=name)

# 	ax.set_xlim(right=50, left=0)
# 	ax.set_ylim(top=10.0, bottom=1e-12)
# 	ax.set_xlabel('Number of iterations of ' + problem._Krylov + ' solver')
# 	ax.set_ylabel('Relative residue')
# 	ax.legend()
# 	fig.tight_layout()
# 	fig.savefig(folder+problem._Krylov+'residueNL2'+'.pdf')
