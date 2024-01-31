from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/paper/spacetime/'
if not os.path.isdir(folder): os.mkdir(folder)

extension = '.dat'
dataExist = True
FIG_CASE = 2
NONLIN_CASE = 1 # 0, 1, 2 or 3
if NONLIN_CASE==3: c = 0.01 # or 0.01
if NONLIN_CASE<3 : c = 0.05 # or 0.001
degree, cuts = 4, 4

def conductivityProperty(args, nlcase=NONLIN_CASE):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			if nlcase==0: Kprop[i, j, :] = Kref[i, j]
			if nlcase==1: Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
			if nlcase>=2: Kprop[i, j, :] = Kref[i, j]*(1+0.5*np.exp(-0.25*np.abs(temperature))*(np.sin(temperature)))
	return Kprop 

def conductivityDersProperty(args, nlcase=NONLIN_CASE):
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

def capacityProperty(args, nlcase=NONLIN_CASE):
	temperature = args['temperature']
	if nlcase==0: Cprop = np.ones(len(temperature))
	if nlcase==1: Cprop = (1.0 + np.exp(-np.abs(temperature)))
	if nlcase>=2: Cprop = 1 + np.exp(-0.1*temperature**2)+0.25*np.sin(10*temperature)
	return Cprop

def capacityDersProperty(args, nlcase=NONLIN_CASE):
	temperature = args['temperature']
	if nlcase==0: Cprop = np.zeros(len(temperature))
	if nlcase==1: Cprop = -np.sign(temperature)*np.exp(-np.abs(temperature))
	if nlcase>=2: Cprop = 2.5*np.cos(10*temperature)-0.2*np.exp(-0.1*temperature**2)*temperature
	return Cprop

def exactTemperature(qpPhy, nlcase=NONLIN_CASE):
	x = qpPhy[0, :]; y = qpPhy[1, :]; t = qpPhy[2, :]
	if nlcase<=2: u = c*(-5*x + 6*y + 45)*(5*x + 6*y - 45)*x*(x-6)*np.sin(np.pi*t)
	if nlcase>2:  u = c*(-6*x + y + 10)*(6*x + y - 10)*x*(2*x-3)*np.sin(np.pi*t)
	return u

def powerDensity(args:dict, nlcase=NONLIN_CASE):
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

def simulate(degree, cuts, quadArgs, uguess=None, problemArgs={}, nlcase=NONLIN_CASE):
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

	# Create time span
	if nlcase<=2: crv = createUniformCurve(degree, 2**cuts, 1.)
	if nlcase>2:  crv = createUniformCurve(degree, 2**cuts, 2.)
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
	
	if uguess is None: uguess = np.zeros(np.prod(stnbctrlpts))
	# if uguess is None: uguess = np.random.uniform(-2, 5, np.prod(stnbctrlpts))

	uguess[boundary.thdod] = 0.0
	problem._itersNL = 11
	isfull = problemArgs.get('isfull', True); isadaptive = problemArgs.get('isadaptive', True)
	output = problem.solveFourierSTHeatProblem(uguess, Fext, isfull=isfull, isadaptive=isadaptive)
	return problem, output

if not dataExist:

	if FIG_CASE == 1:
		degree_list = np.array([1, 2])
		cuts_list   = np.arange(6, 7)
		for quadrule, quadtype in zip(['iga'], ['leg']):
			sufix = '_' + quadrule + '_' + quadtype + '_' + str(NONLIN_CASE)
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			L2errorTable = np.zeros((len(degree_list)+1, len(cuts_list)+1))
			L2relerrorTable = np.zeros((len(degree_list)+1, len(cuts_list)+1))
			L2errorTable[0, 1:] = cuts_list; L2relerrorTable[0, 1:] = cuts_list
			L2errorTable[1:, 0] = degree_list; L2relerrorTable[1:, 0] = degree_list

			for i, degree in enumerate(degree_list):
				for j, cuts in enumerate(cuts_list):
					nbels = 2**cuts_list
					problem, output = simulate(degree, cuts, quadArgs)
					displacement = output['Solution'][-1]
					L2errorTable[i+1, j+1], L2relerrorTable[i+1, j+1] = problem.normOfError(displacement, 
																	normArgs={'type':'L2', 
																	'exactFunction':exactTemperature},)

					np.savetxt(folder+'L2error_meshpar_test'+sufix+extension, L2errorTable)
					np.savetxt(folder+'L2relerror_meshpar_test'+sufix+extension, L2relerrorTable)

	elif FIG_CASE == 2:

		quadArgs = {'quadrule': 'iga', 'type': 'leg'}
		meshpartext = str(NONLIN_CASE) + '_' + str(degree) + '_' + str(cuts) + '/'
		subfolderfolder = folder + meshpartext 
		if not os.path.isdir(subfolderfolder): os.mkdir(subfolderfolder)

		for [i, isadaptive], prefix1 in zip(enumerate([False, True]), ['exact', 'inexact']):
			for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
				prefix = prefix1 + '_' + prefix2 + '_'
				problemArgs = {'isfull':isfull, 'isadaptive':isadaptive}
				blockPrint()
				problem, output = simulate(degree, cuts, quadArgs, problemArgs=problemArgs)
				displacements = output['Solution']
				resKrylovs    = output['KrylovRes']
				resNewtons    = output['NewtonRes']
				L2error, L2relerror = [], []

				for displacement in displacements:
					err, relerr  = problem.normOfError(displacement, normArgs={'type':'L2', 
																	'exactFunction':exactTemperature})
					L2error.append(err); L2relerror.append(relerr)
				resKrylovclean = np.array([]); counter_list = [0]
				for _ in resKrylovs: 
					resKrylovclean = np.append(resKrylovclean, _[np.nonzero(_)])
					counter_list.append(counter_list[-1] + len(_[np.nonzero(_)]))
				enablePrint()
				print(len(counter_list), len(L2error), len(resNewtons))
				np.savetxt(subfolderfolder+prefix+'CumulKrylovRes'+extension, resKrylovclean)
				np.savetxt(subfolderfolder+prefix+'Inner_loops'+extension, counter_list)
				np.savetxt(subfolderfolder+prefix+'NewtonRes'+extension, resNewtons)
				np.savetxt(subfolderfolder+prefix+'L2error'+extension, L2error)
				np.savetxt(subfolderfolder+prefix+'L2relerror'+extension, L2relerror)

else:
	if FIG_CASE == 1:
		normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
		onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
		onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}
		plotoptions = [normalPlot, onlyMarker1, onlyMarker2]

		figname = folder + 'SPTNonLinearConvergenceL2'+str(NONLIN_CASE)+'.pdf'
		filenames = ['L2relerror_meshpar_iga_leg_']
		if FIG_CASE == 1:
			normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
			fig, ax = plt.subplots(figsize=(8, 6))

			for filename, plotops in zip(filenames, plotoptions):
				quadrule = filename.split('_')[2]
				table = np.loadtxt(folder+filename+str(NONLIN_CASE)+extension)	
				nbels   = 2**(table[0, 1:])
				degrees = table[1:, 0]
				errors  = table[1:, 1:]
				for i, degree in enumerate(degrees):
					color = COLORLIST[i]
					if quadrule == 'iga': 
						ax.loglog(nbels, errors[i, :], label='IGA-GL deg. '+str(int(degree)), color=color, marker=plotops['marker'],
									markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])
						
						slope = np.polyfit(np.log10(nbels[2:-1]),np.log10(errors[i, 2:-1]), 1)[0]
						slope = round(slope, 1)
						annotation.slope_marker((nbels[-2], errors[i, -2]), slope, 
										poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)			
					else: 
						ax.loglog(nbels, errors[i, :], color=color, marker=plotops['marker'], markerfacecolor='w',
								markersize=plotops['markersize'], linestyle=plotops['linestyle'])
							
					fig.savefig(figname)

			# ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
			# 				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ 2")
			# ax.loglog([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
			# 		markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label="IGA-WQ 4")

			ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$')
			ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
			ax.set_ylim(top=1e1, bottom=1e-11)
			ax.legend(loc='lower left')
			fig.tight_layout()
			fig.savefig(figname)

	elif FIG_CASE == 2:

		meshpartext = str(NONLIN_CASE) + '_' + str(degree) + '_' + str(cuts) + '/'
		subfolderfolder = folder + meshpartext 

		fig1, ax1 = plt.subplots(figsize=(8, 6))
		fig2, ax2 = plt.subplots(figsize=(8, 6))
		fig3, ax3 = plt.subplots(figsize=(8, 6))
		fig4, ax4 = plt.subplots(figsize=(8, 6))
		figs = [fig1, fig2, fig3, fig4]; axs  = [ax1, ax2, ax3, ax4]
		linestyle_list = ['-', '--', '-', '--']
		marker_list = ['o', 'o', 's', 's']

		for [i, isadaptive], prefix1 in zip(enumerate([True, False]), ['inexact', 'exact']):
			for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
				l = j + i*2
				legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
				prefix = prefix1 + '_' + prefix2 + '_'
				nbInnerLoops = np.loadtxt(subfolderfolder+prefix+'Inner_loops'+extension)
				newtonRes = np.loadtxt(subfolderfolder+prefix+'NewtonRes'+extension)
				L2relerror = np.loadtxt(subfolderfolder+prefix+'L2relerror'+extension)
				newtonRes = newtonRes/newtonRes[0]
				
				for caseplot, fig, ax in zip(range(1, 5), figs, axs):
					if caseplot == 1:
						yy = L2relerror; xx = nbInnerLoops[:len(L2relerror)]
						xlim = 10*np.ceil(np.max(nbInnerLoops)/10); ylim = [2, 0.2e-6]
						ylabel = r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$'
						xlabel = 'Number of inner iterations'
					elif caseplot == 2:
						yy = newtonRes; xx = nbInnerLoops[:len(newtonRes)]
						xlim = 10*np.ceil(np.max(nbInnerLoops)/10); ylim = [2, 5e-9]
						ylabel = 'Relative norm of outer residue'
						xlabel = 'Number of inner iterations'
					elif caseplot == 3:
						yy = newtonRes; xx = np.arange(0, len(newtonRes))
						xlim = 10; ylim = [2, 5e-9]
						ylabel = 'Relative norm of outer residue'
						xlabel = 'Number of outer iterations'
					elif caseplot == 4:
						yy = L2relerror; xx = np.arange(0, len(newtonRes))
						xlim = 10; ylim = [2, 0.2e-6]
						ylabel = r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$'
						xlabel = 'Number of outer iterations'

					ax.semilogy(xx, yy, label=legendname, marker=marker_list[l], linestyle=linestyle_list[l])
					ax.set_xlim(right=xlim, left=0)
					# ax.set_ylim(top=ylim[0], bottom=ylim[1])
					ax.set_xlabel(xlabel)
					ax.set_ylabel(ylabel)
					ax.legend()
					fig.tight_layout()
					fig.savefig(folder+'NLConvergence_iters'+str(NONLIN_CASE)+'_'+str(degree)+str(cuts)+str(caseplot)+'.pdf')

	elif FIG_CASE == 3:
		...
		# degree, cuts = 4, 4
		# quadArgs = {'quadrule': 'iga', 'type': 'leg'}
		# stnbctrlpts = [int(2**cuts+degree) for i in range(3)]
		# uguess =  np.random.random(np.prod(stnbctrlpts))
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
