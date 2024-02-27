from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/spacetime/TR/'
if not os.path.isdir(folder): os.mkdir(folder)

# 1: 0.000001, 2: 0.0001, 3: 0.001
extension = '.dat'
FIG_CASE  = 1
DATAEXIST = True
ISLINEAR  = False
c = 0.01
degree, cuts = 2, 5

def conductivityProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			if ISLINEAR: Kprop[i, j, :] = Kref[i, j]
			else: Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
	return Kprop 

def conductivityDersProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			if ISLINEAR: Kprop[i, j, :] = 0
			else: Kprop[i, j, :] = -Kref[i, j]*2.0*np.sign(temperature)*np.exp(-np.abs(temperature))
	return Kprop 

def capacityProperty(args):
	temperature = args['temperature']
	if ISLINEAR: Cprop = np.ones(len(temperature))
	else: Cprop = (1.0 + np.exp(-np.abs(temperature)))
	return Cprop

def capacityDersProperty(args):
	temperature = args['temperature']
	if ISLINEAR: Cprop = np.zeros(len(temperature))
	else: Cprop = -np.sign(temperature)*np.exp(-np.abs(temperature))
	return Cprop

def exactTemperature(qpPhy):
	x = qpPhy[0, :]; y = qpPhy[1, :]; t = qpPhy[2, :]
	u = c*(-5*x + 6*y + 45)*(5*x + 6*y - 45)*np.sin(np.pi*x)*np.sin(np.pi*t)
	return u

def powerDensity(args:dict):
	position = args['Position']; timespan = args['Time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); f = np.zeros((nc_sp, nc_tm))
	if ISLINEAR:
		for i in range(nc_tm):
			t = timespan[i]
			f[:, i] = (4*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(5*x + 6*y - 45) 
						- 16*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45) 
						- 94*c*np.sin(np.pi*t)*np.sin(np.pi*x) 
						+ c*np.pi*np.cos(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
						+ c*np.pi**2*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
			)
	else: 
		for i in range(nc_tm):
			t = timespan[i]
			u = c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
			f[:, i] = (
				(2*np.exp(-np.abs(u)) + 1)*(
					50*c*np.sin(np.pi*t)*np.sin(np.pi*x) 
					- 10*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45) 
					+ 10*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(5*x + 6*y - 45) 
					+ c*np.pi**2*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45)*(5*x + 6*y - 45)) 
				- 2*(np.exp(-np.abs(u)) + 1/2)*(
					6*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45) 
					+ 6*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(5*x + 6*y - 45)) 
				+ 2*np.exp(-np.abs(u))*np.sign(u)*(
					5*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45) 
					- 5*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(5*x + 6*y - 45) 
					+ c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45)*(5*x + 6*y - 45))**2 
				+ 4*np.exp(-np.abs(u))*np.sign(u)*(
					6*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45) 
					+ 6*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(5*x + 6*y - 45))**2 
				+ 2*np.exp(-np.abs(u))*np.sign(u)*(
					6*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45) 
					+ 6*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(5*x + 6*y - 45))*(
						5*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45) 
						- 5*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(5*x + 6*y - 45) 
						+ c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45)*(5*x + 6*y - 45)) 
				- 72*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(4*np.exp(-np.abs(u)) + 2) 
				+ c*np.pi*np.cos(np.pi*t)*np.sin(np.pi*x)*(np.exp(-np.abs(u)) + 1)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
			)
	return np.ravel(f, order='F')

def simulate(degree, cuts, quadArgs, uguess=None, problemArgs={}):
	# Create model 
	geoArgs = {'name': 'TP', 'degree': degree*np.ones(3, dtype=int), 
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
	isfull = problemArgs.get('isfull', True); isadaptive = problemArgs.get('isadaptive', True)
	output = problem.solveFourierSTHeatProblem(uguess, Fext, isfull=isfull, isadaptive=isadaptive)
	return problem, output

if not DATAEXIST:

	if FIG_CASE == 1:
		lastsufix = 'linear' if ISLINEAR else 'nonlin'
		degree_list = np.array([1, 2, 3, 4])
		cuts_list   = np.arange(1, 6)
		for quadrule, quadtype in zip(['iga'], ['leg']):
			sufix = '_' + quadrule + '_' + quadtype + '_' + lastsufix
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			L2errorTable = np.zeros((len(degree_list)+1, len(cuts_list)+1))
			L2relerrorTable = np.zeros((len(degree_list)+1, len(cuts_list)+1))
			L2errorTable[0, 1:] = cuts_list; L2relerrorTable[0, 1:] = cuts_list
			L2errorTable[1:, 0] = degree_list; L2relerrorTable[1:, 0] = degree_list
			filename1 = folder+'L2error_meshpar'+sufix+extension
			filename2 = folder+'L2relerror_meshpar'+sufix+extension
			# if os.path.exists(filename1): raise Warning('File exist')
			# if os.path.exists(filename2): raise Warning('File exist')
			for j, cuts in enumerate(cuts_list):
				for i, degree in enumerate(degree_list):
					nbels = 2**cuts_list
					blockPrint()
					problem, output = simulate(degree, cuts, quadArgs)
					enablePrint()
					displacement = output['Solution'][-1]
					L2errorTable[i+1, j+1], L2relerrorTable[i+1, j+1] = problem.normOfError(displacement, 
																	normArgs={'type':'L2', 
																	'exactFunction':exactTemperature},)

					np.savetxt(filename1, L2errorTable)
					np.savetxt(filename2, L2relerrorTable)

	elif FIG_CASE == 2:

		quadArgs = {'quadrule': 'iga', 'type': 'leg'}
		meshpartext = '_' + str(degree) + '_' + str(cuts) + '/'
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
				thresholds    = output['Threshold']
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
				np.savetxt(subfolderfolder+prefix+'CumulKrylovRes'+extension, resKrylovclean)
				np.savetxt(subfolderfolder+prefix+'Inner_loops'+extension, counter_list)
				np.savetxt(subfolderfolder+prefix+'NewtonRes'+extension, resNewtons)
				np.savetxt(subfolderfolder+prefix+'L2error'+extension, L2error)
				np.savetxt(subfolderfolder+prefix+'L2relerror'+extension, L2relerror)
				np.savetxt(subfolderfolder+prefix+'threshold'+extension, np.array(thresholds))
else:
	if FIG_CASE == 1:
		normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
		onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
		onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}
		plotoptions = [normalPlot, onlyMarker1, onlyMarker2]

		lastsufix = 'linear' if ISLINEAR else 'nonlin'
		figname = folder + 'SPTNonLinearConvergenceL2'+lastsufix+'.pdf'
		filenames = ['L2relerror_meshpar_iga_leg_']

		normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
		fig, ax = plt.subplots(figsize=(8, 6))

		for filename, plotops in zip(filenames, plotoptions):
			quadrule = filename.split('_')[2]
			table = np.loadtxt(folder+filename+lastsufix+extension)	
			nbels   = 2**(table[0, 1:])
			degrees = table[1:, 0]
			errors  = table[1:, 1:]
			for i, degree in enumerate(degrees):
				color = COLORLIST[i]
				if quadrule == 'iga': 
					ax.loglog(nbels, errors[i, :], label='IGA-GL deg. '+str(int(degree)), color=color, marker=plotops['marker'],
								markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
					# slope = np.polyfit(np.log10(nbels[2:]),np.log10(errors[i, 2:]), 1)[0]
					# slope = round(slope, 1)
					# annotation.slope_marker((nbels[-1], errors[i, -1]), slope, 
					# 				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)			
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
		ax.set_xlim(left=1, right=100)
		ax.set_ylim(top=1e1, bottom=1e-7)
		ax.legend(loc='lower left')
		fig.tight_layout()
		fig.savefig(figname)

	elif FIG_CASE == 2:

		meshpartext = '_' + str(degree) + '_' + str(cuts) + '/'
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
						xlim = 10*np.ceil(np.max(nbInnerLoops)/10); ylim = np.power(10, np.floor(np.log10(np.min(L2relerror))))
						ylabel = r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$'
						xlabel = 'Number of inner iterations'
					elif caseplot == 2:
						yy = newtonRes; xx = nbInnerLoops[:len(newtonRes)]
						xlim = 10*np.ceil(np.max(nbInnerLoops)/10); ylim = np.power(10, np.floor(np.log10(np.min(newtonRes))))
						ylabel = 'Relative norm of outer residue'
						xlabel = 'Number of inner iterations'
					elif caseplot == 3:
						yy = newtonRes; xx = np.arange(0, len(newtonRes))
						xlim = 10; ylim = np.power(10, np.floor(np.log10(np.min(newtonRes))))
						ylabel = 'Relative norm of outer residue'
						xlabel = 'Number of outer iterations'
					elif caseplot == 4:
						yy = L2relerror; xx = np.arange(0, len(newtonRes))
						xlim = 10; ylim = np.power(10, np.floor(np.log10(np.min(L2relerror))))
						ylabel = r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$'
						xlabel = 'Number of outer iterations'

					ax.semilogy(xx, yy, label=legendname, marker=marker_list[l], linestyle=linestyle_list[l])
					ax.set_xlim(right=xlim, left=0)
					# ax.set_ylim(top=ylim[0], bottom=ylim[1])
					ax.set_ylim(top=1e1, bottom=ylim)
					ax.set_xlabel(xlabel)
					ax.set_ylabel(ylabel)
					ax.legend()
					fig.tight_layout()
					fig.savefig(folder+'NLConvergence_iters'+'_'+str(degree)+str(cuts)+str(caseplot)+'.pdf')

		fig, ax = plt.subplots(figsize=(8, 6))
		for [i, isadaptive], prefix1 in zip(enumerate([True]), ['inexact']):
			for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
				l = j + i*2
				legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
				prefix = prefix1 + '_' + prefix2 + '_'
				threshold = np.loadtxt(subfolderfolder+prefix+'threshold'+extension)

				yy = threshold; xx = np.arange(0, len(threshold))
				xlim = 10
				ylabel = 'Forcing term'
				xlabel = 'Number of outer iterations'

				ax.semilogy(xx, yy, label=legendname, marker=marker_list[l], linestyle=linestyle_list[l])
				ax.set_xlim(right=xlim, left=0)
				# ax.set_ylim(top=ylim[0], bottom=ylim[1])
				ax.set_ylim(top=1, bottom=1e-4)
				ax.set_xlabel(xlabel)
				ax.set_ylabel(ylabel)
				ax.legend()
				fig.tight_layout()
				fig.savefig(folder+'NLTolerance'+'_'+str(degree)+str(cuts)+'.pdf')

