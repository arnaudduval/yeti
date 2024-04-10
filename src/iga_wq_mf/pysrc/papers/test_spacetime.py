from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
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

TODOSIMU  = True
c = 0.0001
# degree, cuts = 3, 4
degree, cuts = 2, 5

def zeros(qpPhy):
	x = qpPhy[0, :]
	u = np.zeros(len(x))
	return u

def conductivityProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
	return Kprop 

def conductivityDersProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = -Kref[i, j]*2.0*np.sign(temperature)*np.exp(-np.abs(temperature))
	return Kprop 

def capacityProperty(args):
	temperature = args['temperature']
	Cprop = np.ones(len(temperature))
	return Cprop

def capacityDersProperty(args):
	temperature = args['temperature']
	Cprop = np.zeros(len(temperature))
	return Cprop

def exactTemperature(qpPhy):
	x = qpPhy[0, :]; y = qpPhy[1, :]; t = qpPhy[2, :]
	u = c*(-5*x + 6*y + 45)*(5*x + 6*y - 45)*x*(x-6)*np.sin(np.pi/6*x)*np.sin(np.pi*t)
	return u

def powerDensity(args:dict):
	from numpy import exp, sin, cos, pi, abs, sign
	position = args['position']; timespan = args['time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = timespan[i]
		u = c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
		f[:, i] = (
			4*exp(-abs(u))*sign(u)*(
				6*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45) 
				+ 6*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(5*x + 6*y - 45))**2 
			- 2*(exp(-abs(u)) + 1/2)*(
				6*c*x*sin(pi*t)*sin((pi*x)/6)*(6*y - 5*x + 45) 
				+ 6*c*x*sin(pi*t)*sin((pi*x)/6)*(5*x + 6*y - 45) 
				+ 6*c*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45) 
				+ 6*c*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(5*x + 6*y - 45) 
				+ c*x*pi*cos((pi*x)/6)*sin(pi*t)*(x - 6)*(6*y - 5*x + 45) 
				+ c*x*pi*cos((pi*x)/6)*sin(pi*t)*(x - 6)*(5*x + 6*y - 45)) 
			- (2*exp(-abs(u)) + 1)*(
				2*c*sin(pi*t)*sin((pi*x)/6)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
				- 50*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6) 
				+ 10*c*x*sin(pi*t)*sin((pi*x)/6)*(6*y - 5*x + 45) 
				- 10*c*x*sin(pi*t)*sin((pi*x)/6)*(5*x + 6*y - 45) 
				+ 10*c*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45) 
				- 10*c*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(5*x + 6*y - 45) 
				+ (c*pi*cos((pi*x)/6)*sin(pi*t)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/3 
				+ (5*c*x*pi*cos((pi*x)/6)*sin(pi*t)*(x - 6)*(6*y - 5*x + 45))/3 
				- (5*c*x*pi*cos((pi*x)/6)*sin(pi*t)*(x - 6)*(5*x + 6*y - 45))/3 
				+ (c*x*pi*cos((pi*x)/6)*sin(pi*t)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/3 
				- (c*x*pi**2*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/36) 
			+ 2*exp(-abs(u))*sign(u)*(
				5*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45) 
				- 5*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(5*x + 6*y - 45) 
				+ c*x*sin(pi*t)*sin((pi*x)/6)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
				+ c*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
				+ (c*x*pi*cos((pi*x)/6)*sin(pi*t)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/6)**2 
			+ 2*exp(-abs(u))*sign(u)*(6*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45) 
				+ 6*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(5*x + 6*y - 45))*(
					5*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45) 
					- 5*c*x*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(5*x + 6*y - 45) 
					+ c*x*sin(pi*t)*sin((pi*x)/6)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
					+ c*sin(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
					+ (c*x*pi*cos((pi*x)/6)*sin(pi*t)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45))/6) 
			- 72*c*x*sin(pi*t)*sin((pi*x)/6)*(4*exp(-abs(u)) + 2)*(x - 6) 
			+ c*x*pi*cos(pi*t)*sin((pi*x)/6)*(x - 6)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
		)
	return np.ravel(f, order='F')

def simulate(degree, cuts, quadArgs, problemArgs={}):
	# Create model 
	geoArgs = {'name': 'TP', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	timeGeo  = part1D(createUniformOpenCurve(degree, 2**cuts, 1.), {'quadArgs': quadArgs})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 
	material.addConductivityDers(conductivityDersProperty, isIsotropic=False) 
	material.addCapacityDers(capacityDersProperty, isIsotropic=False) 

	# Block boundaries
	dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
	stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timeGeo.nbctrlpts_total])
	boundary = boundaryCondition(stnbctrlpts)
	boundary.add_DirichletConstTemperature(table=dirichlet_table)

	problem = stheatproblem(material, modelPhy, timeGeo, boundary)
	Fext = problem.compute_volForce(powerDensity, 
									{'position':problem.part.qpPhy, 
									'time':problem.time.qpPhy})
	
	uguess = np.zeros(np.prod(stnbctrlpts))
	# uguess[boundary.thdod] = 0.0
	isfull = problemArgs.get('isfull', False); isadaptive = problemArgs.get('isadaptive', True)
	output = problem.solveFourierSTHeatProblem(uguess, Fext, isfull=isfull, isadaptive=isadaptive)
	return problem, output

meshpartext = '_' + str(degree) + '_' + str(cuts) + '/'
subfolderfolder = folder + meshpartext 
if not os.path.isdir(subfolderfolder): os.mkdir(subfolderfolder)

if TODOSIMU:

	quadArgs = {'quadrule': 'iga', 'type': 'leg'}

	for [i, isadaptive], prefix1 in zip(enumerate([True]), ['inexact']):
		for [j, isfull], prefix2 in zip(enumerate([False, True]), ['picard', 'newton']):
	# for [i, isadaptive], prefix1 in zip(enumerate([False, True]), ['exact', 'inexact']):
	# 	for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
			prefix = prefix1 + '_' + prefix2 + '_'
			problemArgs = {'isfull':isfull, 'isadaptive':isadaptive}
			blockPrint()
			problem, output = simulate(degree, cuts, quadArgs, problemArgs=problemArgs)
			displacements = output['Solution']
			deltadisps 	  = output['Delta']
			resKrylovs    = output['KrylovRes']
			resNewtons    = output['NewtonRes']
			thresholds    = output['Threshold']
			L2error, L2relerror, L2normDelta = [], [], []
			enablePrint()
			print("-------------------")

			# for displacement, deltadisp in zip(displacements, deltadisps):
			# 	err, relerr  = problem.normOfError(displacement, normArgs={'type':'L2', 
			# 													'exactFunction':exactTemperature})
			# 	normdelta, _ = problem.normOfError(-deltadisp, normArgs={'type':'L2', 
			# 												'exactFunction':zeros})
			# 	L2error.append(err); L2relerror.append(relerr), L2normDelta.append(normdelta)
			# resKrylovclean = np.array([]); counter_list = [0]
			# for _ in resKrylovs: 
			# 	resKrylovclean = np.append(resKrylovclean, _[np.nonzero(_)])
			# 	counter_list.append(counter_list[-1] + len(_[np.nonzero(_)]))
			# enablePrint()
			# np.savetxt(subfolderfolder+prefix+'CumulKrylovRes'+extension, resKrylovclean)
			# np.savetxt(subfolderfolder+prefix+'Inner_loops'+extension, counter_list)
			# np.savetxt(subfolderfolder+prefix+'NewtonRes'+extension, resNewtons)
			# np.savetxt(subfolderfolder+prefix+'L2error'+extension, L2error)
			# np.savetxt(subfolderfolder+prefix+'L2relerror'+extension, L2relerror)
			# np.savetxt(subfolderfolder+prefix+'L2normDelta'+extension, L2normDelta)
			# np.savetxt(subfolderfolder+prefix+'threshold'+extension, np.array(thresholds))

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
				yy = L2relerror; xx = np.arange(0, len(L2relerror))
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

fig, ax = plt.subplots(figsize=(8, 6))
for [i, isadaptive], prefix1 in zip(enumerate([True, False]), ['inexact', 'exact']):
	for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
		l = j + i*2
		legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
		prefix = prefix1 + '_' + prefix2 + '_'
		increment = np.loadtxt(subfolderfolder+prefix+'L2normDelta'+extension)

		yy = increment/increment[0]; xx = np.arange(0, len(yy))
		ylabel = 'Relative norm of increment'
		xlabel = 'Number of outer iterations'

		ax.semilogy(xx, yy, label=legendname, marker=marker_list[l], linestyle=linestyle_list[l])
		ax.set_xlim(right=10, left=0)
		ax.set_ylim(top=10, bottom=1e-11)
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.legend()
		fig.tight_layout()
		fig.savefig(folder+'L2normincrement'+'_'+str(degree)+str(cuts)+'.pdf')

