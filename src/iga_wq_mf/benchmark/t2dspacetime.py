from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem

def conductivityProperty(args):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.0, 0.0],[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
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
	position = args['Position']; timespan = args['Time']
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

def simulate(degree, cuts, quadArgs):
	# Create model 
	geoArgs = {'name': 'sq', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Create time span
	crv = createUniformCurve(degree, 2**cuts, 1.0)
	timespan = part1D(crv, {'quadArgs':quadArgs})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
	stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timespan.nbctrlpts])
	boundary = boundaryCondition(stnbctrlpts)
	boundary.add_DirichletConstTemperature(table=dirichlet_table)

	problem = stheatproblem(material, modelPhy, timespan, boundary)
	# External heat force
	Fext = problem.compute_volForce(powerDensity, {'Position':problem.part.qpPhy, 'Time':problem.time.qpPhy})
	u_guess = np.zeros(np.prod(stnbctrlpts)); u_guess[boundary.thdod] = 0.0
	output = problem.solveFourierSTHeatProblem(u_guess, Fext, isfull=False, isadaptive=True)
	return problem, output

# ---------------------
# Transient model
# ---------------------
normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}

degree_list = np.array([1, 2, 3, 4])
cuts_list   = np.arange(1, 7)

fig, ax = plt.subplots(figsize=(8, 6))
# for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 1, 2], [normalPlot, onlyMarker1, onlyMarker2]):
for quadrule, quadtype, plotpars in zip(['iga'], ['leg'], [normalPlot]):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	error_list = np.ones(len(cuts_list))

	for i, degree in enumerate(degree_list):
		color = COLORLIST[i]
		for j, cuts in enumerate(cuts_list):
			nbels = 2**cuts_list
			problem, output = simulate(degree, cuts, quadArgs)
			displacement = output['Solution'][-1]
			_, error_list[j] = problem.normOfError(displacement, normArgs={'type':'L2', 
												'exactFunction':exactTemperature,})
			
		if quadrule == 'iga': 
			ax.loglog(nbels, error_list, label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
						markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			
			slope = np.polyfit(np.log10(nbels),np.log10(error_list), 1)[0]
			slope = round(slope, 1)
			annotation.slope_marker((nbels[-2], error_list[-2]), slope, 
							poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
		else: 
			ax.loglog(nbels, error_list, color=color, marker=plotpars['marker'], markerfacecolor='w',
					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
		
		fig.savefig(folder + 'FigSPTLinearConvergenceL2' + '.pdf')

# ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
# 				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ 2")
# ax.loglog([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
# 		markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label="IGA-WQ 4")

ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$')
ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
ax.set_ylim(top=1e1, bottom=1e-11)
ax.legend()
fig.tight_layout()
fig.savefig(folder + 'FigSPTLinearConvergenceL2' + '.pdf')