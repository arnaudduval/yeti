"""
.. Test of steady heat transfer 2D
.. The geometry is a quart of a annulus
.. All the boundary conditions are considered Dirichlet-like 
.. The heat source is computed by f = -grad(k grad T) where T is a given function
.. We compute the relative L2 norm using the exact solution T
.. The convergence curves are traced for IGA-Legendre and IGA-WQ
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import heatproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2transferheat/'
if not os.path.isdir(folder): os.mkdir(folder)

def conductivityProperty(args:dict):
	P = args.get('position')
	Kref  = np.array([[1, 0.5],[0.5, 2]])
	Kprop = np.zeros((2, 2, np.size(P, axis=1)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j] 
	return Kprop 

def powerDensity_quartCircle(args:dict):
	""" u = sin(pi*x)*sin(pi*y)*sin(pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	P = args.get('position')
	x = P[0, :]
	y = P[1, :]

	f = (3*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- 16*y**2*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 8*x*y*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 4) 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1) 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 4) 
	- 8*x**2*np.sin(np.pi*x)*np.sin(np.pi*y)
	)

	return f

def exactTemperature_quartCircle(P: list):
	""" u = sin(pi*x)*sin(pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]

	u = np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 -1)*(x**2 + y**2 - 4)
	return u

def exactDiffTemperature_quartCircle(P: list):
	""" u = sin(pi*x)*sin(pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
		uders = grad(u)
	"""
	x = P[0, :]
	y = P[1, :]

	uders = np.zeros((1, 2, np.size(P, axis=1)))

	uders[0, 0, :] = (2*x*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
				+ 2*x*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4) 
				+ np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
				)

	uders[0, 1, :] = (2*y*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1) 
				+ 2*y*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 4)
				+ np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
				)

	return uders

# Set global variables
geoName = 'QA'
degree_list = np.array([1, 2, 3, 4, 5])
cuts_list   = np.arange(1, 9)

normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}

fig, ax = plt.subplots(figsize=(8, 7))
figname = folder + 'FigHeatConvergenceAllL2' + '.pdf'

for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 2, 1], [normalPlot, onlyMarker1, onlyMarker2]):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	error_list = np.ones(len(cuts_list))
	for i, degree in enumerate(degree_list):
		meshparam = np.ones(len(cuts_list))
		color = COLORLIST[i]
		for j, cuts in enumerate(cuts_list):
			geoArgs = {'name': geoName, 'degree': degree*np.ones(3, dtype=int), 
						'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
						'extra':{'Rin':1.0, 'Rex':2.0}
			}
			blockPrint()
			material = heatmat()
			material.addConductivity(conductivityProperty, isIsotropic=False)				
			modelGeo = Geomdl(geoArgs)
			modelIGA = modelGeo.getIGAParametrization()
			modelPhy = part(modelIGA, quadArgs=quadArgs)
			meshparam[j] = modelPhy.compute_global_mesh_parameter()	

			# Set Dirichlet boundaries
			boundary = boundaryCondition(modelPhy.nbctrlpts)
			boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))
			enablePrint()

			# Solve elastic problem
			start=time.time()
			problem = heatproblem(material, modelPhy, boundary)
			Fext = problem.compute_volForce(powerDensity_quartCircle)
			temperature = np.zeros(problem.part.nbctrlpts_total)
			problem.solveFourierSteadyProblem(temperature, Fext=Fext)
			stop = time.time()
			error_list[j], _ = problem.normOfError(temperature, normArgs={'type':'L2',
												'exactFunction':exactTemperature_quartCircle,
												'exactFunctionDers':exactDiffTemperature_quartCircle})

		if quadrule == 'iga': 
			ax.loglog(meshparam, error_list, label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
		else: 
			ax.loglog(meshparam, error_list, color=color, marker=plotpars['marker'], markerfacecolor='w',
				markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
		fig.savefig(figname)

ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ 4")
ax.loglog([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
		markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label="IGA-WQ 2")

ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Omega)}$')
# ax.set_ylabel(r'$\displaystyle ||u - u^h||_{H^1(\Omega)}$')
ax.set_xlabel('Mesh parameter ' + r'$h_{max}$')
ax.set_ylim(top=1e1, bottom=1e-10)
ax.set_xlim(left=1e-2, right=5)
ax.legend()
fig.tight_layout()
fig.savefig(figname)