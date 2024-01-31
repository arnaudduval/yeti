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

def conductivityProperty(P:list):
	Kref  = np.array([[1, 0.5],[0.5, 2]])
	Kprop = np.zeros((2, 2, np.size(P, axis=1)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j] 
	return Kprop 

def powerDensity_quartCircle(P: list):
	""" u = sin(pi*x)*sin(pi*y)*sin(pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
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
degree_list = np.array([1, 2, 3, 4, 6, 8])
cuts_list   = np.arange(2, 9)

for quadrule, quadtype in zip(['wq', 'wq', 'iga'], [1, 2, 'leg']):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	error_list = np.ones(len(cuts_list))
	fig, ax = plt.subplots(figsize=(8, 4))

	for i, degree in enumerate(degree_list):
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

			# Set Dirichlet boundaries
			boundary = boundaryCondition(modelPhy.nbctrlpts)
			boundary.add_DirichletConstTemperature(table=np.ones((2, 2), dtype=int))
			enablePrint()

			# Solve elastic problem
			problem = heatproblem(material, modelPhy, boundary)
			Fext = problem.compute_volForce(powerDensity_quartCircle)
			temperature = problem.solveSteadyHeatProblem(Fext=Fext)[0]
			error_list[j], _ = problem.normOfError(temperature, normArgs={'type':'semiH1',
												'exactFunction':exactTemperature_quartCircle,
												'exactFunctionDers':exactDiffTemperature_quartCircle})

		nbctrlpts_list = (2**cuts_list+degree)**2
		ax.loglog(nbctrlpts_list, error_list, marker=MARKERLIST[i], label='degree '+r'$p=\,$'+str(degree))

		if quadrule == 'iga':
			slope = np.polyfit(np.log10(nbctrlpts_list[1:7]),np.log10(error_list[1:7]), 1)[0]
			slope = round(slope, 1)
			annotation.slope_marker((nbctrlpts_list[-3], error_list[-3]), slope, 
									poly_kwargs={'facecolor': (0.73, 0.8, 1)})
			
		ax.set_ylabel(r'$||u - u^h||_{L^2(\Omega)}$')
		# ax.set_ylabel(r'$|u - u^h|_{H^1(\Omega)}$')
		ax.set_xlabel('Total number of DOF')
		ax.set_ylim(top=1e1, bottom=1e-15)
		ax.set_xlim(left=10, right=1e5)

		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigConvergenceSemiH1' +  geoName + '_' + quadrule + str(quadtype) +'.pdf')
