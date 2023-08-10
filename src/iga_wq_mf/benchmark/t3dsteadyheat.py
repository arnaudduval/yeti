"""
.. Test of steady heat transfer 3D
.. In this case we impose a function as a temperature over a surface
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import thermomat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import heatproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d2steadyheat/'
if not os.path.isdir(folder): os.mkdir(folder)

def conductivityProperty(P:list):
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((3, 3, np.size(P, axis=1)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j] 
	# Kprop[0, 0, :] += 0.75*np.cos(np.pi*y)
	# Kprop[1, 1, :] += 2*np.exp(-(z-0.5)**2)
	# Kprop[2, 2, :] += 2.5*np.cos(np.pi*x)**2
	return Kprop 

def powerDensity_rotRing(P: list):
	""" u = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]

	# # Isotropy
	# f = (8*x*y**4*np.sin(np.pi*z) + 8*x**3*y**2*np.sin(np.pi*z) 
	#     + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 1) + 16*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
	#     + 2*x*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	#     - x*y**2*np.pi**2*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
	# )

	# Anisotropy
	f = (16*x*y**4*np.sin(np.pi*z) 
	+ 8*x**2*y**3*np.sin(np.pi*z) 
	+ 8*x**3*y**2*np.sin(np.pi*z) 
	+ 2*y**3*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
	+ 2*y**3*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
	+ 26*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
	+ 4*x**2*y*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
	+ 26*x*y**2*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
	+ 4*x**2*y*np.sin(np.pi*z)*(x**2 + y**2 - 4) 
	+ 4*x*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	+ 2*y*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	+ (y**2*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4))/5 
	+ x*y**3*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 1) 
	+ x*y**3*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 4) 
	+ (2*x**2*y**2*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 1))/5 
	+ (2*x**2*y**2*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 4))/5 
	- 3*x*y**2*np.pi**2*np.sin(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4) 
	+ x*y*np.pi*np.cos(np.pi*z)*(x**2 + y**2 - 1)*(x**2 + y**2 - 4)
	)

	return f

def exactTemperature_rotRing(P: list):
	" T = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*sin(pi*z) "
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	u = -(x**2 + y**2 - 1)*(x**2 + y**2 - 4)*x*(y**2)*np.sin(np.pi*z)
	return u

# Set global variables
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-15}
degree_list = np.array([2, 3, 4, 6, 8])
cuts_list   = np.arange(2, 9)

for quadrule, quadtype in zip(['wq', 'iga'], [1, 'leg']):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	error_list  = np.ones(len(cuts_list))
	fig, ax = plt.subplots(figsize=(8, 4))

	for i, degree in enumerate(degree_list):
		for j, cuts in enumerate(cuts_list):
			geoArgs = {'name': 'RQA', 'degree': degree*np.ones(3, dtype=int), 
						'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			}
			blockPrint()
			material = thermomat()
			material.addConductivity(conductivityProperty, isIsotropic=False)				
			modelGeo = Geomdl(geoArgs)
			modelIGA = modelGeo.getIGAParametrization()
			modelPhy = part(modelIGA, quadArgs=quadArgs)

			# Set Dirichlet boundaries
			boundary = boundaryCondition(modelPhy.nbctrlpts)
			boundary.add_DirichletConstTemperature(table=np.ones((3, 2), dtype=int))
			enablePrint()

		# 	# Solve elastic problem
		# 	problem = heatproblem(material, modelPhy, boundary)
		# 	problem.addSolverConstraints(solverArgs=solverArgs)
		# 	Fext = problem.eval_volForce(powerDensity_rotRing)
		# 	temperature = problem.solveSteadyHeatProblemFT(Fext=Fext)[0]
		# 	error_list[j] = problem.L2NormOfError(exactTemperature_rotRing, temperature)

		# nbctrlpts_list = (2**cuts_list+degree)**2
		# ax.loglog(nbctrlpts_list, error_list, marker=markerSet[i], label='degree '+r'$p=\,$'+str(degree))

		# if str(quadArgs['quadrule']) == 'iga':
		# 	slope = np.polyfit(np.log10(nbctrlpts_list[1:7]),np.log10(error_list[1:7]), 1)[0]
		# 	slope = round(slope, 1)
		# 	annotation.slope_marker((nbctrlpts_list[-3], error_list[-3]), slope, 
		# 							poly_kwargs={'facecolor': (0.73, 0.8, 1)})
			
		# ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{L_2(\Omega)}}{||u||_{L_2(\Omega)}}$')
		# ax.set_xlabel('Total number of DOF')
		# ax.set_ylim(top=1e0, bottom=1e-15)
		# ax.set_xlim(left=10, right=1e5)

		# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		# fig.tight_layout()
		# fig.savefig(folder + 'FigThickRing_' + str(quadArgs['quadrule']) +'.png')
