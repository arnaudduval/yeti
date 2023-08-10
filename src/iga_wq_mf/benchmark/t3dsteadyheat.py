"""
.. Test of steady heat transfer 3D
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
folder = os.path.dirname(full_path) + '/results/d3steadyheat/'
if not os.path.isdir(folder): os.mkdir(folder)

def conductivityProperty(P:list):
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((3, 3, np.size(P, axis=1)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j] 
	return Kprop 

def powerDensity_thickRing(P: list):
	""" u = sin(pi*x)*sin(pi*y)*sin(pi*z)*(x**2+y**2-1)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]; y = P[1, :]; z = P[2, :] 

	f = (6*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
    - np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
	- (np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*z)*np.sin(np.pi*y)*(x**2 + y**2 - 1))/5 
	- (np.pi**2*np.cos(np.pi*y)*np.cos(np.pi*z)*np.sin(np.pi*x)*(x**2 + y**2 - 1))/2 
	- 4*x*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 
	- 2*x*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*np.sin(np.pi*z) 
	- (2*x*np.pi*np.cos(np.pi*z)*np.sin(np.pi*x)*np.sin(np.pi*y))/5 
	- 2*y*np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 
	- 8*y*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*np.sin(np.pi*z) 
	- y*np.pi*np.cos(np.pi*z)*np.sin(np.pi*x)*np.sin(np.pi*y) 
	- 6*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)
	)

	return f

def heatFlux_thickRing(P:list):
	x = P[0, :]; y = P[1, :]; z = P[2, :]; nnz = np.size(P, axis=1)
	r_square = x**2 + y**2
	theta = np.arcsin(y/np.sqrt(r_square))

	flux = np.zeros((3, nnz))
	flux[0, :] = (2*x*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 
				+ y*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 
				+ np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
				+ (np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*np.sin(np.pi*z)*(x**2 + y**2 - 1))/2 
				+ (np.pi*np.cos(np.pi*z)*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1))/10
	)
	flux[1, :] = (x*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 
				+ 4*y*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z) 
				+ (np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)*(x**2 + y**2 - 1))/2 
				+ 2*np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*np.sin(np.pi*z)*(x**2 + y**2 - 1) 
				+ (np.pi*np.cos(np.pi*z)*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1))/4
	)
	flux[2, :] = ((x*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z))/5 
				+ (y*np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z))/2 
				+ (np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)*(x**2 + y**2 - 1))/10 
				+ (np.pi*np.cos(np.pi*y)*np.sin(np.pi*x)*np.sin(np.pi*z)*(x**2 + y**2 - 1))/4 
				+ 3*np.pi*np.cos(np.pi*z)*np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 - 1)
	)

	F = np.zeros(nnz)
	F = flux[0, :]*np.cos(theta) + flux[1, :]*np.sin(theta)

	return F

def exactTemperature_thickRing(P: list):
	" T = sin(pi*x)*sin(pi*y)*sin(pi*z)*(x**2+y**2-1) "
	x = P[0, :]
	y = P[1, :]
	z = P[2, :]
	u = np.sin(np.pi*x)*np.sin(np.pi*y)*np.sin(np.pi*z)*(x**2 + y**2 - 1)
	return u

# Set global variables
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-15}
degree_list = np.array([2, 3, 4, 6])
cuts_list   = np.arange(2, 7)

for quadrule, quadtype in zip(['wq', 'iga'], [1, 'leg']):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	error_list = np.ones(len(cuts_list))
	fig, ax = plt.subplots(figsize=(8, 4))

	for i, degree in enumerate(degree_list):
		for j, cuts in enumerate(cuts_list):
			geoArgs = {'name': 'TR', 'degree': degree*np.ones(3, dtype=int), 
						'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			}
			blockPrint()
			material = thermomat()
			material.addConductivity(conductivityProperty, isIsotropic=False)				
			modelGeo = Geomdl(geoArgs)
			modelIGA = modelGeo.getIGAParametrization()
			modelPhy = part(modelIGA, quadArgs=quadArgs)

			# Set Dirichlet boundaries
			table = np.ones((3, 2), dtype=int); table[0, 1] = 0
			boundary = boundaryCondition(modelPhy.nbctrlpts)
			boundary.add_DirichletConstTemperature(table=table)
			enablePrint()

			# Solve elastic problem
			problem = heatproblem(material, modelPhy, boundary)
			problem.addSolverConstraints(solverArgs=solverArgs)
			Fsurf = problem.compute_surfForce(heatFlux_thickRing, nbFacePosition=1)
			Fvol  = problem.compute_volForce(powerDensity_thickRing)
			Fext  = Fvol + Fsurf
			temperature = problem.solveSteadyHeatProblemFT(Fext=Fext)[0]
			error_list[j] = problem.L2NormOfError(exactTemperature_thickRing, temperature)

		nbctrlpts_list = (2**cuts_list+degree)**3
		ax.loglog(nbctrlpts_list, error_list, marker=markerSet[i], label='degree '+r'$p=\,$'+str(degree))
			
		ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{L_2(\Omega)}}{||u||_{L_2(\Omega)}}$')
		ax.set_xlabel('Total number of DOF')
		ax.set_ylim(top=1e0, bottom=1e-15)
		ax.set_xlim(left=1e2, right=1e6)

		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'FigThickRing2_' + str(quadArgs['quadrule']) +'.png')
