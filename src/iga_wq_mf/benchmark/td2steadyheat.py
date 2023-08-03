"""
.. Test of elastoplasticity 2D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa (200e3)
..      - Length : mm
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

def setKprop(P:list):
	x = P[0, :]
	y = P[1, :]
	Kref  = np.array([[1, 0.5],[0.5, 2]])
	Kprop = np.zeros((2, 2, len(x)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j] 
	return Kprop 

def powden_thickring(P: list):
	""" u = sin(pi*x)*sin(pi*y)*sin(pi*z)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]

	# Anisotropy
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

def temperature_exact(P: list):
	""" u = sin(5*pi*x)*sin(5*pi*y)*(x**2+y**2-1)*(x**2+y**2-4)
		f = -div(lambda * grad(u))
	"""
	x = P[0, :]
	y = P[1, :]

	# Anisotropy
	u = np.sin(np.pi*x)*np.sin(np.pi*y)*(x**2 + y**2 -1)*(x**2 + y**2 - 4)
	return u

# Set global variables
energy_th = 83.43708849059514
# quadArgs = {'quadrule': 'iga', 'type': 'leg'}
quadArgs = {'quadrule': 'wq', 'type': 2}
solverArgs = {'nbIterationsPCG':150, 'PCGThreshold':1e-16}

# Create model 
degree_list = np.array([2, 3, 4, 6])
cuts_list   = np.arange(2, 9)
error = np.ones(len(cuts_list))

fig, ax  = plt.subplots(figsize=(8, 4))
for i, degree in enumerate(degree_list):
	for j, cuts in enumerate(cuts_list):
		geoArgs = {'name': 'QA', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
					'extra':{'Rin':1.0, 'Rex':2.0}
		}
		blockPrint()
		material = thermomat()
		material.addConductivity(setKprop, isIsotropic=False)				
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)

		# Set Dirichlet boundaries
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		boundary.add_DirichletTemperature(table=np.ones((2, 2), dtype=int))
		enablePrint()

		# Solve elastic problem
		problem = heatproblem(material, modelPhy, boundary)
		problem.addSolverConstraints(solverArgs=solverArgs)
		Fn      = problem.eval_volForce(powden_thickring, indi=boundary.thdof)
		tmp     = problem.solveSteadyHeatProblemFT(b=Fn)[0]
		temperature = np.zeros(modelPhy.nbctrlpts_total)
		temperature[boundary.thdof] = tmp
		KT = problem.eval_mfConductivity(temperature, args=problem.part.qpPhy, table=np.zeros((2, 2), dtype=int))
		# error[j] = np.abs(energy_th - np.dot(KT, temperature))/energy_th*100
		error[j] = problem.L2NormOfError(temperature_exact, temperature)*100

	nbctrlpts = (2**cuts_list+degree)**2
	ax.loglog(nbctrlpts, error, marker=markerSet[i], label='degree p='+str(degree))

	if str(quadArgs['quadrule']) == 'iga':
		slope = np.polyfit(np.log10(nbctrlpts[:4]),np.log10(error[:4]), 1)[0]
		slope = round(slope, 1)
		annotation.slope_marker((nbctrlpts[2], error[2]), slope, 
								poly_kwargs={'facecolor': (0.73, 0.8, 1)})
	
	# ax.set_ylabel('Relative error of energy ' + r'$\frac{|U-U^{h}|}{|U|}$')
	ax.set_ylabel('Relative L2 norm error')
	ax.set_xlabel('Total number of DOF')
	ax.set_ylim(top=1e1, bottom=1e-14)
	ax.set_xlim(left=10, right=1e5)

	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.tight_layout()
	fig.savefig(folder + 'FigThickRing2_' + str(quadArgs['quadrule']) +'.png')
