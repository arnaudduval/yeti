"""
.. Test of 1D transient heat transfer 
.. Author: Joaquin Cornejo
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_1d import heatproblem1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/tansient1d/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
extension = '.dat'
FIG_CASE  = 1
ISLINEAR  = True
c = 0.001

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: y = np.ones(len(temperature))
	else: y = (1.0 + 2.0*np.exp(-np.abs(temperature)))
	return y

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: y = np.ones(len(temperature))
	else: y = (1.0 + np.exp(-np.abs(temperature)))
	return y

def exactTemperature(args:dict):
	x = args.get('position')
	t = args.get('time')
	u = c*np.sin(2*np.pi*x)*np.sin(2*np.pi*t)
	return u

def powerDensity(args:dict):
	x = args['position']; t = args['time']
	if ISLINEAR:
		f = (2*c*np.pi*np.cos(2*np.pi*t)*np.sin(2*np.pi*x) 
			+ 4*c*np.pi**2*np.sin(2*np.pi*t)*np.sin(2*np.pi*x)
			)
	else: 
		u = c*np.sin(2*np.pi*t)*np.sin(2*np.pi*x)
		f = (4*c*np.pi**2*np.sin(2*np.pi*t)*np.sin(2*np.pi*x)*(2*np.exp(-np.abs(u)) + 1) 
			+ 2*c*np.pi*np.cos(2*np.pi*t)*np.sin(2*np.pi*x)*(np.exp(-np.abs(u)) + 1) 
			+ 8*c**2*np.pi**2*np.exp(-np.abs(u))*np.sign(u)*np.cos(2*np.pi*x)**2*np.sin(2*np.pi*t)**2
			)
	return f

def simulate(degree, cuts, cuts_time=None):

	# Create geometry
	length = 1.0
	nbel   = int(2**cuts)
	geometry = createUniformCurve(degree, nbel, length)
	modelPhy = part1D(geometry, kwargs={'quadArgs': {'quadrule': 'iga'}})

	if cuts_time is None: cuts_time = np.copy(cuts)
	timespan = 1.0
	nbsteps  = int(2**cuts_time)
	time_list = np.linspace(0, timespan, nbsteps+1)

	# Add material
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False)
	material.addCapacity(capacityProperty, isIsotropic=False)

	# Add boundary condition
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	boundary.add_DirichletConstTemperature(table=np.array([[1, 1]]))
	temperature = np.zeros((modelPhy.nbctrlpts_total, len(time_list)))

	# Create heat problem
	problem = heatproblem1D(material, modelPhy, boundary)

	# Define external force	
	Fext_list = np.zeros((modelPhy.nbctrlpts_total, len(time_list)))
	for i, t in enumerate(time_list):
		Fext_list[:, i] = problem.compute_volForce(powerDensity, args={'position':problem.part.qpPhy, 'time':t})

	# Solve
	problem.solveFourierTransientProblem(temperature, Fext_list, time_list[:3], isLumped=False, alpha=1.0)
	return problem, time_list, temperature

if FIG_CASE == 1:
	lastsufix = 'linear' if ISLINEAR else 'nonlin'
	degree_list = np.array([1, 2, 3, 4])
	cuts_list   = np.arange(1, 6)

	error_list = np.zeros((len(degree_list), len(cuts_list), 2**np.max(cuts_list)))
	for j, cuts in enumerate(cuts_list):
		for i, degree in enumerate(degree_list):
			problem, time_list, output = simulate(degree, cuts, cuts_time=np.max(cuts_list))
			for k, t in enumerate(time_list[1:-1]):
				_, error_list[i, j, k] = problem.normOfError(output[:, k+1], 
															normArgs={'type':'L2',
																	'exactFunction':exactTemperature,
																	'exactExtraArgs':{'time':t}})
	np.save(folder + 'incrementalheat', error_list)

	error_list = np.load(folder + 'incrementalheat.npy')
	for k in range(np.size(error_list, axis=2)):
		fig, ax = plt.subplots(figsize=(9, 6))
		for i, degree in enumerate(degree_list):
			color = COLORLIST[i]
			ax.loglog(2**cuts_list, error_list[i, :, k], color=color, marker='o', markerfacecolor='w',
						markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))
		ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{L_2(\Omega)}}{||u||_{L_2(\Omega)}}$')
		ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
		ax.set_ylim(top=1e1, bottom=1e-5)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(folder + 'steps/FigConvergenceIncrHeat' + str(k+1) +  '.pdf')
		plt.close(fig)