
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_job import heatproblem
from pysrc.lib.lib_1d import heatproblem1D
from numpy import pi, sin, cos, abs, exp, sign

IS1DIM = False
ISLINEAR = True
CST = 100
CUTS_TIME = 7

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	capacity = np.ones(shape=np.shape(temperature))
	return capacity

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: Kprop1d = 2*np.ones(shape=np.shape(temperature))
	else: Kprop1d = (1.0 + 2.0*exp(-abs(temperature)))
	# else: y = (1.0 + 2.*np.exp(-(np.sin(0.01*temperature))**2))

	if IS1DIM: 
		return Kprop1d
	else:
		identity = np.array([[1., 0.0],[0.0, 1.0]])
		Kprop2d  = np.zeros((2, 2, len(temperature)))
		for i in range(2): 
			for j in range(2):
				Kprop2d[i, j, :] = identity[i, j]*Kprop1d
		return Kprop2d

def conductivityDersProperty(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: y = np.zeros(shape=np.shape(temperature))
	else: y = -2.0*sign(temperature)*exp(-abs(temperature))
	# else: y = ...
	return y

def exactTemperature(args:dict):
	t = args['time']
	if IS1DIM: x = args['position']
	else: x = args['position'][0, :]
	u = CST*sin(2*pi*x)*sin(pi/2*t)
	return u

def powerDensity(args:dict):
	t = args['time']
	if IS1DIM: x = args['position']
	else: x = args['position'][0, :]

	if ISLINEAR:
		f = (
			(CST*pi*cos((pi*t)/2)*sin(2*pi*x))/2 
			+ 8*CST*pi**2*sin((pi*t)/2)*sin(2*pi*x)
		)
	else: 
		u = CST*sin((pi*t)/2)*sin(2*pi*x)
		f = (
			(CST*pi*cos((pi*t)/2)*sin(2*pi*x))/2 
			+ 4*CST*pi**2*sin((pi*t)/2)*sin(2*pi*x)*(2*exp(-2*abs(u)) + 1) 
			+ 16*CST**2*pi**2*sign(u)*exp(-2*abs(u))*cos(2*pi*x)**2*sin((pi*t)/2)**2
		)
		# u = ...
		# f = (
		# )
	return f

def simulate_incremental(degree, cuts, powerdensity=None, is1dim=False):

	# Create geometry
	if is1dim:
		geometry = createUniformCurve(degree, int(2**cuts), 1.0)
		modelPhy = part1D(geometry, kwargs={'quadArgs':{'quadrule': 'iga'}})
	else:
		geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': np.array([cuts, 1, 1])}
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs={'quadrule': 'iga'})

	time_inc = np.linspace(0, 1.0, int(2**CUTS_TIME)+1) 

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	dirichlet_table = np.zeros((2, 2)); dirichlet_table[0, :] = 1
	boundary_inc = boundaryCondition(modelPhy.nbctrlpts)
	boundary_inc.add_DirichletConstTemperature(table=dirichlet_table)

	# Transient model
	if is1dim: problem_inc = heatproblem1D(material, modelPhy, boundary_inc)
	else: problem_inc = heatproblem(material, modelPhy, boundary_inc)

	# Add external force 
	Fext_list = np.zeros((problem_inc.part.nbctrlpts_total, len(time_inc)))
	for i, t in enumerate(time_inc):
		Fext_list[:, i] = problem_inc.compute_volForce(powerdensity, 
							args={'position':problem_inc.part.qpPhy, 'time':t})

	# Solve
	Tinout = np.zeros((modelPhy.nbctrlpts_total, len(time_inc)))
	problem_inc._itersNL = 100; problem_inc._thresNL = 1e-7
	problem_inc.solveFourierTransientProblem(Tinout=Tinout, Fext_list=Fext_list, 
											time_list=time_inc, alpha=0.5)
	return problem_inc, time_inc, Tinout