
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_job import heatproblem
from pysrc.lib.lib_stjob import stheatproblem
from pysrc.lib.lib_1d import heatproblem1D
from numpy import pi, sin, cos, abs, exp, sign

IS1DIM = False
ISLINEAR = False
CST = 100
CUTS_TIME = 7

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	capacity = np.ones(shape=np.shape(temperature))
	return capacity

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: Kprop1d = 2*np.ones(shape=np.shape(temperature))
	# else: Kprop1d = 1.0 + 2.0*exp(-abs(temperature))
	else: Kprop1d = 1.0 + 2.0*exp(-(sin(0.1*temperature))**2)

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
	# else: y = -2.0*sign(temperature)*exp(-abs(temperature))
	else: y = -0.2*exp(-(sin(0.1*temperature))**2)*sin(0.2*temperature)
	return y

def exactTemperature_inc(args:dict):
	t = args['time']
	if IS1DIM: x = args['position']
	else: x = args['position'][0, :]
	u = CST*sin(2*pi*x)*sin(pi/2*t)
	return u

def exactTemperature_spt(qpPhy):
	x = qpPhy[0, :]; y = qpPhy[1, :]; t = qpPhy[2, :]
	u = CST*sin(2*pi*x)*sin(pi/2*t)
	return u

def powerDensity_inc(args:dict):
	t = args['time']
	if IS1DIM: x = args['position']
	else: x = args['position'][0, :]

	if ISLINEAR:
		f = (
			(CST*pi*cos((pi*t)/2)*sin(2*pi*x))/2 
			+ 8*CST*pi**2*sin((pi*t)/2)*sin(2*pi*x)
		)
	else: 
		# u = CST*sin((pi*t)/2)*sin(2*pi*x)
		# f = (
		# 	(CST*pi*cos((pi*t)/2)*sin(2*pi*x))/2 
		# 	+ 4*CST*pi**2*sin((pi*t)/2)*sin(2*pi*x)*(2*exp(-abs(u)) + 1) 
		# 	+ 8*CST**2*pi**2*sign(u)*exp(-abs(u))*cos(2*pi*x)**2*sin((pi*t)/2)**2

		# )

		u = CST*sin((pi*t)/2)*sin(2*pi*x)
		f = ((CST*pi*cos((pi*t)/2)*sin(2*pi*x))/2 
		+ 4*CST*pi**2*sin((pi*t)/2)*sin(2*pi*x)*(2*exp(-sin((u)/10)**2) + 1) 
		+ (8*CST**2*pi**2*cos((u)/10)*sin((u)/10)*exp(-sin((u)/10)**2)*cos(2*pi*x)**2*sin((pi*t)/2)**2)/5
		)
	return f

def powerDensity_spt(args:dict):
	timespan = args['time']
	position = args['position']
	nc_sp = np.size(position, axis=1); nc_tm = np.size(timespan); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = timespan[i]
		f[:, i] = powerDensity_inc(args={'time':t, 'position':position})
	return np.ravel(f, order='F')

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

def simulate_spacetime(degree, cuts, powerdensity=None):
	geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, 2, 1])}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs={'quadrule': 'iga'})
	time_spt = part1D(createUniformCurve(degree, int(2**CUTS_TIME), 1.0), {'quadArgs': {'quadrule': 'iga'}})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	tmptable = np.zeros((3, 2)); tmptable[0, :] = 1; tmptable[-1, 0] = 1
	sptnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], time_spt.nbctrlpts_total])
	boundary_spt = boundaryCondition(sptnbctrlpts)
	boundary_spt.add_DirichletConstTemperature(table=tmptable)

	# Transient model
	problem_spt = stheatproblem(material, modelPhy, time_spt, boundary_spt)

	# Add external force
	Fext = problem_spt.compute_volForce(powerdensity, 
									{'position':problem_spt.part.qpPhy, 
									'time':problem_spt.time.qpPhy})
	
	# Solve
	Tinout = np.zeros(np.prod(sptnbctrlpts))
	problem_spt._itersNL = 100; problem_spt._thresNL = 1e-7
	problem_spt.solveFourierSTHeatProblem(Tinout=Tinout, Fext=Fext, isfull=False, isadaptive=True)

	return problem_spt, time_spt, Tinout