
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve, createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_job import heatproblem
from pysrc.lib.lib_stjob import stheatproblem
from pysrc.lib.lib_1djob import heatproblem1D
from numpy import pi, sin, cos, abs, exp, sign, tanh

IS1DIM = False
ISLINEAR = False
NONLINCASE = 2
CST = 100
CUTS_TIME = 7

def nonlinearfunc(args:dict):
	temperature = args.get('temperature')
	if NONLINCASE==1: Kprop1d = 1.0 + 2.0*exp(-abs(temperature))
	if NONLINCASE==2: Kprop1d = 3.0 + 2.0*tanh(temperature/50)
	return np.atleast_2d(Kprop1d)

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	capacity = np.ones(shape=np.shape(temperature))
	return capacity

def capacityDersProperty(args:dict):
	temperature = args.get('temperature')
	capacity = np.zeros(shape=np.shape(temperature))
	return capacity

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	if ISLINEAR: Kprop1d = 2*np.ones(shape=np.shape(temperature))
	else: 
		if NONLINCASE==1: Kprop1d = 1.0 + 2.0*exp(-abs(temperature))
		if NONLINCASE==2: Kprop1d = 3.0 + 2.0*tanh(temperature/50)

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
	else: 
		if NONLINCASE==1: y = -2.0*sign(temperature)*exp(-abs(temperature))
		if NONLINCASE==2: y = (2.0/50)/(np.cosh(temperature))**2
	return y

def exactTemperatureSquare_inc(args:dict):
	t = args['time']
	if IS1DIM: x = args['position']
	else: x = args['position'][0, :]
	u = CST*sin(2*pi*x)*sin(pi/2*t)
	return u

def exactTemperatureQuad_inc(args:dict):
	t = args['time']
	if IS1DIM: raise Warning('Try higher dimension')
	x = args['position'][0, :]
	y = args['position'][1, :]
	u = CST*sin(pi*y)*sin(pi*(y+0.75*x-0.5)*(-y+0.75*x-0.5))*sin(5*pi*x)*sin(pi/2*t)
	return u

def exactTemperatureRing_inc(args:dict):
	t = args['time']
	if IS1DIM: raise Warning('Try higher dimension')
	x = args['position'][0, :]
	y = args['position'][1, :]
	u = -CST*sin(0.5*pi*(x**2+y**2-1.))*sin(pi*(x**2+y**2-0.25**2))*sin(x*y)*sin(pi/2*t)
	return u

def exactTemperatureSquare_spt(qpPhy):
	x = qpPhy[0, :]; t = qpPhy[2, :]
	u = CST*sin(2*pi*x)*sin(pi/2*t)
	return u

def exactTemperatureQuad_spt(qpPhy):
	x = qpPhy[0, :]; y=qpPhy[1, :]; t = qpPhy[2, :]
	u = CST*sin(pi*y)*sin(pi*(y+0.75*x-0.5)*(-y+0.75*x-0.5))*sin(5*pi*x)*sin(pi/2*t)
	return u

def powerDensitySquare_inc(args:dict):
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
		if NONLINCASE==1:
			f = (
				(CST*pi*cos((pi*t)/2)*sin(2*pi*x))/2 
				+ 4*CST*pi**2*sin((pi*t)/2)*sin(2*pi*x)*(2*exp(-abs(u)) + 1) 
				+ 8*CST**2*pi**2*sign(u)*exp(-abs(u))*cos(2*pi*x)**2*sin((pi*t)/2)**2

			)
		if NONLINCASE==2:
			f = ((CST*pi*cos((pi*t)/2)*sin(2*pi*x))/2 
				+ 4*CST*pi**2*sin((pi*t)/2)*sin(2*pi*x)*(2*tanh((u)/50) + 3) 
				+ (4*CST**2*pi**2*cos(2*pi*x)**2*sin((pi*t)/2)**2*(tanh((u)/50)**2 - 1))/25
			)

	return f

def powerDensitySquare_spt(args:dict):
	time = args['time']
	position = args['position']
	nc_sp = np.size(position, axis=1); nc_tm = np.size(time); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = time[i]
		f[:, i] = powerDensitySquare_inc(args={'time':t, 'position':position})
	return np.ravel(f, order='F')

def exportTimeDependentMaterial(time_list, temperature=None, fields=None, geoArgs=None, folder=None):
	assert fields is not None, 'Add material fields'
	if geoArgs is None: geoArgs = {'name': 'SQ', 'degree': 4*np.ones(3, dtype=int), 
						'nb_refinementByDirection': 5*np.ones(3, dtype=int)}
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs={'quadrule': 'iga'})
	extraArgs = {}
	extraArgs['position'] = modelPhy.interpolateMeshgridField()[0]
	for i, tm in enumerate(time_list):
		extraArgs['time'] = tm
		if temperature is not None:	extraArgs['temperature'] = temperature(extraArgs)
		modelPhy.exportResultsCP(fields=fields, extraArgs=extraArgs, 
						folder=folder, name='out_'+str(i))
	return

def createAsymmetricalKnotvector(level, xasym=0.1):

	def discretize(array, level):
		n = len(array) 
		newarray = np.copy(array)
		if level == 2:
			for i in range(n-1): 
				newarray = np.append(newarray, (array[i]+array[i+1])/2)
			newarray = np.sort(newarray)
		else:
			for _ in range(1, level): 
				tmp = discretize(newarray, 2)
				newarray = np.copy(tmp)
		return newarray
	
	assert level > 0, 'Not possible. Try higher level'
	assert xasym < 0.25, 'Try lower value'
	knotvector = np.array([0.0, xasym, 0.5-xasym/2, 0.5+xasym/2, 1.0-xasym, 1.0]) # This is level 1
	if level > 1: 
		tmp = discretize(knotvector, level)
		knotvector = np.copy(tmp)

	return knotvector

def createAsymmetricalCurve(degree, level, length, xasym=0.05):
	crv = BSpline.Curve()
	crv.degree  = degree
	crv.ctrlpts = [[i*length/degree, 0.0] for i in range(degree+1)]
	crv.knotvector = createUniformKnotvector_Rmultiplicity(degree, 1)
	knotList = createAsymmetricalKnotvector(level, xasym=xasym)
	for knot in knotList[1:-1]:
		operations.insert_knot(crv, [knot], [1])
	return crv

def simulate_incremental(degree, cuts, dirichlet_table=None, powerdensity=None, is1dim=False, geoArgs=None):

	# Create geometry
	if is1dim:
		geometry = createUniformCurve(degree, int(2**cuts), 1.0)
		modelPhy = part1D(geometry, kwargs={'quadArgs':{'quadrule': 'iga'}})
	else:
		if geoArgs is None: 
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
	boundary_inc = boundaryCondition(modelPhy.nbctrlpts)
	if dirichlet_table is None: dirichlet_table = np.zeros((2, 2)); dirichlet_table[0, :] = 1
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
	problem_inc._itersNL = 50; problem_inc._thresNL = 1e-8
	problem_inc.solveFourierTransientProblem(Tinout=Tinout, Fext_list=Fext_list, 
											time_list=time_inc, alpha=0.5)
	return problem_inc, time_inc, Tinout

def simulate_spacetime(degree, cuts, dirichlet_table=None, powerdensity=None, degree_spt=None, 
					isfull=False, isadaptive=False, geoArgs=None, outputArgs=None):
	if geoArgs is None: 
		geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': np.array([cuts, 1, 1])}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs={'quadrule': 'iga'})
	if degree_spt is None: degree_spt = 2
	time_spt = part1D(createUniformCurve(degree_spt, int(2**CUTS_TIME), 1.0), {'quadArgs': {'quadrule': 'iga'}})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	sptnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], time_spt.nbctrlpts_total])
	boundary_spt = boundaryCondition(sptnbctrlpts)
	if dirichlet_table is None: dirichlet_table = np.zeros((3, 2)); dirichlet_table[0, :] = 1; dirichlet_table[-1, 0] = 1
	boundary_spt.add_DirichletConstTemperature(table=dirichlet_table)

	# Transient model
	problem_spt = stheatproblem(material, modelPhy, time_spt, boundary_spt)

	# Add external force
	Fext = problem_spt.compute_volForce(powerdensity, 
									{'position':problem_spt.part.qpPhy, 
									'time':problem_spt.time.qpPhy})
	
	# Solve
	Tinout = np.zeros(np.prod(sptnbctrlpts))
	problem_spt._itersNL = 50; problem_spt._thresNL = 1e-8
	output=problem_spt.solveFourierSTHeatProblem(Tinout=Tinout, Fext=Fext, isfull=isfull, isadaptive=isadaptive)
	outputArgs=deepcopy(output)
	return problem_spt, time_spt, Tinout