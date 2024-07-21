from thesis.Incremental.__init__ import *

# Select folder
folder = FOLDER2SAVE + '/welding/'
if not os.path.isdir(folder): os.mkdir(folder)

def conductivityProperty(args:dict):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.0],[0.0, 1.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			Kprop[i, j, :] = Kref[i, j]
	return Kprop 

def capacityProperty(args:dict):
	temperature = args['temperature']
	Cprop = np.ones(len(temperature))
	return Cprop

def powerDensity(args:dict):
	POWER = 20; RADIUS = 0.25; VELOCITY = 0.5
	position = args['position']; t = args['time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); f = np.zeros(nc_sp)
	if t > 16: return f
	rsquared = (x - VELOCITY*t)**2 + y**2
	f = POWER*np.exp(-rsquared/(RADIUS**2))
	return f

def simulate(degree, cuts, quadArgs, cuts_time=None):
	geoArgs = {'name': 'TP', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int),
				'extra':{'XY':np.array([[0.0, 0.0], [10.0, 0.0], [10.0, 5.0], [0.0, 5.0]])}}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	if cuts_time is None: cuts_time = np.copy(cuts)
	timespan  = 20
	nbsteps   = 2**cuts_time
	time_inc  = np.linspace(0, timespan, nbsteps+1) 

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	tmptable = np.zeros((2, 2)); tmptable[1, -1] = 1
	boundary_inc = boundaryCondition(modelPhy.nbctrlpts)
	boundary_inc.add_DirichletConstTemperature(table=tmptable)

	# Transient model
	problem_inc = heatproblem(material, modelPhy, boundary_inc)

	# Add external force 
	Fext_list = np.zeros((problem_inc.part.nbctrlpts_total, len(time_inc)))
	for i, t in enumerate(time_inc):
		Fext_list[:, i] = problem_inc.compute_volForce(powerDensity, args={'position':problem_inc.part.qpPhy, 'time':t})

	# Solve
	Tinout = np.zeros((modelPhy.nbctrlpts_total, len(time_inc)))
	problem_inc.solveFourierTransientProblem(Tinout=Tinout, Fext_list=Fext_list, time_list=time_inc, alpha=0.5)

	return problem_inc, time_inc, Tinout

quadArgs = {'quadrule': 'iga', 'type': 'leg'}
degree, cuts = 4, 6
problem_inc, time_inc, output = simulate(degree, cuts, quadArgs, cuts_time=7)
for k, i in enumerate(range(0, np.size(output, axis=1), 4)):
	problem_inc.part.postProcessingPrimal(fields={'temp':np.atleast_2d(output[:, i])}, 
									name='out_'+str(k), folder=folder)
	
run(folder=folder, filename='out_', nbFiles=k)