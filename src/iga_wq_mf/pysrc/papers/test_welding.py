from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import heatproblem
from pysrc.lib.lib_stjob import stheatproblem
from pyevtk.vtk import VtkGroup

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/welding/'
if not os.path.isdir(folder): os.mkdir(folder)

def run(folder=None):
	assert folder is not None, 'Folder unknown'
	print("Running group...")
	g = VtkGroup(folder)
	for i in range(65):
		g.addFile(filepath = folder + "out_"+str(i)+".vts", sim_time = i)
	g.save()

ISLINEAR = True

def conductivityProperty(args:dict):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.0],[0.0, 1.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			if ISLINEAR: Kprop[i, j, :] = Kref[i, j]
			else: Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
	return Kprop 

def capacityProperty(args:dict):
	temperature = args['temperature']
	if ISLINEAR: Cprop = np.ones(len(temperature))
	else: Cprop = (1.0 + np.exp(-np.abs(temperature)))
	return Cprop

def powerDensity(args:dict):
	POWER = 20; RADIUS = 0.5; VELOCITY = 0.5
	position = args['position']; t = args['time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); f = np.zeros(nc_sp)
	if t < 2 or t > 18: return f
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
	time_crv  = part1D(createUniformCurve(1, nbsteps, timespan), {'quadArgs': quadArgs})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	dirichlet_table = np.zeros((2, 2)); dirichlet_table[1, -1] = 1
	boundary_inc = boundaryCondition(modelPhy.nbctrlpts)
	boundary_inc.add_DirichletConstTemperature(table=dirichlet_table)

	dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
	stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], time_crv.nbctrlpts_total])
	boundary_st = boundaryCondition(stnbctrlpts)
	boundary_st.add_DirichletConstTemperature(table=dirichlet_table)

	# Transient model
	problem_inc = heatproblem(material, modelPhy, boundary_inc)
	problem_st  = stheatproblem(material, modelPhy, time_crv, boundary_st)

	# Add external force 
	Fext_list = np.zeros((problem_inc.part.nbctrlpts_total, len(time_inc)))
	for i, t in enumerate(time_inc):
		Fext_list[:, i] = problem_inc.compute_volForce(powerDensity, args={'position':problem_inc.part.qpPhy, 'time':t})

	# Solve
	Tinout = np.zeros((modelPhy.nbctrlpts_total, len(time_inc)))
	problem_inc.solveFourierTransientProblem(Tinout=Tinout, Fext_list=Fext_list, time_list=time_inc, alpha=0.5)

	return problem_inc, problem_st, Tinout

quadArgs = {'quadrule': 'iga', 'type': 'leg'}
degree, cuts = 4, 5
problem_inc, problem_st, output = simulate(degree, cuts, quadArgs, cuts_time=6)
for i in range(np.size(output, axis=1)):
	problem_inc.part.exportResultsCP(fields={'temp':np.atleast_2d(output[:, i])}, 
									name='out_'+str(i), folder=folder)
	
run(folder)