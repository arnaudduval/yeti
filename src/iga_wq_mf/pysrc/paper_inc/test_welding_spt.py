from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_stjob import stheatproblem
from pyevtk.vtk import VtkGroup

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/welding_spt/'
if not os.path.isdir(folder): os.mkdir(folder)

def run(folder=None):
	assert folder is not None, 'Folder unknown'
	print("Running group...")
	g = VtkGroup(folder)
	for i in range(33):
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
	POWER = 20; RADIUS = 0.25; VELOCITY = 0.5
	position = args['position']; time = args['time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); nc_tm = np.size(time); f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = time[i]
		if t <= 16: 
			rsquared = (x - VELOCITY*t)**2 + y**2
			f[:, i] = POWER*np.exp(-rsquared/(RADIUS**2))	
	return np.ravel(f, order='F')

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
	time_spt  = part1D(createUniformOpenCurve(degree, nbsteps, timespan), {'quadArgs': quadArgs})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	tmptable = np.zeros((3, 2)); tmptable[1, -1] = 1; tmptable[-1, 0] = 1
	sptnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], time_spt.nbctrlpts_total])
	boundary_spt = boundaryCondition(sptnbctrlpts)
	boundary_spt.add_DirichletConstTemperature(table=tmptable)

	# Transient model
	problem_spt  = stheatproblem(material, modelPhy, time_spt, boundary_spt)

	# Add external force
	Fext = problem_spt.compute_volForce(powerDensity, 
									{'position':problem_spt.part.qpPhy, 
									'time':problem_spt.time.qpPhy})
	
	# Solve
	Tinout = np.zeros(np.prod(sptnbctrlpts))
	problem_spt.solveFourierSTHeatProblem(Tinout=Tinout, Fext=Fext, isfull=False, isadaptive=False)

	return problem_spt, time_spt, Tinout

quadArgs = {'quadrule': 'iga', 'type': 'leg'}
degree, cuts = 4, 6
problem_spt, time_spt, temp_spt = simulate(degree, cuts, quadArgs, cuts_time=7)
output = np.reshape(temp_spt, order='F', 
		newshape=(problem_spt.part.nbctrlpts_total, time_spt.nbctrlpts_total),
		)
for k, i in enumerate(range(0, np.size(output, axis=1), 2)):
	problem_spt.part.postProcessingPrimal(fields={'temp':np.atleast_2d(output[:, i])}, 
									name='out_'+str(k), folder=folder)
	
run(folder)