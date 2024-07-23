
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import stmechaproblem1D
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition

# Global variables
NBSTEPS = 11
YOUNG, LENGTH  = 1e3, 1

def forceVol_inc(args:dict):
	position = args['position']; time = args.get('time', 1)
	force = 2.5*np.ones(len(position))*time
	return force

def forceVol_spt(args:dict):
	time = args['time']
	position = args['position']
	nc_sp = len(position); nc_tm = np.size(time)
	f = np.zeros((nc_sp, nc_tm))
	for i in range(nc_tm):
		t = time[i]
		f[:, i] = forceVol_inc(args={'time':t, 'position':position})
	return np.ravel(f, order='F')

# Define geometry
degree, nbel = 1, 10
quadArgs = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
geometry = createUniformOpenCurve(degree, nbel, LENGTH)
modelPhy = part1D(geometry, quadArgs)
timespan = part1D(createUniformOpenCurve(1, NBSTEPS-1, 1.), kwargs={'quadArgs': {'quadrule': 'iga', 'type': 'leg'}})

# Define plasticity model
MECHAMATERIAL = mechamat({'elastic_modulus':YOUNG, 
			'elastic_limit':1, 
			'poisson_ratio':0.3,
			'viscoparameter':1e-2,
			'isoHardLaw': {'name':'linear', 'Eiso':1e2}})

# Set boundary conditions
boundary = boundaryCondition(np.array([modelPhy.nbctrlpts_total, timespan.nbctrlpts_total, 1]))
boundary.add_DirichletConstTemperature(table=np.array([[1, 0], [1, 0]], dtype=int))

# Define problem
problem = stmechaproblem1D(mechanical_material=MECHAMATERIAL, part=modelPhy, tspan=timespan, boundary=boundary)

# Add external force
Fext_list = np.ravel(problem.compute_volForce(forceVol_spt))
displacement = np.zeros(np.shape(Fext_list))
problem.solveViscoPlasticityProblem(displacement, Fext_list)
