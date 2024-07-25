
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import stmechaproblem1D
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition

# Global variables
DEGREE, NBEL = 2, 21
NBSTEPS = 21
YOUNG, LENGTH = 1e6, 1

# Define plasticity model
MECHAMATERIAL = mechamat({'elastic_modulus':YOUNG, 
			'elastic_limit':1.e1, 
			'poisson_ratio':0.3,
			'isoHardLaw': {'name':'linear', 'Eiso':1e3}})

def forceVol_inc(args:dict):
	position = args['position']; time = args.get('time', 1)
	force = 1e2*np.ones(len(position))*time
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
quadArgs = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
geometry = createUniformOpenCurve(DEGREE, NBEL, LENGTH)
modelPhy = part1D(geometry, quadArgs)
timespan = part1D(createUniformOpenCurve(2, NBSTEPS, 1.), kwargs={'quadArgs': {'quadrule': 'iga', 'type': 'leg'}})

# Set boundary conditions
boundary = boundaryCondition(np.array([modelPhy.nbctrlpts_total, timespan.nbctrlpts_total, 1]))
boundary.add_DirichletConstTemperature(table=np.array([[1, 1], [1, 0]], dtype=int))

# Define problem
problem = stmechaproblem1D(mechanical_material=MECHAMATERIAL, part=modelPhy, tspan=timespan, boundary=boundary)

# Add external force
Fext_list = np.ravel(problem.compute_volForce(forceVol_spt))
displacement = np.zeros(np.shape(Fext_list))
disp_cp = problem.solveViscoPlasticityProblem(displacement, Fext_list, init=1e-5)

print(np.max(disp_cp))