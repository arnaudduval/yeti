from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import mechaproblem1D
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

# Define geometry
quadArgs = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
geometry = createUniformOpenCurve(DEGREE, NBEL, LENGTH)
modelPhy = part1D(geometry, quadArgs)

# Set boundary conditions
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.array([[1, 0]]))

# Define problem
problem = mechaproblem1D(mechanical_material=MECHAMATERIAL, part=modelPhy, boundary=boundary)

# Add external force
TIME_LIST = np.linspace(0, 1, NBSTEPS)
Fref = np.atleast_2d(problem.compute_volForce(forceVol_inc)).transpose()
Fext_list = np.kron(Fref, TIME_LIST)
# for i, t in enumerate(TIME_LIST): Fext_list[-1, i] += 1e7*t 

# Solve problem
displacement = np.zeros(np.shape(Fext_list))
blockPrint()
Allstrain, Allstress, Allplseq, AllCep = problem.solvePlasticityProblem(displacement, Fext_list)
enablePrint()
print(Allstress[:, -1])

