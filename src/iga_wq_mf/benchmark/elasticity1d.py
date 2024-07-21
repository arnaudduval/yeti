from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import mechaproblem1D
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition

# Global variables
NBSTEPS = 201
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
YOUNG, LENGTH  = 2e11, 1

def forceVol(P:dict):
	x = P['position']
	force = 1e7*(x - 1/10*x**2)
	return force

# Define geometry
degree, nbel = 2, 128
quadArgs = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
geometry = createUniformOpenCurve(degree, nbel, LENGTH)
modelPhy = part1D(geometry, quadArgs)

# Define plasticity model
MECHAMATERIAL = mechamat({'elastic_modulus':YOUNG, 
			'elastic_limit':1e6, 
			'poisson_ratio':0.3,
			'isoHardLaw': {'name':'linear', 'Eiso':YOUNG/10}})

# Set boundary conditions
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.array([[1, 1]]))

# Define problem
problem = mechaproblem1D(mechanical_material=MECHAMATERIAL, part=modelPhy, boundary=boundary)

# Add external force
Fref = np.atleast_2d(problem.compute_volForce(forceVol)).transpose()
Fext_list = np.kron(Fref, np.sin(TIME_LIST))

# Solve problem
displacement = np.zeros(np.shape(Fext_list))
problem.solvePlasticityProblem(displacement, Fext_list)
