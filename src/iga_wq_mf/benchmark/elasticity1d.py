"""
.. Test of elastoplasticity 1D
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import mechaproblem1D
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition

# Global variables
NBSTEPS = 201
TIME_LIST = np.linspace(0, np.pi, NBSTEPS)
YOUNG, CST, LENGTH  = 2e11, 4.e7, 1
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':1e6, 'poisson_ratio':0.3,
		'isoHardLaw': {'name':'linear', 'Eiso':YOUNG/10}}
MECHAMATERIAL = mechamat(MATARGS)

def forceVol(P:dict):
	x = P['position']
	force = CST*(x - 1/10*x**2)
	return force

degree, nbel = 2, 128
quadArgs = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
geometry = createUniformOpenCurve(degree, nbel, LENGTH)
modelPhy = part1D(geometry, quadArgs)
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.array([[1, 1]]))
problem = mechaproblem1D(mechanical_material=MECHAMATERIAL, part=modelPhy, boundary=boundary)
model2return = deepcopy(problem)
Fref = np.atleast_2d(problem.compute_volForce(forceVol)).transpose()
Fext_list = np.kron(Fref, np.sin(TIME_LIST))
displacement = np.zeros(np.shape(Fext_list))
problem.solvePlasticityProblem(displacement, Fext_list)
