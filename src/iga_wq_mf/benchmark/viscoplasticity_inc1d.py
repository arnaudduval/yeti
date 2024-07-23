from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import mechaproblem1D
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition

# Global variables
NBSTEPS = 5
TIME_LIST = np.linspace(0, 1, NBSTEPS)
YOUNG, LENGTH  = 2e11, 1

def forceVol(P:dict):
	position = P['position']
	force = 2.5e6*np.ones(len(position))
	return force

# Define geometry
degree, nbel = 1, 4
quadArgs = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
geometry = createUniformOpenCurve(degree, nbel, LENGTH)
modelPhy = part1D(geometry, quadArgs)

# Define plasticity model
MECHAMATERIAL = mechamat({'elastic_modulus':YOUNG, 
			'elastic_limit':1e6, 
			'poisson_ratio':0.3,
			'viscoparameter':1e0,
			'isoHardLaw': {'name':'linear', 'Eiso':1e6}})

# Set boundary conditions
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.array([[1, 0]]))

# Define problem
problem = mechaproblem1D(mechanical_material=MECHAMATERIAL, part=modelPhy, boundary=boundary)

# Add external force
Fref = np.atleast_2d(problem.compute_volForce(forceVol)).transpose()
Fext_list = np.kron(Fref, TIME_LIST)
# for i, t in enumerate(TIME_LIST): Fext_list[-1, i] += 1e7*t 

# Solve problem
displacement = np.zeros(np.shape(Fext_list))
# blockPrint()
# Allstrain, Allstress, Allplseq, AllCep = problem.solvePlasticityProblem(displacement, Fext_list)
# enablePrint()
# print(Allstress[:, -1])

# blockPrint()
Allstrain, Allstress, Allplseq, AllCep = problem.solveViscoPlasticityProblem(displacement, Fext_list, TIME_LIST)
enablePrint()
print(Allstress[:, -1])
print(Allstrain[:, -1])
print(Allplseq[:, -1])