"""
.. Test of mecanical displacement 1D
.. Author: Fabio MADIE
.. Joaquin Cornejo added some corrections 28 nov. 2024
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import mechaproblem1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import mechamat

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

YOUNG, RHO, LENGTH = 0.5e9, 500, 1
MECHAMATERIAL = mechamat({'elastic_modulus':YOUNG, 
						'elastic_limit':1e6, 
						'poisson_ratio':0.3,
						'isoHardLaw': {'name':'None'}})
MECHAMATERIAL.addDensity(RHO, isIsotropic=True)
WAVEVEL = np.sqrt(YOUNG/RHO)
COEF = 100 # (20*np.pi/WAVEVEL)**2

def computeIitialDisplacement(args:dict):
	x = args['position']
	u_ini = np.exp(-COEF/2*(x-LENGTH/2)**2)
	return u_ini

def computeIitialVelocity(args:dict):
	x = args['position']
	v_ini = COEF*WAVEVEL*(x-LENGTH/2)*np.exp(-COEF/2*(x-LENGTH/2)**2)
	return v_ini

# Create geometry
degree, nbel = 2, 32
crv = createUniformOpenCurve(degree, nbel, LENGTH)
modelIGA = part1D(crv, kwargs={'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}) #iga : quadrature de Gauss # wq quadrature pondérée

# Create boundary condition
boundary = boundaryCondition(modelIGA.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.array([[1 , 0]]))

# Set mechanical problem
problem = mechaproblem1D(MECHAMATERIAL, modelIGA, boundary)
eigs = np.sqrt(problem.compute_EigenvalueProblem()[0])
Nmin = int(np.ceil(LENGTH/WAVEVEL*eigs.max()/2))

# Create external force	
timeList = np.linspace(0, LENGTH/WAVEVEL, 1000)
FextList = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
displacement = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
velocity = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
acceleration = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))

# Compute initial values
displacement[:, 0] = problem.L2projectionCtrlpts(computeIitialDisplacement({'position': problem.part.qpPhy}))
velocity[:,0] = problem.L2projectionCtrlpts(computeIitialVelocity({'position': problem.part.qpPhy}))
problem.solveNewmarkresolution(displacement, velocity, acceleration, FextList, timeList)

fig, ax = plt.subplots()
disp_app, pos = problem.interpolateMeshgridField(displacement[:, 0])
ax.plot(pos, disp_app, '--', label='Init')
disp_app, pos = problem.interpolateMeshgridField(displacement[:, -1])
ax.plot(pos, disp_app, '.', label='Final')
disp_ex = computeIitialDisplacement({'position': pos})
ax.plot(pos, disp_ex, label='Exact', alpha=0.5)
ax.legend()
fig.savefig(folder+'displacement')
