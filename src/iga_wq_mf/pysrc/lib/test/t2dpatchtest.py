"""
.. Test of elasticity 2D
.. A square plate is under uniaxial traction (following x or y). 
.. It has been imposed symetry conditions, so the strains has to be the same everywhere
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat, computeSymTensorNorm4All
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Set global variables
sampleSize   = 2500
degree, cuts = 1, 8

# Create model 
geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'XY':np.array([[0.0, 0.0], [1.e3, 0.0], [1.e3, 1.e3], [0.0, 1.e3]])}}
quadArgs  = {'quadrule': 'wq'}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
model    = part(modelIGA, quadArgs=quadArgs)

# Add material 
matArgs  = {'elastic_modulus':2e5, 'elastic_limit':100, 'poisson_ratio': 0.3}
material = mechamat(matArgs)

# Set Dirichlet boundaries
boundary = boundaryCondition(model.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
table[0, 0, 0] = 1
table[1, 0, 1] = 1
boundary.add_DirichletDisplacement(table=table)

# Elasticity problem
problem = mechaproblem(material, model, boundary)

# # TRACTION FOLLOWING X
# # --------------------
# def forceSurfFun(P:list):
# 	ref  = np.array([4e1, 0.0])
# 	prop = np.zeros((2, np.size(P, axis=1)))
# 	for i in range(2): prop[i, :] = ref[i] 
# 	return prop
# Fext = problem.eval_surfForce(forceSurfFun, nbFacePosition=1)

# TRACTION FOLLOWING Y
# --------------------
def forceSurfFun(P:list):
	ref  = np.array([0.0, 4e1])
	prop = np.zeros((2, np.size(P, axis=1)))
	for i in range(2): prop[i, :] = ref[i] 
	return prop
Fext = problem.compute_surfForce(forceSurfFun, nbFacePosition=3)[0]

# Solve in fortran 
start = time.time()
disp_cp = problem.solveElasticityProblem(Fext=Fext)[0]
stop = time.time()
print('CPU time: %5e' %(stop-start))

strain = problem.interpolate_strain(disp_cp)
stress = problem.mechamaterial.evalElasticStress(strain, problem.part.dim)
stress_vm = computeSymTensorNorm4All(stress, problem.part.dim)
print('Von misses max:%.4e, min:%.4e' %(stress_vm.max(), stress_vm.min()))
print('Difference: %.4e' %(abs(stress_vm.max()-stress_vm.min())))
