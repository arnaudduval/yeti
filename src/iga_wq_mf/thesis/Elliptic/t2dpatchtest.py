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
from pysrc.lib.lib_job3d import mechaproblem

# Set global variables
degree, cuts = 2, 1

# Create model 
geoArgs = {
			'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'XY':np.array([[0.0, 0.0], [1., 0.0], [1., 1.], [0.0, 1.]])}
		}
quadArgs  = {'quadrule': 'iga'}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
modelPhy = part(modelIGA, quadArgs=quadArgs)

# Add material 
material = mechamat({'elastic_modulus':1e2, 'elastic_limit':100, 'poisson_ratio': 0.3})

# Set Dirichlet boundaries
boundary = boundaryCondition(modelPhy.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
table[0, 0, 0] = 1; table[1, 0, 1] = 1
boundary.add_DirichletDisplacement(table=table)

# Elasticity problem
problem = mechaproblem(material, modelPhy, boundary)

# TRACTION FOLLOWING X
# --------------------
def forceSurfFun(P:list):
	ref  = np.array([4e1, 0.0])
	prop = np.zeros((2, np.size(P, axis=1)))
	for i in range(2): prop[i, :] = ref[i] 
	return prop
Fext = problem.compute_surfForce(forceSurfFun, nbFacePosition=1)[0]

# # TRACTION FOLLOWING Y
# # --------------------
# def forceSurfFun(P:list):
# 	ref  = np.array([0.0, 1e1])
# 	prop = np.zeros((2, np.size(P, axis=1)))
# 	for i in range(2): prop[i, :] = ref[i] 
# 	return prop
# Fext = problem.compute_surfForce(forceSurfFun, nbFacePosition=3)[0]

# Solve system
start = time.process_time()
disp_cp, residual = problem._solveLinearizedElasticityProblem(Fext=Fext)
print(len(residual[residual>0]))
stop = time.process_time()
print('CPU time: %5e' %(stop-start))

# Post processing
straintmp = problem.interpolate_strain(disp_cp)
strain_qp = np.zeros((6, problem.part.nbqp_total))
strain_qp[0:2, :] = straintmp[0:2, :]; strain_qp[3, :] = straintmp[-1, :]
stress_qp = problem.mechamaterial.evalElasticStress(strain_qp)
VMstress_qp = computeSymTensorNorm4All(stress_qp)
print('Von misses max:%.4e, min:%.4e' %(VMstress_qp.max(), VMstress_qp.min()))
print('Difference: %.4e' %(abs(VMstress_qp.max()-VMstress_qp.min())))
