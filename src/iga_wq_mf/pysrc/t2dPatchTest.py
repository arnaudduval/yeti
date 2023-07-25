"""
.. Test of elastoplasticity 2D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa (200e3)
..      - Length : mm
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.lib_geomdl import Geomdl
from lib.lib_part import part
from lib.lib_material import (mechamat, array2symtensorForAll, evalMultiTraceVgt, 
					computeMultiVMStressVgt, symtensor2arrayForAll)
from lib.lib_boundary import boundaryCondition
from lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t2delasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
sampleSize   = 2500
degree, cuts = 2, 8
name = 'SQ'

# Create model 
geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'Rin':5.e2, 'Rex':1.e3, 
			'XY':np.array([[0.0, 0.0], [1.e3, 0.0], [1.e3, 1.e3], [0.0, 1.e3]])}}
quadArgs  = {'quadrule': 'iga', 'type': 'leg'}

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
Fext = problem.eval_surfForce(forceSurfFun, nbFacePosition=3)

# -------------
# ELASTICITY
# -------------
# Solve in fortran 
start = time.time()
disp_cp = problem.solveElasticityProblemFT(Fext=Fext)[0]
stop = time.time()
print('CPU time: %5e' %(stop-start))
# model.exportResultsCP(u_ctrlpts=disp_cp, nbDOF=2, folder=folder)

strain  = problem.compute_strain(disp_cp)
TenStrain = array2symtensorForAll(strain, 2)
trStrain  = evalMultiTraceVgt(strain, 2)
devStrain = TenStrain
for i in range(2): devStrain[i, i, :] -= 1.0/3.0*trStrain
Tstress = 2*problem.material.lame_mu*devStrain
stress    = symtensor2arrayForAll(Tstress, 2)
stress_vm = computeMultiVMStressVgt(stress, 2)
print('Von misses max:%.4e, min:%.4e' %(stress_vm.max(), stress_vm.min()))
print('Difference: %.4e' %(abs(stress_vm.max()-stress_vm.min())))
