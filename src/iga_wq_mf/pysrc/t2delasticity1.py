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
from lib.lib_material import (mechamat, array2symtensorForAll, evalTraceForAll, 
					computeVMStressForAll, symtensor2arrayForAll)
from lib.lib_load import forceSurf
from lib.lib_boundary import boundaryCondition
from lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t2delasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 6, 6
name = 'QA'

# Create model 
geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'Rin':1.e1, 'Rex':1.e2, 
			'XY':np.array([[0.0, 0.0], [1.e3, 0.0], [1.e3, 1.e3], [0.0, 1.e3]])}}
quadArgs  = {'quadrule': 'wq', 'type': 1}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
model    = part(modelIGA, quadArgs=quadArgs)

# Add material 
matArgs  = {'elastic_modulus':2e5, 'elastic_limit':100, 'poisson_ratio': 0.3}
material = mechamat(matArgs)

# Set Dirichlet boundaries
boundary = boundaryCondition(model.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
if name == 'SQ':
	table[0, 0, 0] = 1
	table[0, 0, 1] = 1
elif name == 'QA':
	# table[1, 0, 0] = 1
	table[1, 1, 0] = 1
	table[1, 0, 1] = 1
else: raise Warning('Not possible')
boundary.add_DirichletDisplacement(table=table)

# Elasticity problem
problem = mechaproblem(material, model, boundary)
Fext = problem.eval_surfForce(forceSurf, nbFacePosition=1)

# -------------
# ELASTICITY
# -------------
displacement = problem.solveElasticityProblemFT(Fext=Fext)[0]

strain_qp  = problem.compute_strain(displacement)
Tstrain_qp = array2symtensorForAll(strain_qp, 2)
traceStrain_qp = evalTraceForAll(strain_qp, 2)
devStrain_qp   = Tstrain_qp
for i in range(2): devStrain_qp[i, i, :] -= 1.0/3.0*traceStrain_qp

Tstress_qp = 2*problem.material.lame_mu*devStrain_qp
stress_qp  = symtensor2arrayForAll(Tstress_qp, 2)
stress_vm_qp = computeVMStressForAll(stress_qp, 2)
print(np.max(stress_vm_qp))
stress_vm_cp = problem.L2projectionCtrlpts(datafield=stress_vm_qp)
model.exportResultsCP(fields={'disp':displacement, 'svm':stress_vm_cp}, folder=folder)

# disp_interp = problem.part.interpolateMeshgridField(u_ctrlpts=displacement, sampleSize=2500)[-1]
# disp_norm = np.sqrt(disp_interp[0, :]**2+disp_interp[1, :]**2)
# disp_interp = np.vstack([disp_interp, disp_norm])
# np.save(folder+'disp_interp_ref.npy', disp_interp)

# fig, ax = plt.subplots()
# # Solve in fortran 
# for i, [methodPCG, label] in enumerate(zip(['WP', 'C', 'JMC'], 
# 							['w.o. preconditioner', 'Fast diag. (FD)', 'This work'])):
#     problem.addSolverConstraints(solverArgs={'PCGmethod': methodPCG})
#     displacement, resPCG = problem.solveElasticityProblemFT(Fext=Fext)
#     resPCG = resPCG[resPCG>0]
#     ax.semilogy(np.arange(len(resPCG)), resPCG, '-', label=label, marker=markerSet[i])

# ax.set_ybound(lower=1e-12, upper=1e1)
# ax.legend()
# ax.set_xlabel('Number of iterations of BiCGSTAB solver')
# ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
# fig.tight_layout()
# fig.savefig(folder + name + 'ElasRes.png')