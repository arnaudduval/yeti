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
degree, cuts = 4, 5
name = 'QA'

# Create model 
geoArgs = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
			'extra':{'Rin':1.0, 'Rex':4.0}
}
quadArgs  = {'quadrule': 'wq', 'type': 1}

modelGeo = Geomdl(geoArgs)
modelIGA = modelGeo.getIGAParametrization()
model    = part(modelIGA, quadArgs=quadArgs)

# Add material 
matArgs  = {'elastic_modulus':1e3, 'elastic_limit':1e10, 'poisson_ratio': 0.3}
material = mechamat(matArgs)

# Set Dirichlet boundaries
boundary = boundaryCondition(model.nbctrlpts)
table = np.zeros((2, 2, 2), dtype=int)
table[1, 1, 0] = 1
table[1, 0, 1] = 1
boundary.add_DirichletDisplacement(table=table)

# Elasticity problem
problem = mechaproblem(material, model, boundary)
Fext = problem.eval_surfForce(forceSurf, nbFacePosition=1)

# -------------
# ELASTICITY
# -------------
displacement, _, stress_qp = problem.solveElasticityProblemFT(Fext=Fext)
stress_cp = problem.L2projectionCtrlpts(datafield=stress_qp)
model.exportResultsCP(fields={'disp':displacement, 'S':stress_cp}, folder=folder)

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