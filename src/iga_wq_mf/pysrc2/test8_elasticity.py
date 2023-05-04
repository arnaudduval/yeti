"""
.. Test of elasticity 3D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : Pa (210e9)
..      - Length : m
..      - Force  : N
..      - Mass   : kg 
..      - Density: kg/m^3 (7.8e3)
..      - Gravity: m/s^2 (9.8)
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.lib_geomdl import Geomdl
from lib.lib_model import part
from lib.lib_material import mechamat
from lib.lib_step import step
from lib.lib_job import heatproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test8/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 2, 3
name = 'CB'

# Create model 
inputs = {'name': name, 'degree':degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
modelGeo = Geomdl(**inputs)
modelIGA = modelGeo.getIGAParametrization()
model    = part(modelIGA)

# Add material 
kwargs = {'density': 7800, 'elastic_modulus': 1e9, 'poisson_ratio': 0.3, 'elastic_limit': 500e9}
material    = mechamat(**kwargs)

# Set Dirichlet boundaries
boundary = step(model._nbctrlpts)
table = np.zeros((3, 2, 3), dtype=int)
table[0, 0, 0] = 1
table[1, 0, 1] = 1
table[2, 0, 2] = 1
boundary.add_DirichletDisplacement(table=table)

# # Set Neumann boundaries
# forces = [[0 for i in range(3)] for j in range(6)]
# forces[1] = [1e8, 2e8, 0.0]
# modelPhy._set_neumann_condition({'mechanical': forces})
# Fsurf = modelPhy.eval_force_surf()

# # -------------
# # ELASTICITY
# # -------------
# # fig, ax = plt.subplots(nrows=1, ncols=1)

# # # Solve in fortran 
# # for methodPCG, label in zip(['WP', 'C', 'JMC'], ['w.o. preconditioner', 'Fast diag. (FD)', 'This work']):
# #     displacement, resPCG = modelPhy.MFelasticity_fortran(indi=dod, Fext=Fsurf, nbIterPCG= 100, methodPCG=methodPCG)
# #     resPCG = resPCG[resPCG>0]
# #     ax.semilogy(np.arange(len(resPCG)), resPCG, label=label)

# # ax.set_ybound(lower=1e-8, upper=10)
# # ax.legend()
# # ax.set_xlabel('Number of iterations of BiCGSTAB solver')
# # ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
# # fig.tight_layout()
# # fig.savefig(folder + geoName + 'ElasRes.png')

# displacement, resPCG = modelPhy.MFelasticity_fortran(indi=dod, Fext=Fsurf, nbIterPCG=100)
# disp_app = np.ndarray.flatten(displacement)[dof_gen]
# relerror = relativeError(disp_app, disp_th)
# print(relerror)

# # # strain = modelPhy.compute_strain(displacement)
# # # strain_cp = []
# # # for _ in range(6):
# # # 	strain_cp.append(modelPhy.interpolate_ControlPoints(datafield=strain[_, :]))
# # # strain_cp = np.array(strain_cp)
# # # modelPhy.export_results(u_ctrlpts=strain_cp, nbDOF=6)