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
from lib.D3viscoplasticity import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.base_functions import relativeError

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test8/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 2, 3
geoName = 'CB'

# Create model 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel(geoName, **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
											np.array([cuts, cuts, cuts]))
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'density': 7800, 'young': 1e9, 'poisson': 0.3, 'sigmaY': 500e9, 'hardening':50e9, 'betahard':0.5}
modelPhy._set_material(material)

# Set Dirichlet boundaries
table_Dir = np.zeros((3, 2, 3), dtype=int)
table_Dir[0, 0, 0] = 1
table_Dir[1, 0, 1] = 1
table_Dir[2, 0, 2] = 1
modelPhy._set_dirichlet_boundaries({'mechanical': table_Dir})
dod = modelPhy._mechanical_dod
dof = modelPhy._mechanical_dof
dof_gen = []
for i in range(3):
    dof_gen.extend(np.array(dof[i]) + i*modelPhy._nb_ctrlpts_total)

# Set Neumann boundaries
forces = [[0 for i in range(3)] for j in range(6)]
forces[1] = [1e8, 2e8, 0.0]
modelPhy._set_neumann_condition({'mechanical': forces})
Fsurf = modelPhy.eval_force_surf()

A = modelPhy.eval_stiffness_matrix()[dof_gen, :][:, dof_gen]
b = np.ndarray.flatten(Fsurf)[dof_gen]
disp_th = sclin.solve(A.todense(), b)

# -------------
# ELASTICITY
# -------------
# fig, ax = plt.subplots(nrows=1, ncols=1)

# # Solve in fortran 
# for methodPCG, label in zip(['WP', 'C', 'JMC'], ['w.o. preconditioner', 'Fast diag. (FD)', 'This work']):
#     displacement, resPCG = modelPhy.MFelasticity_fortran(indi=dod, Fext=Fsurf, nbIterPCG= 100, methodPCG=methodPCG)
#     resPCG = resPCG[resPCG>0]
#     ax.semilogy(np.arange(len(resPCG)), resPCG, label=label)

# ax.set_ybound(lower=1e-8, upper=10)
# ax.legend()
# ax.set_xlabel('Number of iterations of BiCGSTAB solver')
# ax.set_ylabel('Relative residue ' + r'$\displaystyle\frac{||r||_\infty}{||b||_\infty}$')
# fig.tight_layout()
# fig.savefig(folder + geoName + 'ElasRes.png')

displacement, resPCG = modelPhy.MFelasticity_fortran(indi=dod, Fext=Fsurf, nbIterPCG=100)
disp_app = np.ndarray.flatten(displacement)[dof_gen]
relerror = relativeError(disp_app, disp_th)
print(relerror)

# # strain = modelPhy.compute_strain(displacement)
# # strain_cp = []
# # for _ in range(6):
# # 	strain_cp.append(modelPhy.interpolate_ControlPoints(datafield=strain[_, :]))
# # strain_cp = np.array(strain_cp)
# # modelPhy.export_results(u_ctrlpts=strain_cp, nbDOF=6)