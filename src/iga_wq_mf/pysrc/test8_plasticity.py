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

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test8/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 4, 4

# Create model 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('CB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'density': 7800, 'young': 210e9, 'poisson': 0.3, 'sigmaY': 500e6, 'hardening':50e9, 'betahard':0.5}
modelPhy._set_material(material)

# Set Dirichlet boundaries
table_Dir = np.zeros((3, 2, 3), dtype=int)
table_Dir[0, 0, 0] = 1
table_Dir[1, 0, 1] = 1
table_Dir[2, 0, 2] = 1
modelPhy._set_dirichlet_boundaries({'mechanical': table_Dir})
dod = modelPhy._mechanical_dod

# Set Neumann boundaries
forces = [[0 for i in range(3)] for j in range(6)]
forces[1] = [1e8, 2e8, 0.0]
modelPhy._set_neumann_condition({'mechanical': forces})
Fsurf = modelPhy.eval_force_surf()

# --------------
# PLASTICITY
# --------------
nbStep = 6; dt = 1/nbStep
Fext = np.zeros((*np.shape(Fsurf), nbStep+1))
for i in range(1, nbStep+1): Fext[:, :, i] = i*dt*Fsurf

# Solve in fortran
displacement, stress_vm = modelPhy.MFplasticity_fortran(Fext=Fext, indi=dod)

# # Interpolate displacement
# modelPhy.export_results(u_ctrlpts=displacement[:,:,-1], nbDOF=3, folder=folder)

# # Interpolate Von Mises field
# stress_ctrlpts = modelPhy.interpolate_ControlPoints(datafield=stress_vm[:, -1])
# modelPhy.export_results(u_ctrlpts=stress_ctrlpts, nbDOF=1)