"""
.. Test of plasticity 3D
.. We test how plasticity module works
.. Unities : MPa, mm2, kg
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.D3viscoplasticity import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test7/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 4, 4
isElastic = True

# Create geometry 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('CB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create physical model
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'density': 7.8e-6, 'young': 210e3, 'poisson': 0.3, 'sigmaY': 500, 'hardening':50e3, 'betahard':0.5}
modelPhy._set_material(material)

# Set Dirichlet boundaries
table_Dir = np.zeros((3, 2, 3), dtype=int)
table_Dir[0, 0, :] = 1
# table_Dir[0, 0, 0] = 1
# table_Dir[1, 0, 1] = 1
# table_Dir[2, 0, 2] = 1
Dirichlet = {'mechanical':table_Dir}
modelPhy._set_dirichlet_boundaries(Dirichlet)
Mdod = modelPhy._mechanical_dod

# Set Neumann boundaries
forces = [[0 for i in range(3)] for j in range(6)]
forces[1] = [0.0, 20.0, 0.0]
Neumann = {'mechanical': forces}
modelPhy._set_neumann_condition(Neumann)

# Set external forces
Fsurf = modelPhy.eval_force_surf()

if isElastic:
    # -------------
    # ELASTICITY
    # -------------
    # Compute iterative solution in fortran 
    displacement, resPCG = modelPhy.MFelasticity_fortran(indi=Mdod, Fext=Fsurf)

    # Plot residue
    resPCG = resPCG[resPCG>0]*100

    # Colors
    colorset = ['#377eb8', '#ff7f00', '#4daf4a',
            '#f781bf', '#a65628', '#984ea3',
            '#999999', '#e41a1c', '#dede00']

    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.semilogy(np.arange(len(resPCG)), resPCG)

    # Set properties
    ax.set_xlabel('Number of iterations')
    ax.set_ylabel('Relative residue (%)')
    fig.tight_layout()
    fig.savefig(folder + 'ElasticityRes.png')

    # ---------------------

    # # Interpolate displacement
    # modelPhy.export_results(u_ctrlpts=displacement, nbDOF=3, folder=folder)

    # Compute strain 
    strain = modelPhy.compute_strain(u=displacement)

    # Compute stress
    stress = modelPhy.compute_linear_stress(strain)
    stress_vm = np.zeros(np.shape(stress)[1])
    for k in range(np.shape(strain)[1]):
        stress_vm[k] = compute_stress_vonmises(3, stress[:, k])

    # Interpolate Von Mises field
    stress_ctrlpts = modelPhy.interpolate_ControlPoints(datafield=stress_vm)
    modelPhy.export_results(u_ctrlpts=stress_ctrlpts, nbDOF=1, folder=folder)

else:
    # --------------
    # PLASTICITY
    # --------------
    # Do ramp function (Fext increases linearly)
    nbStep = 6; dt = 1/nbStep
    Fext = np.zeros((*np.shape(Fsurf), nbStep+1))
    for i in range(1, nbStep+1): Fext[:, :, i] = i*dt*Fsurf

    # Solve system in fortran
    displacement, stress_vm = modelPhy.MFplasticity_fortran(Fext=Fext, indi=Mdod)

    # Interpolate displacement
    modelPhy.export_results(u_ctrlpts=displacement[:,:,-1], nbDOF=3, folder=folder)

    # # Interpolate Von Mises field
    # stress_ctrlpts = modelPhy.interpolate_ControlPoints(datafield=stress_vm[:, -1])
    # modelPhy.export_results(u_ctrlpts=stress_ctrlpts, nbDOF=1)