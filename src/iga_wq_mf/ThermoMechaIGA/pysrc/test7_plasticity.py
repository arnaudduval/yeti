"""
.. Test of plasticity 3D
.. We test how plasticity module works
.. Unities : MPa, mm2, kg
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Set global variables
degree, cuts = 4, 4
isPlascticity = True

# Create geometry 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('CB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create physical model
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'density': 7.8e-6, 'young': 210e3, 'poisson': 0.3, 'sigmaY': 25e3, 'hardening':250, 'betahard':0.5}
modelPhy._set_material(material)

# Set Dirichlet boundaries
table_Dir = np.zeros((3, 2, 3), dtype=int)
table_Dir[0, 0, :] = 1
Dirichlet = {'mechanical':table_Dir}
modelPhy._set_dirichlet_boundaries(Dirichlet)
Mdod = modelPhy._mechanical_dod

# Set Neumann boundaries
forces = [[0 for i in range(3)] for j in range(6)]
forces[1] = [0.8, 0.9, 0.0]
Neumann = {'mechanical': forces}
modelPhy._set_neumann_condition(Neumann)

# Set external forces
Fsurf = modelPhy.eval_force_surf()

if not isPlascticity:
    # -------------
    # ELASTICITY
    # -------------
    # Compute iterative solution in python 
    itersol_py = modelPhy.MFelasticity_py(indi=Mdod, Fext=Fsurf)

    # Compute iterative solution in fortran 
    itersol_fr = modelPhy.MFelasticity_fortran(indi=Mdod, Fext=Fsurf)

    error = itersol_py - itersol_fr
    relerror = np.linalg.norm(error)/np.linalg.norm(itersol_py)
    print('Error between python and fortran' %relerror)

else:
    # --------------
    # PLASTICITY
    # --------------
    # Do ramp function (Fext increases linearly)
    nbStep = 6; dt = 1/nbStep
    Fext = np.zeros((*np.shape(Fsurf), nbStep+1))
    for i in range(1, nbStep+1): Fext[:, :, i] = i*dt*Fsurf

    # # Solve system in fortran
    # itersol_fr = modelPhy.MFplasticity_fortran(Fext=Fext, indi=Mdod)

    # Solve system in Python
    itersol_py = modelPhy.MFplasticity_py(Fext=Fext, indi=Mdod)
