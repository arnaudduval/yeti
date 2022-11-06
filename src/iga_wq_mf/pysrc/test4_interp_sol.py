"""
.. Test of matrix free / interpolation
.. We test if matrix-Free fortran subroutines work correctly
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.fortran_mf_iga import fortran_mf_iga
from lib.physics import powden_rotring, temperature_rotring
from lib.base_functions import relativeError

# Set global variables
isIGA   = False
degree  = 6
cuts    = 4

if isIGA: cfortran = fortran_mf_iga
else: cfortran = fortran_mf_wq

# Create model 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('RQA', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# ----------------------
# By interpolation
# ----------------------
modelPhy = cfortran(modelIGA)
u_interp = modelPhy.interpolate_ControlPoints(funfield=temperature_rotring)
output = modelPhy.interpolate_field(u_ctrlpts=u_interp, nbDOF=1)
qp_sample, u_interp_sample = output[1], output[-1]

# Compare results
u_exact_sample = temperature_rotring(qp_sample)
error_1 = relativeError(u_interp_sample, u_exact_sample)
print("Error using interpolation : %.3e %%" %(error_1,))

# Add material 
material = {'capacity':1, 'conductivity':np.eye(3)}
modelPhy._set_material(material)

# Block boundaries
Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
modelPhy._set_dirichlet_boundaries(Dirichlet)
dof = modelPhy._thermal_dof
dod = modelPhy._thermal_dod 

# ----------------------
# By direct method
# ----------------------
ud = u_interp[dod]
K = modelPhy.eval_conductivity_matrix()
Knn = K[dof, :][:, dof]
Knd = K[dof, :][:, dod]
Fn = modelPhy.eval_source_vector(powden_rotring, dof) - Knd @ ud
un = sclin.solve(Knn.todense(), Fn)
usol = np.zeros(modelPhy._nb_ctrlpts_total)
usol[dof] = un; usol[dod] = ud

# Compare solutions 
u_interp_sample2 = modelPhy.interpolate_field(u_ctrlpts=usol, nbDOF=1)[-1]
error_2 = relativeError(u_interp_sample, u_interp_sample2)
print("Error interpolation/direct solution : %.3e %%" %(error_2,))

# ----------------------
# By iterative solver
# ----------------------
method_list = ["WP", "C", "TDS", "JMS", "TDC", "JMC"]
for method_name in method_list: 
    inputs = [Fn, 50, 1e-15, method_name]   
    un_t = modelPhy.MFsteadyHeat(*inputs)[0]
    usol_t = np.zeros(modelPhy._nb_ctrlpts_total)
    usol_t[dof] = un_t; usol_t[dod] = ud
    error_3 = relativeError(usol_t, usol)
    print("Error direct/iterative solution : %.3e %%" %(error_3,))