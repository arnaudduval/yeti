"""
.. Test of transient heat solver
.. ATTENTION: IT ONLY WORKS IN 'ISOTROPIC' MATERIALS
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.base_functions import create_table_properties, sigmoid

def setKpop(T, prop=1.0):
    y = prop*np.ones(len(T))
    return y

def setCpop(T, prop=1.0):
    y = prop*np.ones(len(T))
    return y

# Set global variables
degree, cuts = 4, 3
conductivity, capacity = 0.1, 1

# Set time simulation
N = 10
time_list = np.linspace(0, 10, N)

# Create material
table_Kprop = create_table_properties(setKpop, prop=conductivity)
table_Cprop = create_table_properties(setCpop, prop=capacity)

# Create geometry 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('CB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create model
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'capacity':capacity, 'conductivity':conductivity*np.eye(3)}
modelPhy._set_material(material)

# Block boundaries
Dirichlet = {'thermal':np.array([[1, 1], [0, 0], [0, 0]])}
modelPhy._set_dirichlet_boundaries(Dirichlet)

# Add constant temperature
modelPhy._add_thermal_IBC(np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)

# Create source
Fend = modelPhy.eval_source_vector(powdentest)

# ---------------------
# Transient model
# ---------------------
# Interpolate temperature on boundaries over time 
GBound = np.zeros((len(modelPhy._thermal_dod), len(time_list)))
for i in range(len(time_list)): GBound[:, i] = modelPhy._get_thermal_IBC()

# Compute 'velocity'
VBound = modelPhy._compute_velocity(GBound, time_list)

# Compute 'true' vector
Tsol1 = np.zeros((modelPhy._nb_ctrlpts_total, len(time_list)))
Tsol1[modelPhy._thermal_dod, :] = GBound
Kt0 = np.zeros(np.shape(Tsol1))
for i in range(len(time_list)):
    Kt0[:, i] = modelPhy.eval_Ku(Tsol1[:, i])

Vsol1 = np.zeros((modelPhy._nb_ctrlpts_total, len(time_list)))
Vsol1[modelPhy._thermal_dod, :] = VBound
Cv0 = np.zeros(np.shape(Vsol1))
for i in range(len(time_list)):
    Cv0[:, i] = modelPhy.eval_Cu(Vsol1[:, i])

# Add external force (transient)
Fendt = np.atleast_2d(Fend).reshape(-1, 1)
Fext  = np.kron(Fendt, sigmoid(time_list))
Fext  = Fext - Kt0 - Cv0
Fext  = Fext[modelPhy._thermal_dof, :]

# Solve transient problem at internal control points
T_ctrlpts = modelPhy.MFtransientHeatNL(F=Fext, time_list=time_list, table_Kprop=table_Kprop, table_Cprop=table_Cprop)
print('Well done')
# Tsol1[modelPhy._thermal_dof] = T_ctrlpts

#  modelPhy.export_results(u_ctrlpts=Tsol1[:, -1], nbDOF=1)

# # ---------------------
# # Steady model
# # ---------------------
# # Add material 
# material = {'capacity':1, 'conductivity':np.eye(3)}
# modelPhy._set_material(material)

# # Solve steady problem
# Fn = Fend[dof]
# Tn = modelPhy.MFsteadyHeat(Fn)[0]
# Tsol2 = np.zeros(modelPhy._nb_ctrlpts_total); Tsol2[dof] = Tn

# # Compare 
# Terror = Tsol1[:, -1] - Tsol2
# error = np.linalg.norm(Terror, np.inf)/np.linalg.norm(Tsol2, np.inf)*100
# print('Relative error : %.5f %%' %(error,))