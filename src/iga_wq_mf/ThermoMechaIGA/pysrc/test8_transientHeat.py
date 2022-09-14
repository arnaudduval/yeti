"""
.. Test of transient heat solver
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.base_functions import create_table_properties

def setKpop(T):
    y = 0.1*np.ones(len(T))
    return y

def setCpop(T):
    y = 1.0*np.ones(len(T))
    return y

# Set global variables
degree, cuts = 4, 3

# Set time simulation
N, n = 10, 0
tt = np.linspace(0, 10, N)
ones = np.ones(len(tt))
time_list = np.zeros(N+n)
time_list[:N] = tt
time_list[N:] = [tt[-1] + 1*(i+1) for i in range(n)]

# Create material
table_Kprop = create_table_properties(setKpop)
table_Cprop = create_table_properties(setCpop, uref=np.linspace(-20, 20, 21))

# Create geometry 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('CB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create model
modelPhy = fortran_mf_wq(modelIGA)

# Block boundaries
Dirichlet = {'thermal':np.array([[1, 1], [0, 0], [0, 0]])}
modelPhy._set_dirichlet_boundaries(Dirichlet)
modelPhy.add_thermal_initial_dirichlet_condition(np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)
dod = modelPhy._thermal_dod
dof = modelPhy._thermal_dof

# Create source
Fend = modelPhy.eval_source_vector(powdentest)

# ---------------------
# Transient model
# ---------------------
# Add boundaries temperature
Gbound = np.zeros((len(dod), len(time_list)))
for i in range(len(time_list)):
    Gbound[:, i] = modelPhy.get_thermal_initial_dirichlet_condition()

# Add external force (transient)
Fendt = np.atleast_2d(Fend).reshape(-1, 1)
Fext = np.zeros((len(Fendt), len(time_list)))
Fext[:,:len(tt)] = np.kron(Fendt, tt)/np.max(tt)
for i in range(len(tt), len(time_list)): Fext[:,i] = Fendt[:, 0]
del Fendt

# Solve transient problem
Tsol1 = modelPhy.MFtransientHeatNL(F=Fext, G=Gbound, time_list=time_list, 
                                table_Kprop=table_Kprop, table_Cprop=table_Cprop)
modelPhy.export_results(u_ctrlpts=Tsol1[:, -1], nbDOF=1)

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