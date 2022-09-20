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
N = 20
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

# ---------------------
# Transient model
# ---------------------
# Interpolate temperature on boundaries over time 
GBound = np.zeros((len(modelPhy._thermal_dod), len(time_list)))
for i in range(len(time_list)): GBound[:, i] = modelPhy._get_thermal_IBC()

# Add external force (transient)
Fend = modelPhy.eval_source_vector(powdentest)
Fendt = np.atleast_2d(Fend).reshape(-1, 1)
Fext  = np.kron(Fendt, sigmoid(time_list))

# Solve transient problem at internal control points
Tsol = modelPhy.MFtransientHeatNL(F=Fext, G=GBound, time_list=time_list,
                                table_Kprop=table_Kprop, table_Cprop=table_Cprop)
