"""
.. Test of transient heat solver substitution
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.base_functions import create_table_properties, sigmoid

def setKpop(T, prop=0.1):
    y = prop*np.ones(len(T))
    return y

def setCpop(T, prop=1.0):
    y = prop*np.ones(len(T))
    return y

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
degree, cuts = 4, 4
conductivity, capacity = 1.0, 1.0
isLinear = True

if isLinear: geoname = 'VB'
else: geoname = 'CB'

# Set time simulation
N = 10
time_list = np.linspace(0, 10, N)

# Create material
table_Kprop = create_table_properties(setKpop, prop=conductivity)
table_Cprop = create_table_properties(setCpop, prop=capacity)

# Create geometry 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel(geoname, **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))

# Create model
modelPhy = fortran_mf_wq(modelIGA)

if isLinear:

    # Add material 
    conductivity = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
    material = {'capacity':10.0, 'conductivity': conductivity}
    modelPhy._set_material(material)

    # Block boundaries
    Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
    modelPhy._set_dirichlet_boundaries(Dirichlet)

    # Create source
    Fext = modelPhy.eval_source_vector(powden_prism)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    for itmethod in ['FDC', 'JMC', 'JMS']:
        # Solve linear transient problem
        sol, residue = modelPhy.MF_THLSubs(F=Fext, method=itmethod)
        newresidue = residue[residue>0]*100
        ax.semilogy(np.arange(len(newresidue )), newresidue, label=itmethod)

    # Set properties
    ax.legend()
    ax.set_xlabel('Number of iterations')
    ax.set_ylabel('Relative error (%)')
    ax.set_ylabel('Relative residue (%)')
    fig.tight_layout()
    fig.savefig(folder + 'Transient3D.png')

else:
    # Block boundaries
    Dirichlet = {'thermal':np.array([[1, 1], [0, 0], [0, 0]])}
    modelPhy._set_dirichlet_boundaries(Dirichlet)
    modelPhy._add_thermal_IBC(np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)
    dod = modelPhy._thermal_dod
    dof = modelPhy._thermal_dof

    # Create source
    Fend = modelPhy.eval_source_vector(powdentest)

    # ---------------------
    # Transient model
    # ---------------------
    # Add boundaries temperature
    Gbound = np.zeros((len(dod), len(time_list)))
    for i in range(len(time_list)): Gbound[:, i] = modelPhy._get_thermal_IBC()

    # Add external force (transient)
    Fendt = np.atleast_2d(Fend).reshape(-1, 1)
    Fext  = np.kron(Fendt, sigmoid(time_list))
    del Fendt

    # Solve transient problem
    Tsol1 = modelPhy.MF_THNonLSubs(F=Fext, G=Gbound, time_list=time_list, 
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