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
from lib.physics import setCprop, setKprop, powden
from lib.base_functions import relativeError

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test8/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
dataExist    = True
theta        = 1.0
degree, cuts = 5, 5
conductivity = 0.1
capacity     = 1.0

# Create model 
geometry = {'degree':[degree, degree, degree]}
modelGeo = geomdlModel('CB', **geometry)
modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                            np.array([cuts, cuts, cuts]))
modelPhy = fortran_mf_wq(modelIGA)

# Add material 
material = {'capacity':capacity, 'conductivity':conductivity*np.eye(3)}
modelPhy._set_material(material)

# Block boundaries
Dirichlet = {'thermal':np.array([[1, 1], [0, 0], [0, 0]])}
modelPhy._set_dirichlet_boundaries(Dirichlet)

# Add external force
N = 100
time_list = np.linspace(0, 20, N)
Fend  = modelPhy.eval_source_vector(powden)
Fendt = np.atleast_2d(Fend).reshape(-1, 1)
Fext  = np.kron(Fendt, sigmoid(time_list))

# Add constant temperature
modelPhy._add_thermal_IBC(np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)

if not dataExist:
    
    table_Kprop = create_table_properties(setKprop, prop=conductivity)
    table_Cprop = create_table_properties(setCprop, prop=capacity)

    GBound = np.zeros((len(modelPhy._thermal_dod), len(time_list)))
    for i in range(len(time_list)): GBound[:, i] = modelPhy._get_thermal_IBC()

    # Solve
    Tsol = modelPhy.MFtransientHeatNL(F=Fext, G=GBound, time_list=time_list,
                                    table_Kprop=table_Kprop, table_Cprop=table_Cprop, 
                                    methodPCG='JMC', theta=theta)[0]
    np.savetxt(folder + 'data3D.dat', Tsol)

else:
    # --------------
    # Post-treatment
    # --------------
    # Temperature compared to steady  
    Ku = modelPhy.eval_Ku(modelPhy._DirichletBound)[0]
    F  = Fext[:,-1] - Ku
    F  = F[modelPhy._thermal_dof]
    
    Tref  = modelPhy.MFsteadyHeat(b=F, methodPCG='JMC')[0]
    Tsol  = np.loadtxt(folder + 'data3D.dat')
    Tapp  = Tsol[:, -1][modelPhy._thermal_dof]
    error = relativeError(Tapp, Tref) 
    print('Relative error in steady-transient: %.5e %%' %error)

    # Temperature of mid-point
    samplesize = 61
    pos = int((samplesize-1)/2)
    Tpoint_list = []
    for i in range(np.shape(Tsol)[1]): 
        Tinterp = modelPhy.interpolate_field(u_ctrlpts=Tsol[:, i], nbDOF=1, samplesize=samplesize)[-1]
        Tpoint = Tinterp[pos + pos*samplesize + pos*samplesize**2]
        Tpoint_list.append(Tpoint)

    fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
    ax1.plot(time_list, Tpoint_list)
    ax1.set_ylim(top=0.6)

    # Get 1D data
    datapoint1D = np.loadtxt(folder + 'data1D.dat')
    ax2.semilogy(abs(Tpoint_list - datapoint1D[:, 1]), 'o', 
                nonpositive='mask')
    ax2.set_ylim([1e-12, 1e-3])

    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Temperature (K)')
    ax2.set_xlabel('Step')
    ax2.set_ylabel(r'$||T_{1D} - T_{3D}||$')
    fig.tight_layout()
    fig.savefig(folder + 'EvolTemp_midP31.png')
