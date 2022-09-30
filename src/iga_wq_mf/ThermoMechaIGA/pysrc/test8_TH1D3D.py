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

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test8/'
if not os.path.isdir(folder): os.mkdir(folder)

def setKpop(T, prop=0.1):
    y = prop*np.ones(len(T))
    return y

def setCpop(T, prop=1.0):
    # y = 1 + prop*np.exp(-0.1*abs(T))
    y = prop*np.ones(len(T))
    return y

def powdentest(P:list):
    x = P[0, :]
    f = 0.0*np.sin(np.pi*x)
    return f

# Set global variables
degree, cuts = 4, 4
conductivity, capacity = 0.1, 1.0
newmark = 1.0

# Set time simulation
N = 100
time_list = np.linspace(0, 30, N)

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
Tsol, resPCG = modelPhy.MFtransientHeatNL(F=Fext, G=GBound, time_list=time_list,
                                table_Kprop=table_Kprop, table_Cprop=table_Cprop, 
                                methodPCG='JMS', newmark=newmark)

# modelPhy.export_results(u_ctrlpts=Tsol[:, -1], folder=folder, nbDOF=1)

# --------------
# Post-treatment
# --------------

# Temperature of mid-point
# -------------------------
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
datapoint1D = np.loadtxt(folder+'data1D.dat')
ax2.semilogy(abs(Tpoint_list - datapoint1D[:, 1]), 'o', 
            nonpositive='mask')

# Set properties
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Temperature (K)')
ax2.set_xlabel('Step')
ax2.set_ylabel(r'$||T_{1D} - T_{3D}||$')
fig.tight_layout()
fig.savefig(folder + 'EvolTemp_3D.png')

# Residue
# -------
resPCG = resPCG[:, resPCG[0, :]>0]

# Colors
colorset = ['#377eb8', '#ff7f00', '#4daf4a',
            '#f781bf', '#a65628', '#984ea3',
            '#999999', '#e41a1c', '#dede00']

fig, ax = plt.subplots(nrows=1, ncols=1)
for _ in range(np.shape(resPCG)[1]):
    step = resPCG[0, _]; iterNL = resPCG[1, _]
    newresidue = resPCG[2:, _]
    newresidue = newresidue[newresidue>0]*100
    ax.semilogy(np.arange(len(newresidue)), newresidue, 
                color=colorset[int(step%len(colorset))], alpha=1.0/iterNL)

# Set properties
ax.set_xlabel('Number of iterations')
ax.set_ylabel('Relative residue (%)')
fig.tight_layout()
fig.savefig(folder + 'TransientNL_CB.png')