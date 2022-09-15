from lib.__init__ import *
from lib.base_functions import (eval_basis_python,
                                iga_find_positions_weights,
                                create_knotvector, sigmoid
)
from lib.D1transientheat import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test/'
if not os.path.isdir(folder): os.mkdir(folder)

def conductivity(T):
    K = 0.1*np.ones(np.shape(T))
    return K

def capacity(T):
    # C = 1 + 0.2*np.exp(-0.1*abs(T))
    C = 5*np.ones(np.shape(T))
    return C

def source(qp):
    f = 0*np.sin(np.pi*qp)
    return f

# Define some properties to solver
alpha, JJ = 0.5, 1
properties = [JJ, conductivity, capacity, alpha]

# Create geometry
degree, nbel = 5, 32
nb_ctrlpts = degree + nbel
ctrlpts = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)

# Get basis and weights in IGA analysis
qp_cgg, weight_cgg = iga_find_positions_weights(degree, knotvector)
basis_cgg = eval_basis_python(degree, knotvector, qp_cgg)

# Define time discretisation
N = 100
time_list = np.linspace(0, 20, N)

# Compute volumetric heat source and external force
Fprop = source(qp_cgg)
FFend = compute_volsource_1D(JJ, basis_cgg, weight_cgg, Fprop)
FFend = np.atleast_2d(FFend).reshape(-1, 1)
Fext = np.kron(FFend, sigmoid(time_list))

# Define boundaries conditions
dod = [0, -1]
dof = np.arange(1, nb_ctrlpts-1, dtype=int)
temperature = np.zeros(np.shape(Fext))
temperature[0, :] = 0
temperature[-1,:] = sigmoid(time_list)

# Solve transient heat problem
solve_transient_heat_1D(properties, DB=basis_cgg, W=weight_cgg, Fext=Fext, 
                        time_list=time_list, dof=dof, dod=dod, Tinout=temperature)

# ------------------
# Post-treatement
# ------------------
# Create eval points
knots = np.linspace(0, 1, 101)
DB = eval_basis_python(degree, knotvector, knots)
temperature_interp = DB[0].T @ temperature

# Create mesh space-time
XX, TIME = np.meshgrid(knots*JJ, time_list)

# Plot figure
fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10,4))
for ax, variable in zip([ax1], [temperature_interp]):
    im = ax.contourf(XX, TIME, variable.T, 20)
    cbar = plt.colorbar(im)
    cbar.set_label('Temperature (K)')

    ax.grid(None)
    ax.set_ylabel('Time (s)')
    ax.set_xlabel('Position (m)')

fig.tight_layout()
fig.savefig(folder + 'Transient_heat_1D.png')

