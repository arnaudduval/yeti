from lib.__init__ import *
from lib.base_functions import (eval_basis_python,
                                iga_find_positions_weights,
                                create_knotvector, sigmoid
)
from lib.D1transientheat import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test8/'
if not os.path.isdir(folder): os.mkdir(folder)

def setKprop(T):
    y = 0.1 + 0.1*np.exp(-0.01*abs(T))
    return y

def setCprop(T):
    y = 1.0 + 1.0*np.exp(-0.01*abs(T))
    # y = 1.0*np.ones(np.shape(T))
    return y

def powdentest(qp):
    f = 0*np.sin(np.pi*qp)
    return f

# Define some properties to solver
newmark, JJ = 1.0, 1.0
properties = [JJ, setKprop, setCprop, newmark]

# Create geometry
degree, nbel = 5, 64
nb_ctrlpts = degree + nbel
ctrlpts = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)

# Get basis and weights in IGA analysis
qp_cgg, weight_cgg = iga_find_positions_weights(degree, knotvector)
basis_cgg = eval_basis_python(degree, knotvector, qp_cgg)

# Define time discretisation
N = 100
time_list = np.linspace(0, 30, N)

# Compute volumetric heat source and external force
Fprop = powdentest(qp_cgg)
FFend = compute_volsource_1D(JJ, basis_cgg, weight_cgg, Fprop)
FFend = np.atleast_2d(FFend).reshape(-1, 1)
Fext = np.kron(FFend, sigmoid(time_list))

# Define boundaries conditions
dod = [0, -1]
dof = np.arange(1, nb_ctrlpts-1, dtype=int)
temperature = np.zeros(np.shape(Fext))
temperature[0, :] = 0.0
# temperature[-1,:] = sigmoid(time_list)
temperature[-1,:] = 1.0

# Solve transient heat problem
solve_transient_heat_1D(properties, DB=basis_cgg, W=weight_cgg, Fext=Fext, 
                        time_list=time_list, dof=dof, dod=dod, Tinout=temperature)

# ------------------
# Post-treatement
# ------------------
# Create eval points
nbknots = 101
knots = np.linspace(0, 1, nbknots)
DB = eval_basis_python(degree, knotvector, knots)
temperature_interp = DB[0].T @ temperature

# Create mesh space-time
XX, TIME = np.meshgrid(knots*JJ, time_list)

# Plot figure
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,4))
im = ax.contourf(XX, TIME, temperature_interp.T, 20)
cbar = plt.colorbar(im)
cbar.set_label('Temperature (K)')

ax.grid(None)
ax.set_ylabel('Time (s)')
ax.set_xlabel('Position (m)')
fig.tight_layout()
fig.savefig(folder + 'TransientHeat_1D.png')

# -----------------------

# Save data to compare with 3D
pos = int((nbknots-1)/2)
midpoint_temp = temperature_interp[pos, :]
data2save = np.column_stack((time_list, midpoint_temp))
np.savetxt(folder + 'data1D.dat', data2save)