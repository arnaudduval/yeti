from lib.__init__ import *
from lib.base_functions import (eval_basis_python,
								iga_find_positions_weights,
								create_knotvector, sigmoid
)
from lib.physics import setCprop, setKprop, powden
from lib.D1transientheat import *

# Select folder
full_path = os.path.realpath(__file__)
# folder = os.path.dirname(full_path) + '/results/test7/' !!!!!!!!!!!!
folder = os.path.dirname(full_path) + '/results/stability/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
theta, JJ    = 1, 1.0
degree, nbel = 10, 10
multiplicity = degree

knotvector = create_knotvector(degree, nbel, multiplicity=multiplicity)
nb_ctrlpts = len(knotvector) - degree - 1

# Get basis and weights in IgA
qp_cgg, weight_cgg = iga_find_positions_weights(degree, knotvector)
basis_cgg 	 = eval_basis_python(degree, knotvector, qp_cgg, multiplicity=multiplicity)
properties   = [JJ, setKprop, setCprop, theta]

# Define boundaries conditions	
N = 1000
time_list = np.linspace(0, 0.02, N)
time_step = time_list.max()/N
print('Time step: %3e' %time_step)
dod = [0, -1]
dof = np.arange(1, nb_ctrlpts-1, dtype=int)
Fprop     = powden(qp_cgg)
FFend     = compute_volsource_1D(JJ, basis_cgg, weight_cgg, Fprop)
FFend     = np.atleast_2d(FFend).reshape(-1, 1)
Fext      = np.kron(FFend, sigmoid(time_list))

temperature = np.zeros(np.shape(Fext))
temperature[0, :] = 0.0
temperature[-1,:] = 1.0

# Solve
solve_transient_heat_1D(properties, DB=basis_cgg, W=weight_cgg, Fext=Fext, 
						time_list=time_list, dof=dof, dod=dod, Tinout=temperature)

# ------------------
# Post-treatement
# ------------------
nbknots = 201
knots = np.linspace(0, 1, nbknots)
DB    = eval_basis_python(degree, knotvector, knots, multiplicity=multiplicity)
temp_interp = DB[0].T @ temperature
print(temp_interp.min())

# Plot
XX, TIME = np.meshgrid(knots*JJ, time_list)
fig, ax  = plt.subplots(nrows=1, ncols=1, figsize=(10,4))
levels   = np.array([-0.2]); levels = np.append(levels, np.linspace(0, 1, 9))
norm 	 = mpl.colors.BoundaryNorm(levels, len(levels))
colors   = list(plt.cm.Greys(np.linspace(0, 1, len(levels)-1))); colors[0] = "red"
cmap = mpl.colors.ListedColormap(colors,"", len(colors))
im   = ax.contourf(XX, TIME, temp_interp.T, norm=norm, cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label('Temperature (K)')

ax.grid(False)
ax.set_ylabel('Time (s)')
ax.set_xlabel('Position (m)')
fig.tight_layout()
fig.savefig(folder + 'TransientHeat_1D.png')

# # Save data to compare with 3D
# pos = int((nbknots-1)/2)
# midpoint_temp = temp_interp[pos, :]
# data2save = np.column_stack((time_list, midpoint_temp))
# np.savetxt(folder + 'data1D.dat', data2save)