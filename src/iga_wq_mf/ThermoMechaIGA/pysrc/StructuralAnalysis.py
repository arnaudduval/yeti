from lib.timoshenko import *
from lib.base_functions import eval_basis_python

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test/'
if not os.path.isdir(folder): os.mkdir(folder)

# ----------------
#  MODELIZE
# ----------------
# Define material (isotropic)
nu, E = 0.0, 1e8
length = 2.0

# Define B-spline 
degree, nbel = 3, 10 
knotvector = create_knotvector(degree, nbel)

# Create Timoshenko object
geoname = 'rectangle'
properties = {'E': E, 'nu': nu, 'degree': degree, 'knotvector': knotvector, 'length': length}
model = Timoshenko(name=geoname, properties=properties)

# ----------------
#  ANALYSIS
# ----------------
# Block boundaries
model.block_boundary(0, values=[0.0, 0.0, 0.0], table=[True, True, True])
model.block_boundary(-1, values=[0.0, 0.0, 0.0], table=[True, True, True])

# Compute external force vector of uniform force field
Fu, Qw = 0.0, 1.0
Fext = model.compute_timoshenko_force(Fu, Qw)

# Solve Timoshenko problem
u, w, theta = model.solve_timoshenko(Fext)

# ----------------
#  RESULTS  
# ----------------
# Create eval points
nbknots = 101
knots = np.linspace(0, 1, nbknots)
DB = eval_basis_python(degree, knotvector, knots)
u_interp = DB[0].T @ u
w_interp = DB[0].T @ w

# Plot figure
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(knots*length+u_interp, w_interp)
ax.set_ylabel('Deflection')
ax.set_xlabel('Position ')
fig.tight_layout()
fig.savefig(folder + 'Timoshenko.png')

