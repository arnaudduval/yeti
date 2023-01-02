from lib.StructAnalysis import *
from lib.base_functions import eval_basis_python

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/examples/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
nu, E, JJ    = 0.3, 1e8, 1.0
degree, nbel = 4, 10 

# Create Timoshenko beam
knotvector = create_knotvector(degree, nbel)
geoname    = 'rectangle'
properties = {'E': E, 'nu': nu, 'degree': degree, 'knotvector': knotvector, 'length': JJ}
model      = Timoshenko(name=geoname, properties=properties)

# Block boundaries
model.block_boundary(0, values=[0.0, 0.0, 0.0], table=[True, True, True])
# model.block_boundary(-1, values=[0.0, 0.0, 0.0], table=[False, True, False])

# Compute external force vector of uniform force field
Fu, Qw = 0.0, -4.0
Fext = model.compute_timoshenko_force(Fu, Qw)

# Solve 
u, w, theta = model.solve_timoshenko(Fext)

# ----------------
#  RESULTS  
# ----------------
knots = np.linspace(0, 1, 101)
DB    = eval_basis_python(degree, knotvector, knots)
u_interp = DB[0].T @ u
w_interp = DB[0].T @ w
beam = knots*JJ + u_interp

# Plot
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(beam, w_interp)
ax.set_ylabel('Deflection (m)')
ax.set_xlabel('Position (m)')
fig.tight_layout()
fig.savefig(folder + 'Timoshenko.png')

