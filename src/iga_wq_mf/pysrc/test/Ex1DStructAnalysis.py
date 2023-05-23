from lib.structanalysislib import *
from lib.lib_base import createKnotVector

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1plasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
length       = 1.0
degree, nbel = 4, 10 
knotvector   = createKnotVector(degree, nbel)

# Create Timoshenko beam
quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'wq'}
geoArgs   = {'length': length, 'section': 'rectangle'}
args      = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
model     = Timoshenko(args)

# Add material
matArgs   = {'elastic_modulus':1e8, 'elastic_limit':506, 'poisson_ratio': 0.3}
model.activate_mechanical(matArgs)

# Add boundary condition
model.add_DirichletCondition(0, values=[0.0, 0.0, 0.0], table=[True, True, True])
model.add_DirichletCondition(-1, values=[0.0, 0.0, 0.0], table=[True, True, True])

# Compute external force vector of uniform force field
Fu, Qw = 0.0, 4.0
Fext   = model.compute_volforce(Fu, Qw)

# Solve 
u, w, theta = model.solve(Fext)

# ------------------
# Post-treatement
# ------------------
u_interp = model.interpolate_sampleField(u)[0]
w_interp, x_interp = model.interpolate_sampleField(w)
beam = x_interp + u_interp

fig, ax = plt.subplots()
ax.plot(beam, w_interp)
ax.set_ylabel('Deflection (m)')
ax.set_xlabel('Position (m)')
fig.tight_layout()
fig.savefig(folder + 'Timoshenko.png')