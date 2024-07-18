"""
	In this file we develop the ideas in structural analysis
	FILE DEPRECATED
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.others.lib_structanalysis import Timoshenko

full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)
# raise Warning('VERIFY FIRST TIMOSHENKO LIBRARY')

# Set global variables
length       = 1.0
degree, nbel = 4, 10 
geometry = createUniformOpenCurve(degree, nbel, length)
modelPhy = part1D(geometry, kwargs={'quadArgs': {'quadrule': 'wq'}})

# Create Timoshenko beam
args    = {'geoArgs': {'section': 'rectangle'}}
problem = Timoshenko(modelPhy, args)

# Add material
matArgs   = {'elastic_modulus':1e8, 'elastic_limit':506, 'poisson_ratio': 0.3}
problem.activate_mechanical(matArgs)

# Add boundary condition
problem.add_DirichletCondition(0, values=[0.0, 0.0, 0.0], table=[True, True, True])
problem.add_DirichletCondition(-1, values=[0.0, 0.0, 0.0], table=[True, True, True])

# Compute external force vector of uniform force field
Fu, Qw = 0.0, 4.0
Fext   = problem.compute_volforce(Fu, Qw)

# Solve 
u, w, theta = problem.solve(Fext)

# Post-treatement
u_interp = problem.interpolateMeshgridField(u)[0]
w_interp, x_interp = problem.interpolateMeshgridField(w)
beam = x_interp + u_interp

fig, ax = plt.subplots()
ax.plot(beam, w_interp)
ax.set_ylabel('Deflection (m)')
ax.set_xlabel('Position (m)')
fig.tight_layout()
fig.savefig(folder + 'Timoshenko.png')