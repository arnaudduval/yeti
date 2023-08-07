"""
.. Test of 1D transient heat transfer 
.. Author: Joaquin Cornejo
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformMaxregularKnotvector, sigmoid
from pysrc.lib.thermomecha1D import thermo1D

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1heat/'
if not os.path.isdir(folder): os.mkdir(folder)

def setKprop(T):
	y = np.ones(len(T))
	return y

def setCprop(T):
	y = np.ones(len(T))
	return y

# Set global variables
length       = 1.0
degree, nbel = 6, 100 
knotvector   = createUniformMaxregularKnotvector(degree, nbel)

# Create geometry
quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'wq'}
args      = {'quadArgs': quadArgs, 'geoArgs': {'length': 1.0}}
model 	  = thermo1D(args)

# Add material 
matArgs = {'heattheta': 1.0, 'conductivity': setKprop, 'capacity': setCprop}
model.activate_thermal(matArgs)

# Add boundary condition
model.add_DirichletCondition(table=[1, 1])

# Define external force	
N = 20
time_list = np.linspace(0, 0.02, N)
print('Time step: %3e' %(time_list.max()/N))
Fend = np.zeros((model.nbctrlpts, 1))
Fext = np.kron(Fend, sigmoid(time_list))

temperature = np.zeros(np.shape(Fext))
temperature[0, :] = 0.0
temperature[-1,:] = 1.0

# Solve
model.solve(Fext=Fext, time_list=time_list, Tinout=temperature)
temp_interp, x_interp = model.interpolateMeshgridField(temperature)
print(temp_interp.min())

# ------------------
# Post-treatement
# ------------------
XX, TIME = np.meshgrid(x_interp, time_list)
fig, ax  = plt.subplots(figsize=(10,4))
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
fig.savefig(folder + 'TransientHeat1D.png')
