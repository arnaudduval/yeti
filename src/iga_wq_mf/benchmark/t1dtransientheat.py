"""
.. Test of 1D transient heat transfer 
.. Author: Joaquin Cornejo
density: 7800 kg/m3
capacity: 460 J/(kg.K)
conductivity: 55 W/(m.K)
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve, sigmoid
from pysrc.lib.lib_1d import heattransfer1D

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1transferheat/'
if not os.path.isdir(folder): os.mkdir(folder)

def conductivityProperty(T):
	y = 55*100/(7800*460)*np.ones(len(T))
	return y

def capacityProperty(T):
	y = np.ones(len(T))
	return y

# Set global variables
length = 1.0
degree, nbel = 4, 200 
crv = createUniformCurve(degree, nbel, length)

# Create geometry
args     = {'quadArgs': {'quadrule': 'iga'}}
modelPhy = heattransfer1D(crv, args)

# Add material 
matArgs = {'heattheta': 1.0, 'conductivity': conductivityProperty, 'capacity': capacityProperty}
modelPhy.activate_thermal(matArgs)

# Add boundary condition
modelPhy.add_DirichletCondition(table=[1, 1])

# Define external force	
nbsteps   = 100
time_list = np.linspace(0, 1, nbsteps)
print('Time step: %3e' %(time_list.max()/nbsteps))
Fend = np.zeros((modelPhy.nbctrlpts, 1))
Fext = np.kron(Fend, sigmoid(time_list))

temperature = np.zeros(np.shape(Fext))
temperature[-1, 1:] = 1

# Solve
modelPhy.solve(Fext=Fext, time_list=time_list, Tinout=temperature, isLumped=True)
temp_interp, x_interp = modelPhy.interpolateMeshgridField(temperature)
print(temp_interp.min(), temp_interp.max())

# ------------------
# Post-treatement
# ------------------
XX, TIME = np.meshgrid(x_interp, time_list)
fig, ax  = plt.subplots(figsize=(10, 4))
levels   = np.array([-0.1]); levels = np.append(levels, np.linspace(0, 1, 9))
norm 	 = mpl.colors.BoundaryNorm(levels, len(levels))
colors   = list(plt.cm.Greys(np.linspace(0, 1, len(levels)+1))); colors[0] = "red"
cmap = mpl.colors.ListedColormap(colors,"", len(colors))
im   = ax.contourf(XX, TIME, temp_interp.T, norm=norm, cmap=cmap)
cbar = plt.colorbar(im, )
cbar.set_label(r'$\displaystyle\frac{u-u_{\infty}}{u_{\infty}}$')

ax.grid(False)
ax.set_ylabel(r'$\displaystyle\frac{\tau}{T_{s}}$')
ax.set_xlabel(r'$\displaystyle\frac{x}{L}$')
fig.tight_layout()
fig.savefig(folder + 'TransientHeat1D.png')
