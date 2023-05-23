from lib.__init__ import *
from lib.lib_base import createKnotVector, sigmoid
from lib.thermomecha1D import thermo1D
from lib.lib_load import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1heat/'
if not os.path.isdir(folder): os.mkdir(folder)

def setKprop(T, prop=0.1):
	# y = prop + prop*np.exp(-0.1*abs(T))
	# y = prop + prop*2.0/(1.0 + np.exp(-5*(T-1.0)))
	y   = np.ones(len(T))
	return y

def setCprop(T, prop=1.0):
	# y = prop + prop*np.exp(-2.0*abs(T))
	y   = np.ones(len(T))
	return y

# Set global variables
theta        = 1
degree, nbel = 10, 10
knotvector   = createKnotVector(degree, nbel)
heatprop     = {'conductivity': setKprop, 'capacity': setCprop}
kwargs = {'length': 1.0, 'degree': degree, 'knotvector': knotvector,
        'quadrule': 'wq', 'heattheta': theta, 'property': heatprop}

model = thermo1D(**kwargs)
model.add_DirichletCondition(table=[1, 1])

# Define boundaries conditions	
N = 1000
time_list = np.linspace(0, 0.02, N)
print('Time step: %3e' %(time_list.max()/N))
Fprop     = powden(model._qpPar)
FFend     = model.compute_volForce(Fprop)
FFend     = np.atleast_2d(FFend).reshape(-1, 1)
Fext      = np.kron(FFend, sigmoid(time_list))

temperature = np.zeros(np.shape(Fext))
temperature[0, :] = 0.0
temperature[-1,:] = 1.0

# Solve
model.solve(Fext=Fext, time_list=time_list, Tinout=temperature)
temp_interp, x_interp = model.interpolate_sampleField(temperature)
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
