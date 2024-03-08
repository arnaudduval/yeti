"""
.. Test of 1D transient heat transfer 
.. Author: Joaquin Cornejo
density: 7800 kg/m3
capacity: 460 J/(kg.K)
conductivity: 55 W/(m.K)
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve, sigmoid
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_1d import heatproblem1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1transferheat/'
if not os.path.isdir(folder): os.mkdir(folder)

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	y = 55*100/(7800*460)*np.ones(len(temperature))
	return y

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	y = np.ones(len(temperature))
	return y

def relaxationProperty(args:dict):
	temperature = args.get('temperature')
	y = np.ones(len(temperature))
	return y

# Set global variables
length = 1.0
degree, nbel = 6, 512 
nbsteps = 100
time_list = np.linspace(0, 1, nbsteps)
print('Time step: %3e' %(time_list.max()/nbsteps))

# Create geometry
geometry = createUniformCurve(degree, nbel, length)
modelPhy = part1D(geometry, kwargs={'quadArgs': {'quadrule': 'iga'}})

# Add boundary condition
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.array([[1, 1]]))
temperature = np.zeros((modelPhy.nbctrlpts_total, len(time_list)))
temperature[-1, :] = 1

# Add material
material = heatmat()
material.addConductivity(conductivityProperty, isIsotropic=False)
material.addCapacity(capacityProperty, isIsotropic=False)
material.addRelaxation(relaxationProperty, isIsotropic=False)

# Create heat problem
problem = heatproblem1D(material, modelPhy, boundary)

# Define external force	
Fext_list = np.zeros((modelPhy.nbctrlpts_total, len(time_list)))

# Solve
# problem.solveFourierTransientProblem(temperature, Fext_list, time_list, isLumped=False)
problem.solveCattaneoTransientProblem(temperature, Fext_list, time_list, isLumped=False)
temp_interp, x_interp = problem.interpolateMeshgridField(temperature)
print(temp_interp.min(), temp_interp.max())

# ------------------
# Post-treatement
# ------------------
XX, TIME = np.meshgrid(x_interp, time_list)
fig, ax  = plt.subplots(figsize=(10, 4))
levels   = np.array([-0.1]); levels = np.append(levels, np.linspace(0, 1, 9))
norm 	 = mpl.colors.BoundaryNorm(levels, len(levels))
colors   = list(plt.cm.Greys(np.linspace(0, 1, len(levels)+1))); colors[0] = 'red'
cmap = mpl.colors.ListedColormap(colors, '', len(colors))
im   = ax.contourf(XX, TIME, temp_interp.T, norm=norm, cmap=cmap)
cbar = plt.colorbar(im, )
cbar.set_label(r'$\displaystyle\frac{u - u_0}{u_1 - u_0}$')

ax.grid(False)
ax.set_ylabel(r'$\displaystyle\frac{\tau}{T_{s}}$')
ax.set_xlabel(r'$\displaystyle\frac{x}{L}$')
fig.tight_layout()
fig.savefig(folder + 'TransientHeat1DCattaneo.png')
