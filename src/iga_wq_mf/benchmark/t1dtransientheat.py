"""
.. Test of 1D transient heat transfer 
.. Author: Joaquin Cornejo
density: 7800 kg/m3
capacity: 460 J/(kg.K)
conductivity: 55 W/(m.K)
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve, sigmoid
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_1djob import heatproblem1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1transferheat/'
if not os.path.isdir(folder): os.mkdir(folder)

c = 0.1

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	y = np.ones(len(temperature))
	return y

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	y = np.ones(len(temperature))
	return y

def powerDensity(args:dict):
	x = args['position']; t = args['time']
	f = (2*c*np.pi*np.cos(2*np.pi*t)*np.sin(2*np.pi*x) 
		+ 4*c*np.pi**2*np.sin(2*np.pi*t)*np.sin(2*np.pi*x)
		)
	return f

# Set global variables
length = 1.0
degree, cuts = 4, 5 
cuts_time = 5

# Create geometry
length = 1.0
nbel   = int(2**cuts)
geometry = createUniformOpenCurve(degree, nbel, length)
modelPhy = part1D(geometry, kwargs={'quadArgs': {'quadrule': 'iga'}})

if cuts_time is None: cuts_time = np.copy(cuts)
timespan = 1.0
nbsteps  = int(2**cuts_time)
time_list = np.linspace(0, timespan, nbsteps+1)

# Add material
material = heatmat()
material.addConductivity(conductivityProperty, isIsotropic=False)
material.addCapacity(capacityProperty, isIsotropic=False)

# Add boundary condition
boundary = boundaryCondition(modelPhy.nbctrlpts)
boundary.add_DirichletConstTemperature(table=np.array([[1, 1]]))
temperature = np.zeros((modelPhy.nbctrlpts_total, len(time_list)))

# Create heat problem
problem = heatproblem1D(material, modelPhy, boundary)

# Define external force	
Fext_list = np.zeros((modelPhy.nbctrlpts_total, len(time_list)))
for i, t in enumerate(time_list):
	Fext_list[:, i] = problem.compute_volForce(powerDensity, args={'position':problem.part.qpPhy, 'time':t})

# Solve
problem.solveFourierTransientProblem(temperature, Fext_list, time_list, isLumped=False, alpha=1.0)

# Post-processing
temp_interp, x_interp = problem.interpolateMeshgridField(temperature, sampleSize=201)
XX, TIME = np.meshgrid(x_interp, time_list)
fig, ax  = plt.subplots(figsize=(10, 4))
im   = ax.contourf(XX, TIME, temp_interp.T, cmap='viridis')
cbar = plt.colorbar(im)
cbar.set_label(r'$\displaystyle\frac{u - u_0}{u_1 - u_0}$')

ax.grid(False)
ax.set_ylabel(r'$\displaystyle\frac{\tau}{T_{s}}$')
ax.set_xlabel(r'$\displaystyle\frac{x}{L}$')
fig.tight_layout()
fig.savefig(folder + 'TransientHeat1D.png')
