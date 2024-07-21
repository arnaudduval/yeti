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
from pysrc.lib.lib_job1d import heatproblem1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import heatmat

# Set global variables
LENGTH = 1.0
degree, cuts = 4, 5 

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
	f = (2*np.pi*np.cos(2*np.pi*t)*np.sin(2*np.pi*x) 
		+ 4*np.pi**2*np.sin(2*np.pi*t)*np.sin(2*np.pi*x)
		)
	return f

# Create geometry
LENGTH = 1.0
nbel   = int(2**cuts)
geometry = createUniformOpenCurve(degree, nbel, LENGTH)
modelPhy = part1D(geometry, kwargs={'quadArgs': {'quadrule': 'iga'}})

timespan, cuts_time = 1.0, np.copy(cuts)
time_list = np.linspace(0, timespan, int(2**cuts_time)+1)

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
