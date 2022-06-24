"""
.. Test of simulation
.. We test a single example of what is meant to do simulation script
.. Joaquin Cornejo 
"""

# Python libraries
import os
import numpy as np

# My libraries
from lib.physics import (powden_cube, 
                        powden_prism,
                        powden_thickring, 
)
from lib.post_treat_methods import ThermalSimulation, SimulationData, plot_iterative_solver

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test5/'
if not os.path.isdir(folder): os.mkdir(folder)

# Some constants
FileExist = False
DEGREE, CUTS = 3, 5
IS_IGA = False
IS_CG = False
method_list = ["WP", "C", "TDS", "JMS", "TDC", "JMC"]

for GEOMETRY_CASE in ['VB', 'CB', 'TR']:
    
    # Get file name
    if GEOMETRY_CASE == 'CB':   funpow, funtemp = powden_cube, None 
    elif GEOMETRY_CASE == 'VB': funpow, funtemp = powden_prism, None 
    elif GEOMETRY_CASE == 'TR': funpow, funtemp = powden_thickring, None 
    
    # Run simulation
    thermalinputs = {'degree': DEGREE, 'cuts': CUTS, 'case': GEOMETRY_CASE, 'isIGA': IS_IGA, 'isCG': IS_CG, 
                'funPowDen': funpow, 'funTemp': funtemp, 'IterMethods': method_list}
    Simulation = ThermalSimulation(thermalinputs, folder)  
    filename = Simulation._filename

    # Run simulation
    if not FileExist:
        conductivity = np.array([[1, -0.5, -0.1],[-0.5, 1, -0.25], [-0.1, -0.25, 1]])
        properties = {"conductivity": conductivity}
        Simulation.run_simulation(**properties)

    else :
        Data = SimulationData(filename)
        output = Data._dataIter
        plot_iterative_solver(filename, output)

