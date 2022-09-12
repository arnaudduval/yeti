"""
.. Test of simulation
.. We test a single example of what is meant to do simulation script
.. Joaquin Cornejo 
"""

# Python libraries
import os, numpy as np

# My libraries
from lib.physics import (powden_cube, 
                        powden_prism,
                        powden_thickring, 
)
from lib.post_treat_methods import (ThermalSimulation, 
                                    SimulationData, 
                                    plot_iterative_solver
)

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test6/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
FileExist = True
isIGA = False
degree, cuts = 6, 5
method_list = ["WP", "C", "TDS", "TDC", "JMS", "JMC"]

for geometryName in ['VB']:  
    
    # Get file name
    if geometryName   == 'CB': funpow, funtemp = powden_cube, None 
    elif geometryName == 'VB': funpow, funtemp = powden_prism, None 
    elif geometryName == 'TR': funpow, funtemp = powden_thickring, None 

    # Run simulation
    thermalinputs = {'degree': degree, 'cuts': cuts, 'case': geometryName, 'isIGA': isIGA, 
                    'funPowDen': funpow, 'funTemp': funtemp, 'IterMethods': method_list}
    Simulation = ThermalSimulation(thermalinputs, folder)  
    filename = Simulation._filename

    if not FileExist:
        conductivity = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
        # conductivity = np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
        material = {'capacity':1.0, 'conductivity': conductivity}
        Simulation.run_simulation(material=material, Dirichlet=Dirichlet)

    else :
        inputs = SimulationData(filename)._dataSimulation
        plot_iterative_solver(filename, inputs)

