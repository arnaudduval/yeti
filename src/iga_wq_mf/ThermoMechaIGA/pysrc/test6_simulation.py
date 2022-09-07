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

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test6/'
if not os.path.isdir(folder): os.mkdir(folder)

# Some constants
FileExist = True
degree, cuts = 4, 5
is_iga = False
method_list = ["WP", "C", "TDS", "JMS", "TDC", "JMC"]

for geometryName in ['VB']:  
    
    # Get file name
    if geometryName == 'CB':   funpow, funtemp = powden_cube, None 
    elif geometryName == 'VB': funpow, funtemp = powden_prism, None 
    elif geometryName == 'TR': funpow, funtemp = powden_thickring, None 

    # Run simulation
    thermalinputs = {'degree': degree, 'cuts': cuts, 'case': geometryName, 'isIGA': is_iga, 
                'funPowDen': funpow, 'funTemp': funtemp, 'IterMethods': method_list}
    Simulation = ThermalSimulation(thermalinputs, folder)  
    filename = Simulation._filename

    # Run simulation
    if not FileExist:
        conductivity = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
        # conductivity = np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
        material = {"conductivity": conductivity}
        Simulation.run_simulation(material=material, Dirichlet=Dirichlet)

    else :
        Data = SimulationData(filename)
        output = Data._dataIter
        plot_iterative_solver(filename, output)

