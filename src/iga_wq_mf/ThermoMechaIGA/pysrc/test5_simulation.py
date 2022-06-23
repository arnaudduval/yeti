"""
.. Test of simulation
.. We test a single example of what is meant to do simulation script
.. Joaquin Cornejo 
"""

# Python libraries
import os

# My libraries
from lib.physics import (powden_cube, 
                        powden_prism,
                        powden_thickring, 
                        powden_rotring, 
)
from lib.post_treat_methods import ThermalSimulation, SimulationData, plot_iterative_solver

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test5/'
if not os.path.isdir(folder): os.mkdir(folder)

# Some constants
FileExist = False
GEOMETRY_CASE = 'RQA'
DEGREE, CUTS = 3, 6
IS_IGA = False
if IS_IGA: is_cg_list = [True]
else: is_cg_list = [False]
method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]

for IS_CG in is_cg_list:
    
    # Get file name
    if GEOMETRY_CASE == 'CB':   funpow, funtemp = powden_cube, None 
    elif GEOMETRY_CASE == 'VB': funpow, funtemp = powden_prism, None 
    elif GEOMETRY_CASE == 'TR': funpow, funtemp = powden_thickring, None 
    elif GEOMETRY_CASE == 'RQA': funpow, funtemp = powden_rotring, None 
    
    # Run simulation
    thermalinputs = {'degree': DEGREE, 'cuts': CUTS, 'case': GEOMETRY_CASE, 'isIGA': IS_IGA, 'isCG': IS_CG, 
                'funPowDen': funpow, 'funTemp': funtemp, 'IterMethods': method_list}
    Simulation = ThermalSimulation(thermalinputs, folder)  
    filename = Simulation._filename

    # Run simulation
    if not FileExist:
        Simulation.run_simulation()

    else :
        Data = SimulationData(filename)
        output = Data._dataIter
        plot_iterative_solver(filename, output)

