"""
.. Test of simulation
.. We test a single example of what is meant to do simulation script
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.physics import *
from lib.post_treat_methods import (ThermalSimulation, 
                                    SimulationData, 
                                    plot_iterative_solver
)

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test6/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
dataExist    = True
isIGA        = False
degree_list  = np.arange(6, 7)
cuts_list    = np.arange(5, 6)
method_list  = ["WP", "C", "JMC", "TDC"]
geoName_list = ['CB', 'VB', 'TR', 'RQA']

for cuts in cuts_list:
    for degree in degree_list:
        for geoName in geoName_list: 

            if geoName   == 'CB': funpow, funtemp = powden_cube, None 
            elif geoName == 'VB': funpow, funtemp = powden_prism, None 
            elif geoName == 'TR': funpow, funtemp = powden_thickring, None 
            elif geoName == 'RQA': funpow, funtemp = powden_rotring, temperature_rotring 

            thermalinputs = {'degree': degree, 'cuts': cuts, 'name': geoName, 'isIGA': isIGA, 
                            'funPowDen': funpow, 'funTemp': funtemp, 'IterMethods': method_list}
            Simulation = ThermalSimulation(thermalinputs, folder)  
            filename = Simulation._filename

            if not dataExist:
                conductivity = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
                Dirichlet = {'thermal':np.array([[1, 1], [1, 1], [1, 1]])}
                material = {'capacity':1.0, 'conductivity': conductivity}
                Simulation.run_simulation(material=material, Dirichlet=Dirichlet)

            else :
                inputs = SimulationData(filename)._dataSimulation
                plot_iterative_solver(filename, inputs, extension= '.pdf')
