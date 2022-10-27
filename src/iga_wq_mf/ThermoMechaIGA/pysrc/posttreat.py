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
isIGA = False
degree_list, cuts_list = np.arange(3, 7), np.arange(7, 8)
method_list = ["WP", "C", "JMC", "TDC"]

for geometryName in ['CB', 'VB', 'TR', 'RQA']: 
    nbIterList = np.zeros((len(cuts_list), len(degree_list)))
    timeList = np.zeros((len(cuts_list), len(degree_list)))
    for m, cuts in enumerate(cuts_list):
        for n, degree in enumerate(degree_list):
        
            funpow, funtemp = None, None 

            # Run simulation
            thermalinputs = {'degree': degree, 'cuts': cuts, 'case': geometryName, 'isIGA': isIGA, 
                            'funPowDen': funpow, 'funTemp': funtemp, 'IterMethods': method_list}
            Simulation = ThermalSimulation(thermalinputs, folder)  
            filename = Simulation._filename

            inputs = SimulationData(filename)._dataSimulation
            i = inputs['Methods'].index('JMC')
            timeSolver = inputs['TimeIter'][i] - inputs['TimeNoIter'][i]

            # Define inputs
            residue = inputs["Res"]
            residue = np.asarray(residue[i])
            residue = residue[residue>1e-12]

            savename1 = folder + geometryName + 'Iter2' + '.txt'
            savename2 = folder + geometryName + 'time2' + '.txt'
            nbIterList[m, n] = len(residue)
            timeList[m, n] = timeSolver
    np.savetxt(savename1, nbIterList)
    np.savetxt(savename2, timeList)

