"""
.. Test of elastoplasticity 1D
.. Joaquin Cornejo 
"""

import pickle
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformMaxregularKnotvector
from pysrc.lib.thermomecha1D import mechamat1D, plot_results

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceVol(P:list):
	force = 0.4*np.sin(P/1e3)
	return force

# Set global variables
nbSteps = 50
geoArgs = {'length': 1.e3}
matArgs = {'elastic_modulus':2e5, 'elastic_limit':100,
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}

degree, nbel = 8, 1024
knotvector   = createUniformMaxregularKnotvector(degree, nbel, multiplicity=degree)

# Create geometry
quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
modelPhy = mechamat1D(args)

with open(folder + 'refpart.pkl', 'wb') as outp:
    pickle.dump(modelPhy, outp, pickle.HIGHEST_PROTOCOL)

# Add material
modelPhy.activate_mechanical(matArgs)

# Add boundary condition
modelPhy.add_DirichletCondition(table=[1, 0])

# Define boundaries conditions
Fext    = np.zeros((modelPhy.nbctrlpts, 2*nbSteps + 1))
Fextref = modelPhy.compute_volForce(forceVol(modelPhy.qpPhy))
for i in range(0, nbSteps+1): Fext[:, i] = i/nbSteps*Fextref
for i in range(nbSteps+1, 2*nbSteps+1): Fext[:, i] = (2*nbSteps - i)/nbSteps*Fextref

# Solve
disp_cp, strain_qp, stress_qp, plastic_qp, Cep_qp = modelPhy.solve(Fext=Fext)
np.save(folder+'disp', disp_cp)
plastic_cp  = modelPhy.L2projectionCtrlpts(plastic_qp)
stress_cp 	= modelPhy.L2projectionCtrlpts(stress_qp)
plot_results(modelPhy.quadRule, geoArgs['length'], disp_cp, plastic_cp, stress_cp, folder=folder, method=quadArgs['quadrule'])