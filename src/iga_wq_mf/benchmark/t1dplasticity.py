"""
.. Test of elastoplasticity 1D
.. Joaquin Cornejo 
"""

import pickle
from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_1d import mechamat1D

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1elastoplasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceVol(P:list):
	force = 3.2e8*P
	return force

# Create geometry
length = 1
degree, nbel = 2, 256
crv = createUniformCurve(degree, nbel, length)
args = {'quadArgs': {'quadrule': 'iga', 'type': 'leg'}}
modelPhy = mechamat1D(crv, args)

with open(folder + 'refpart.pkl', 'wb') as outp:
    pickle.dump(modelPhy, outp, pickle.HIGHEST_PROTOCOL)

# Add material
matArgs = {'elastic_modulus':2e11, 'elastic_limit':1e8, 'plasticLaw': {'Isoname': 'linear', 'Eiso':2e8}}
modelPhy.activate_mechanical(matArgs)

# Add boundary condition
modelPhy.add_DirichletCondition(table=[1, 1])

# Define boundaries conditions
nbsteps = 151
time_list = np.linspace(0, np.pi, nbsteps)
Fend = np.atleast_2d(modelPhy.compute_volForce(forceVol)).transpose()
Fext = np.kron(Fend, np.sin(time_list))

# Solve
disp_cp, strain_qp, stress_qp, plastic_qp, Cep_qp = modelPhy.solve(Fext=Fext)
np.save(folder+'disp', disp_cp)