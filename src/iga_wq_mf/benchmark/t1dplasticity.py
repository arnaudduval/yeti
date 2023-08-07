"""
.. Test of elastoplasticity 1D
.. We test how elasticity module works
.. SI (Steel) : 
..      - Stress : MPa
..      - Length : mm
..      - Force  : N
..      - Mass   : metric ton 
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformMaxregularKnotvector
from pysrc.lib.thermomecha1D import mechamat1D, plot_results

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1plasticity/'
if not os.path.isdir(folder): os.mkdir(folder)

def forceVol(P:list):
	force = 0.4*np.sin(P/1e3)
	return force

# Set global variables
samplesize = 2500
nbSteps    = 50
matArgs    = {'elastic_modulus':2e5, 'elastic_limit':100,
			'plasticLaw': {'name': 'swift', 'K':2e4, 'exp':0.5}}

degree, nbel = 4, 32
knotvector   = createUniformMaxregularKnotvector(degree, nbel, multiplicity=degree)

# Create geometry
quadArgs  = {'degree': degree, 'knotvector': knotvector, 'quadrule': 'iga', 'type': 'leg'}
geoArgs   = {'length': 1.e3}
args  = {'quadArgs': quadArgs, 'geoArgs': geoArgs}
model = mechamat1D(args)

# Add material
model.activate_mechanical(matArgs)

# Add boundary condition
model.add_DirichletCondition(table=[1, 0])

# Define boundaries conditions
Fext    = np.zeros((model.nbctrlpts, 2*nbSteps + 1))
Fextref = model.compute_volForce(forceVol(model.qpPhy))
for i in range(0, nbSteps+1): Fext[:, i] = i/nbSteps*Fextref
for i in range(nbSteps+1, 2*nbSteps+1): Fext[:, i] = (2*nbSteps - i)/nbSteps*Fextref

# Solve
disp_cp, strain_qp, stress_qp, plastic_qp, Cep_qp = model.solve(Fext=Fext)
plastic_cp  = model.L2projectionCtrlpts(plastic_qp)
stress_cp 	= model.L2projectionCtrlpts(stress_qp)
plot_results(model.quadRule, geoArgs['length'], disp_cp, plastic_cp, stress_cp, folder=folder, method=quadArgs['quadrule'])