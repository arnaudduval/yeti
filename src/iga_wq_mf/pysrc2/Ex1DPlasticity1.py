from lib.__init__ import *
from lib.lib_base import (createKnotVector)
from lib.thermomecha1D import mechamat1D, plot_results
from lib.lib_load import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t1dim/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
length       = 1.0
degree, nbel = 4, 51
knotvector   = createKnotVector(degree, nbel)

mechaprop = {'elastic_modulus':200e3, 'elastic_limit':506, 'Hbar':1445, 'theta':1.0}
# mechaprop = {'elastic_modulus':200e3, 'elastic_limit':506, 'Hbar':1445, 'delta':65.8, 'Kinf': 272, 'theta':1.0}
# mechaprop = {'elastic_modulus':200e3, 'elastic_limit':506, 'K':2e4, 'exp':0.5}

kwargs = {'length': 1.0, 'degree': degree, 'knotvector': knotvector,
        'quadrule': 'wq', 'law': 'linear', 'property': mechaprop}
model = mechamat1D(**kwargs)
model.set_DirichletCondition(table=[1, 0])

# Define boundaries conditions
nbSteps = 101
Fext        = np.zeros((model._nbctrlpts, nbSteps))
Fext[:, -1] = model.compute_volForce(forceVol(model._qpPar))
for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

# Solve using IGA
disp_iga, strain_iga, stress_iga, plastic_iga, Cep_iga = model.solve(Fext=Fext)
strain_cp   = model.interpolate_CPfield(strain_iga)
plastic_cp  = model.interpolate_CPfield(plastic_iga)
stress_cp 	= model.interpolate_CPfield(stress_iga)
plot_results(degree, knotvector, length, disp_iga, plastic_cp,
                stress_cp, folder=folder, method='IGA')

# ------------------
# Post-treatement
# ------------------
fig, [ax0, ax1] = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
ax0.plot(model._qpPar, stress_iga[:, 51])
ax0.set_ylabel('Stress (MPa)')
ax1.plot(model._qpPar, Cep_iga[:, 51])
ax1.set_ylabel('Tangent modulus (MPa)')
for ax in [ax0, ax1]:
    ax.set_ylim(bottom=0.0)
    ax.set_xlim(left=0.0, right=1.0)
    ax.set_xlabel('Quadrature point position')
fig.tight_layout()
fig.savefig(folder+'data_step51')
