from lib.__init__ import *
from lib.lib_base import (createKnotVector)
from lib.thermomecha1D import mechamat1D, plot_results
from lib.lib_load import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/t1dim/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
mechaprop    = {'law': 'swift', 'elastic_modulus':200e3, 'elastic_limit':506, 'K':2e4, 'exp':0.5}
# mechaprop    = {'law': 'linear', 'elastic_modulus':200e3, 'elastic_limit':506, 'Hbar':1445, 'theta':1.0}
# mechaprop    = {'law': 'voce', 'elastic_modulus':200e3, 'elastic_limit':506, 'Hbar':1445, 'delta':65.8, 'Kinf': 272, 'theta':1.0}
length       = 1.0
nbSteps      = 101

# ---------------
# IGA
# ---------------
degree, nbel = 4, 22
knotvector   = createKnotVector(degree, nbel)
kwargs = {'length': length, 'degree': degree, 'knotvector': knotvector,
            'quadrule': 'iga', 'property': mechaprop, 'quadmethod': 'leg'}
model = mechamat1D(**kwargs)
model.set_DirichletCondition(table=[1, 0])
nbqpiga = model._nbqp

# Define boundaries conditions
Fext        = np.zeros((model._nbctrlpts, nbSteps))
Fext[:, -1] = model.compute_volForce(forceVol(model._qpPar))
for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

# Solve 
disp_iga, strain_iga, stress_iga, plastic_iga, Cep_iga = model.solve(Fext=Fext)
strain_cp   = model.interpolate_CPfield(strain_iga)
plastic_cp  = model.interpolate_CPfield(plastic_iga)
stress_cp 	= model.interpolate_CPfield(stress_iga)
plot_results(degree, knotvector, length, disp_iga, plastic_cp,
                stress_cp, folder=folder, method='IGA')
disp_iga_interp = model.interpolate_sample(disp_iga[:, 1:])[0]

# -------------------
# Weighted Quadrature
# -------------------
degree, nbel = 4, 51
knotvector   = createKnotVector(degree, nbel)
kwargs = {'length': length, 'degree': degree, 'knotvector': knotvector,
            'quadrule': 'wq', 'property': mechaprop}
# kwargs = {'length': length, 'degree': degree, 'knotvector': knotvector,
#             'quadrule': 'iga', 'property': mechaprop, 'quadmethod': 'lob'}
model = mechamat1D(**kwargs)
model.set_DirichletCondition(table=[1, 0])
nbqpwq = model._nbqp

# Define boundaries conditions
Fext        = np.zeros((model._nbctrlpts, nbSteps))
Fext[:, -1] = model.compute_volForce(forceVol(model._qpPar))
for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

# Solve
disp_wq, strain_wq, stress_wq, plastic_wq, Cep_wq = model.solve(Fext=Fext)
strain_cp   = model.interpolate_CPfield(strain_wq)
plastic_cp  = model.interpolate_CPfield(plastic_wq)
stress_cp 	= model.interpolate_CPfield(stress_wq)
plot_results(degree, knotvector, length, disp_wq, plastic_cp,
                stress_cp, folder=folder, method='WQ')
disp_wq_interp = model.interpolate_sample(disp_wq[:, 1:])[0]

# ------------------
# Post-treatement
# ------------------
fig, [ax0, ax1] = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
ax0.plot(model._qpPar, stress_wq[:, 51])
ax0.set_ylabel('Stress (MPa)')
ax1.plot(model._qpPar, Cep_wq[:, 51])
ax1.set_ylabel('Tangent modulus (MPa)')
for ax in [ax0, ax1]:
    ax.set_ylim(bottom=0.0)
    ax.set_xlim(left=0.0, right=1.0)
    ax.set_xlabel('Quadrature point position')
fig.tight_layout()
fig.savefig(folder+'data_step51')

error    = disp_iga_interp - disp_wq_interp
relerror = np.linalg.norm(error, axis=0)/np.linalg.norm(disp_iga_interp, axis=0)
fig, ax  = plt.subplots(nrows=1, ncols=1)
ax.semilogy(relerror*100)
ax.set_ylim(bottom=1e-12, top=1e0)
ax.set_ylabel('L2 Relative error (\%)')
ax.set_xlabel('Step')
fig.tight_layout()
fig.savefig(folder + 'Relative_error_distance.png')

# relerror    = (np.abs(stress_iga[int((nbqpiga-1)/2), 1:] - stress_wq[int((nbqpwq-1)/2), 1:])/
#             np.abs(stress_iga[int((nbqpiga-1)/2), 1:])
# )
# fig, ax  = plt.subplots(nrows=1, ncols=1)
# ax.semilogy(relerror*100)
# ax.set_ylim(bottom=1e-12, top=1e0)
# ax.set_ylabel('Relative error (\%)')
# ax.set_xlabel('Step')
# fig.tight_layout()
# fig.savefig(folder + 'Relative_error_stress.png')