from lib.__init__ import *
from lib.base_functions import (wq_find_basis_weights_opt,
								eval_basis_python,
								create_knotvector
)
from lib.D1PlasticityWQ import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/examples/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
E, H, sigma_Y, beta, JJ = 200e3, 50e3, 100, 0.5, 1.0
degree, nbel            = 5, 32

nb_ctrlpts = degree + nbel
ctrlpts    = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)

# Get basis and weights in IgA 
qp_wq, B0, B1, W00, W01, W10, W11 = wq_find_basis_weights_opt(degree, knotvector, 2)
properties  = [JJ, E, H, beta, sigma_Y, len(qp_wq)]
DB = [B0, B1]; DW = [W00, W01, W10, W11]

# Define boundaries conditions
N = 101
time_list   = np.linspace(0, 1, N)
dof         = np.arange(1, nb_ctrlpts, dtype=int)
Fext        = np.zeros((nb_ctrlpts, N))
Fext[-1, :] = 400*time_list

# Solve 
disp, strain, stress, energy, internal = solve_plasticity_1D(properties, DB=DB, DW=DW, Fext=Fext, dof=dof)
strain_cp   = interpolate_controlPoints_1D(DB, DW, strain)
stress_cp   = interpolate_controlPoints_1D(DB, DW, stress)

# fig, ax = plt.subplots(nrows=1, ncols=1)
# ax.loglog(np.arange(1, 1+len(internal)), abs(internal))
# fig.tight_layout()
# fig.savefig(folder + 'Energy.png')

filename = folder + 'disp_iga' + '.dat'
disp_ref = np.loadtxt(filename)
error    = disp_ref - disp
relerror = np.linalg.norm(error, np.inf, axis=0)/np.linalg.norm(disp_ref, np.inf, axis=0)
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.semilogy(relerror*100)
ax.set_ylabel('Relative error')
ax.set_xlabel('Step')
fig.tight_layout()
fig.savefig(folder + 'Relative_error.png')

# ------------------
# RESULTS
# ------------------
knots  = np.linspace(0, 1, 101)
DB     = eval_basis_python(degree, knotvector, knots)
strain = interpolate_strain_1D(JJ, DB, disp) 
displacement   = DB[0].T @ disp
strain_interp  = DB[0].T @ strain_cp
stress_interp  = DB[0].T @ stress_cp

# Plot fields
XX, STEPS = np.meshgrid(knots*JJ, np.arange(N))
names = ['Displacement field', 'Total strain field', 'Stress field']
fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14, 4))
for ax, variable, name in zip([ax1, ax2, ax3], [displacement, strain_interp, stress_interp], names):
	ax.contourf(XX, STEPS, variable.T, 20)

	ax.grid(None)
	ax.set_title(name)
	ax.set_ylabel('Step')
	ax.set_xlabel('Position (m)')

fig.tight_layout()
fig.savefig(folder + 'ElastoPlasticity2.png')

# Plot stress-strain of single point
fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
for ax, pos in zip([ax1, ax2, ax3], [25, 50, 75]):
	ax.plot(strain_interp[pos, :]*100, stress_interp[pos, :])
	ax.set_ylabel('Stress (MPa)')
	ax.set_xlabel('Strain (\%)')

fig.tight_layout()
fig.savefig(folder + 'TractionCurve2.png')