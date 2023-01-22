from lib.__init__ import *
from lib.base_functions import (eval_basis_python,
								iga_find_positions_weights,
								create_knotvector
)
from lib.D1PlasticityIGA import *
from lib.physics import forceVol

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/examples/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
E, H, sigma_Y, beta, JJ = 200e3, 50e3, 100, 0.5, 1.0
degree, nbel            = 5, 128

nb_ctrlpts = degree + nbel
ctrlpts    = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)

# Get basis and weights in IgA 
qp_cgg, weight_cgg = iga_find_positions_weights(degree, knotvector)
basis_cgg          = eval_basis_python(degree, knotvector, qp_cgg)
properties         = [JJ, E, H, beta, sigma_Y, len(qp_cgg)]

# Define boundaries conditions
N = 101
time_list   = np.linspace(0, 1, N)
dof         = np.arange(1, nb_ctrlpts, dtype=int)
Fext        = np.zeros((nb_ctrlpts, N))
Fext[:, -1] = compute_volForce_1D(basis_cgg, weight_cgg, forceVol(qp_cgg))
Fext[-1,-1] += 400
for i in range(1, N-1): Fext[:, i] = i/(N-1)*Fext[:, -1]

# Solve 
disp, strain, stress = solve_plasticity_1D(properties, DB=basis_cgg, W=weight_cgg, Fext=Fext, dof=dof)
strain_cp   = interpolate_controlPoints_1D(basis_cgg, weight_cgg, strain)
stress_cp 	= interpolate_controlPoints_1D(basis_cgg, weight_cgg, stress)

filename = folder + 'disp_iga' + '.dat'
np.savetxt(filename, disp)

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
fig.savefig(folder + 'ElastoPlasticity.png')

# Plot stress-strain of single point
fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
for ax, pos in zip([ax1, ax2, ax3], [25, 50, 75]):
	ax.plot(strain_interp[pos, :]*100, stress_interp[pos, :])
	ax.set_ylabel('Stress (MPa)')
	ax.set_xlabel('Strain (\%)')

fig.tight_layout()
fig.savefig(folder + 'TractionCurve.png')