from lib.__init__ import *
from lib.base_functions import (eval_basis_python,
                                iga_find_positions_weights,
                                create_knotvector
)
from lib.D1viscoplasticity import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test/'
if not os.path.isdir(folder): os.mkdir(folder)

# Define mechanical properties
E, H, sigma_Y, beta, JJ = 200e3, 25e3, 250, 0.5, 1

# Define geometry
degree, nbel = 5, 32
nb_ctrlpts = degree + nbel
ctrlpts = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)

# Get basis and weights in IGA analysis
qp_cgg, weight_cgg = iga_find_positions_weights(degree, knotvector)
basis_cgg = eval_basis_python(degree, knotvector, qp_cgg)
nb_qp = len(qp_cgg)
properties = [JJ, E, H, beta, sigma_Y, nb_qp]

# Define boundaries conditions
N = 500
t = np.linspace(0, 1, N)
dof = np.arange(1, nb_ctrlpts, dtype=int)
Fext = np.zeros((nb_ctrlpts, N))
Fext[-1, :] = 400*t

# Solve elastoplasticity problem
disp, epn, sigma = solve_plasticity_1D(properties, DB=basis_cgg, W=weight_cgg, Fext=Fext, dof=dof)
epn_cp = interpolate_controlPoints_1D(basis_cgg, weight_cgg, epn)
sigma_cp = interpolate_controlPoints_1D(basis_cgg, weight_cgg, sigma)

# ------------------
# Post-treatement
# ------------------
# Create eval points
knots = np.linspace(0, 1, 101)
DB = eval_basis_python(degree, knotvector, knots)
displacement = DB[0].T @ disp
strain = interpolate_strain_1D(JJ, DB, disp) 
plastic_strain = DB[0].T @ epn_cp
stress = DB[0].T @ sigma_cp

# Create space-time mesh
XX, STEPS = np.meshgrid(knots*JJ, np.arange(N))

# Plot figure
names = ['Displacement field', 'Plastic strain field', 'Stress field']
fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
for ax, variable, name in zip([ax1, ax2, ax3], [displacement, plastic_strain, stress], names):
    ax.contourf(XX, STEPS, variable.T, 20)

    ax.grid(None)
    ax.set_title(name)
    ax.set_ylabel('Step')
    ax.set_xlabel('Position (m)')

fig.tight_layout()
fig.savefig(folder + 'ElastoPlasticity.png')

# Plot figure
fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
for ax, pos in zip([ax1, ax2, ax3], [25, 50, 75]):
    ax.plot(strain[pos, :]*100, stress[pos, :])

    ax.set_ylabel('Stress (MPa)')
    ax.set_xlabel('Strain (%)')

fig.tight_layout()
fig.savefig(folder + 'TractionCurve.png')