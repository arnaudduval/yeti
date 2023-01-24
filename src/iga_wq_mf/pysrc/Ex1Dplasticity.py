from lib.__init__ import *
from lib.base_functions import (eval_basis_python,
								iga_find_positions_weights,
								create_knotvector, 
                                wq_find_basis_weights_fortran
)
from lib.D1Plasticity import *
from lib.physics import forceVol

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/examples/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
E, H, sigma_Y, beta, JJ = 200e3, 50e3, 100, 0.5, 1.0
degree, nbel            = 10, 64

method     = 1
nb_ctrlpts = degree + nbel
ctrlpts    = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)

# Get basis and weights in IgA 
qp_cgg, weight_cgg = iga_find_positions_weights(degree, knotvector)
basis_cgg          = eval_basis_python(degree, knotvector, qp_cgg)
properties_iga     = [JJ, E, H, beta, sigma_Y, len(qp_cgg)]

# Get basis and weights in WQ 
qp_wq, B, W, indi, indj = wq_find_basis_weights_fortran(degree, knotvector, method=method)[1:]
indi -= 1; indj -= 1; DB = []; DW = []
for i in range(2): DB.append(sp.csr_matrix((B[:, i], indj, indi), shape=(nb_ctrlpts, len(qp_wq))))
for i in range(4): DW.append(sp.csr_matrix((W[:, i], indj, indi), shape=(nb_ctrlpts, len(qp_wq))))
properties_wq  = [JJ, E, H, beta, sigma_Y, len(qp_wq)]

# Define boundaries conditions
nbSteps = 101
time_list   = np.linspace(0, 1, nbSteps)
dof         = np.arange(1, nb_ctrlpts, dtype=int)
Fext        = np.zeros((nb_ctrlpts, nbSteps))
Fext[:, -1] = compute_IGA_Fvol_1D(basis_cgg, weight_cgg, forceVol(qp_cgg))
Fext[-1,-1] += 400
for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

# Solve using IGA
disp_iga, strain, stress, plastic = solve_IGA_plasticity_1D(properties_iga, DB=basis_cgg, W=weight_cgg, Fext=Fext, dof=dof)
strain_cp   = interpolate_IGA_CP_1D(basis_cgg, weight_cgg, strain)
plastic_cp  = interpolate_IGA_CP_1D(basis_cgg, weight_cgg, plastic)
stress_cp 	= interpolate_IGA_CP_1D(basis_cgg, weight_cgg, stress)
plot_results(degree, knotvector, JJ, disp_iga, plastic_cp, 
                stress_cp, folder=folder, method='IGA', extension='.pdf')

# Solve using WQ
disp_wq, strain, stress, plastic = solve_WQ_plasticity_1D(properties_wq, DB=DB, DW=DW, Fext=Fext, dof=dof)
strain_cp   = interpolate_WQ_CP_1D(DB, DW, strain)
plastic_cp  = interpolate_WQ_CP_1D(DB, DW, plastic)
stress_cp   = interpolate_WQ_CP_1D(DB, DW, stress)
plot_results(degree, knotvector, JJ, disp_wq, plastic_cp, 
                stress_cp, folder=folder, method='WQ')

# ------------------
# RESULTS
# ------------------
error    = disp_iga[:, 1:] - disp_wq[:, 1:]
relerror = np.linalg.norm(error, np.inf, axis=0)/np.linalg.norm(disp_iga[:, 1:], np.inf, axis=0)
fig, ax  = plt.subplots(nrows=1, ncols=1)
ax.semilogy(relerror*100)
ax.set_ylim(bottom=1e-12, top=1e0)
ax.set_ylabel('Relative error (\%)')
ax.set_xlabel('Step')
fig.tight_layout()
fig.savefig(folder + 'Relative_error.png')