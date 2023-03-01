from lib.__init__ import *
from lib.base_functions import (eval_basis_python,
								iga_find_positions_weights,
								create_knotvector, 
                                wq_find_basis_weights_fortran
)
from lib.D1viscoplasticity import *
from lib.physics import forceVol

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/examples/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
mechaprop     = {'young':200e3, 'sigma_Y0':506, 'K':2e4, 'exp':0.5}
mechabehavior = MechaBehavior('swift', kwargs=mechaprop)
length        = 1.0
degree        = 4

# Get basis and weights in IgA 
nbel1       = 51
nb_ctrlpts1 = degree + nbel1
ctrlpts1    = np.linspace(0, 1, nb_ctrlpts1)
knotvector1 = create_knotvector(degree, nbel1)
qp_cgg, weight_cgg = iga_find_positions_weights(degree, knotvector1)
basis_cgg          = eval_basis_python(degree, knotvector1, qp_cgg)
properties_iga     = [length, len(qp_cgg)]

# Get basis and weights in WQ 
method      = 1
# nbel2       = int(nbel1*(degree+1)/2-degree)
nbel2       = nbel1
nb_ctrlpts2 = degree + nbel2
ctrlpts2    = np.linspace(0, 1, nb_ctrlpts2)
knotvector2 = create_knotvector(degree, nbel2)
qp_wq, B, W, indi, indj = wq_find_basis_weights_fortran(degree, knotvector2, method=method)[1:]
indi -= 1; indj -= 1; DB = []; DW = []
for i in range(2): DB.append(sp.csr_matrix((B[:, i], indj, indi), shape=(nb_ctrlpts2, len(qp_wq))))
for i in range(4): DW.append(sp.csr_matrix((W[:, i], indj, indi), shape=(nb_ctrlpts2, len(qp_wq))))
properties_wq  = [length, len(qp_wq)]

# Define boundaries conditions
nbSteps = 101
dof1        = np.arange(1, nb_ctrlpts1, dtype=int)
Fext        = np.zeros((nb_ctrlpts1, nbSteps))
Fext[:, -1] = compute_IGA_Fvol_1D(basis_cgg, weight_cgg, forceVol(qp_cgg))
for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

# Solve using IGA
disp_iga, strain_iga, stress_iga, plastic_iga, Cep_iga = solve_IGA_plasticity_1D(properties_iga, mechabehavior, DB=basis_cgg, W=weight_cgg, Fext=Fext, dof=dof1)
strain_cp   = interpolate_IGA_CP_1D(basis_cgg, weight_cgg, strain_iga)
plastic_cp  = interpolate_IGA_CP_1D(basis_cgg, weight_cgg, plastic_iga)
stress_cp 	= interpolate_IGA_CP_1D(basis_cgg, weight_cgg, stress_iga)
plot_results(degree, knotvector1, length, disp_iga, plastic_cp,
                stress_cp, folder=folder, method='IGA')

fig, [ax0, ax1] = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
ax0.plot(qp_cgg, stress_iga[:, 51])
ax0.set_ylabel('Stress (MPa)')
ax1.plot(qp_cgg, Cep_iga[:, 51])
ax1.set_ylabel('Tangent modulus (MPa)')
for ax in [ax0, ax1]:
    ax.set_ylim(bottom=0.0)
    ax.set_xlim(left=0.0, right=1.0)
    ax.set_xlabel('Quadrature point position')
fig.tight_layout()
fig.savefig(folder+'data_step51')

# Define boundaries conditions
dof2        = np.arange(1, nb_ctrlpts2, dtype=int)
Fext        = np.zeros((nb_ctrlpts2, nbSteps))
Fext[:, -1] = compute_WQ_Fvol_1D(DW, forceVol(qp_wq))
for i in range(1, nbSteps-1): Fext[:, i] = i/(nbSteps-1)*Fext[:, -1]

# Solve using WQ
disp_wq, strain_wq, stress_wq, plastic_wq, Cep_wq = solve_WQ_plasticity_1D(properties_wq, mechabehavior, DB=DB, DW=DW, Fext=Fext, dof=dof2)
strain_cp   = interpolate_WQ_CP_1D(DB, DW, strain_wq)
plastic_cp  = interpolate_WQ_CP_1D(DB, DW, plastic_wq)
stress_cp   = interpolate_WQ_CP_1D(DB, DW, stress_wq)
plot_results(degree, knotvector2, length, disp_wq, plastic_cp,
                stress_cp, folder=folder, method='WQ')

# ------------------
# RESULTS
# ------------------
relerror = np.array([])
for i in range(1, nbSteps):
    tmp = relativeNormL2(degree, knotvector1, knotvector2, disp_iga[:, i], disp_wq[:, i])
    relerror = np.append(relerror, tmp)

fig, ax = plt.subplots(nrows=1, ncols=1)
ax.semilogy(relerror*100)
ax.set_ylim(top=1e0, bottom=1e-16)
ax.set_ylabel('L2 Relative error (\%)')
ax.set_xlabel('Step')
fig.tight_layout()
fig.savefig(folder + 'Relative_error_distance.png')

if (len(qp_cgg)%2==1):
    error    = np.abs(stress_iga[int((len(qp_cgg)-1)/2), 1:] - stress_wq[int((len(qp_wq)-1)/2), 1:])
    fig, ax  = plt.subplots(nrows=1, ncols=1)
    ax.semilogy(np.abs(error)/np.abs(stress_iga[int((len(qp_cgg)-1)/2), 1:])*100)
    ax.set_ylim(bottom=1e-12, top=1e0)
    ax.set_ylabel('Relative error (\%)')
    ax.set_xlabel('Step')
    fig.tight_layout()
    fig.savefig(folder + 'Relative_error_stress.png')
else: 
    print('Even number of quadrature points')