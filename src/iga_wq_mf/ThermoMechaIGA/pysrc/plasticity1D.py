# Python libraries
import numpy as np
from matplotlib import pyplot as plt

# My libraries
from lib.base_functions import (eval_basis_python,
                                iga_find_positions_weights,
                                create_knotvector
)
from lib.D1viscoplasticity import interpolate_controlPoints, solve_plasticity, compute_strain

# Define mechanical properties
E, H, sigma_Y, beta, JJ = 200e3, 25e3, 250, 0.5, 10

# Define geometry
degree, nbel = 5, 32
nb_ctrlpts = degree + nbel
ctrlpts = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)
qp_cgg, W = iga_find_positions_weights(degree, knotvector)
DB = eval_basis_python(degree, knotvector, qp_cgg)
nnz = len(qp_cgg)
properties = [JJ, E, H, beta, sigma_Y, nnz]

# Define boundary and initial conditions
N = 800
t = np.linspace(0, 1, N)
dof = np.arange(1, nb_ctrlpts, dtype=int)
Fext = np.zeros((nb_ctrlpts, N))
Fext[-1, :] = 400*t

# Solve problem
disp, epn, sigma = solve_plasticity(properties, DB=DB, W=W, Fext=Fext, dof=dof)
epn_cp = interpolate_controlPoints(DB, W, epn)
sigma_cp = interpolate_controlPoints(DB, W, sigma)

# Post-treatement
xi = np.linspace(0, 1, 101)
B0, B1 = eval_basis_python(degree, knotvector, xi)
disp_it = B0.T*disp
ep_it = compute_strain(JJ, [B0, B1], disp) 
epn_it = B0.T*epn_cp
sigma_it = B0.T*sigma_cp

# Plot figure
names = ['Displacement field', 'Plastic strain field', 'Stress field']
XX, TT = np.meshgrid(xi*JJ, np.arange(N))
fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
for ax, variable, name in zip([ax1, ax2, ax3], [disp_it, epn_it, sigma_it], names):
    # Plot
    ax.contourf(XX, TT, variable.T, 20)

    # Properties
    ax.set_title(name, fontsize=14)
    ax.set_ylabel('Step', fontsize=12)
    ax.set_xlabel('Position (m)', fontsize=12)
    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)

plt.tight_layout()
plt.savefig('Plasticity.png')


# Plot figure
fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
for ax, pos in zip([ax1, ax2, ax3], [25, 50, 75]):

    # Plot
    sig = sigma_it[pos, :]
    eps = ep_it[pos, :]*100
    ax.plot(eps, sig)

    # Properties
    ax.grid()
    ax.set_title(name, fontsize=14)
    ax.set_ylabel('Stress (MPa)', fontsize=12)
    ax.set_xlabel('Strain (%)', fontsize=12)
    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)

plt.tight_layout()
plt.savefig('Curve.png')