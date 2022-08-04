# Python libraries
import numpy as np
from matplotlib import pyplot as plt

# My libraries
from lib.base_functions import (eval_basis_python,
                                iga_find_positions_weights,
                                create_knotvector
)
from lib.D1viscoplasticity import interpolate_CP, newton_raphson

# Define mechanical properties
E, K, H, sigma_Y = 210e3, 100e3, 50e3, 80
rho, g = 7.8e-6, 9.81e3
L = 100

# Define geometry
degree, nbel = 5, 16
nb_ctrlpts = degree + nbel
ctrlpts = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)
qp_cgg, W = iga_find_positions_weights(degree, knotvector)
DB = eval_basis_python(degree, knotvector, qp_cgg)
nnz = len(qp_cgg)

# Define boundary and initial conditions
N = 1000
t = np.linspace(0, 10, N)
b = L*rho*g*np.zeros(nnz)
dof = np.arange(1, nb_ctrlpts, dtype=int)
displacement = np.zeros((nb_ctrlpts, N))
stress = np.zeros((nb_ctrlpts, N))
vsigma_ext = np.zeros((nb_ctrlpts, N))
vsigma_ext[-1, :] = 300*np.sin(2*np.pi/10*t) 

# Initialize
eps_pn0, alpha_n0, q_n0 = np.zeros(nnz), np.zeros(nnz), np.zeros(nnz)

# Solve 
for i in range(N-1):
    sigma_ext = vsigma_ext[:, i+1] 
    d_n0 = displacement[:, i]
    d_n1, sigma_n1, eps_pn1, alpha_n1, q_n1 = newton_raphson(L, E, K, H, sigma_Y, DB, W, d_n0, b, sigma_ext, dof, eps_pn0, alpha_n0, q_n0)
    displacement[:, i+1] = d_n1
    eps_pn0, alpha_n0, q_n0 = eps_pn1, alpha_n1, q_n1

    # Interpolate 
    sigma_CP = interpolate_CP(DB, W, sigma_n1)
    stress[:, i+1] = sigma_CP

# Post-treatement
xi = np.linspace(0, 1, 101)
B0, B1 = eval_basis_python(degree, knotvector, xi)
disp_it = B0.T*displacement
eps_it = B1.T*displacement/L
sigma_it = B0.T*stress

# Plot figure
names = ['Displacement', 'Strain field', 'Stress field']
XX, TT = np.meshgrid(xi*L, t)
fig, [ax1, ax2, ax3] = plt.subplots(nrows=1, ncols=3, figsize=(14,4))
for ax, variable, name in zip([ax1, ax2, ax3], [disp_it, eps_it, sigma_it], names):
    # Plot
    print(name + ' residual : %.5f' %variable[:,-1].max())
    print(name + ' max value : %.5f' %variable.max())
    variable /= variable.max()
    ax.contourf(XX, TT, variable.T, 100)

    # Properties
    ax.set_title(name, fontsize=14)
    ax.set_ylabel('Time (s)', fontsize=12)
    ax.set_xlabel('Position (m)', fontsize=12)
    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)

fig.tight_layout()
plt.savefig('Plasticity.png')