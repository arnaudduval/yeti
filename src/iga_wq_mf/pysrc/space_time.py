"""
In this file, we test space-time Galerkin method in 1D
Hypotesis: 
	- Isotropic material
	- There is no heat source
    - Conditions: T(x, 0) = 0, T(0, t) = 0, T(1, t) = 1
"""

from lib.__init__ import *
from lib.base_functions import (eval_basis_python,
								iga_find_positions_weights,
								create_knotvector, 
)
from lib.create_geomdl import geomdlModel
from lib.create_model import thermoMechaModel

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/examples/'
if not os.path.isdir(folder): os.mkdir(folder)

# Material and geometry
rho, Cap, Cond = 1.0, 1.0, 1.0
length 		   = 1.0
Tspan          = 0.02

# Space discretization
p_sp, nbel    = 1, 10; nbctrlpts_sp = p_sp + nbel
knotvector_sp = create_knotvector(p_sp, nbel)
qp_sp, wgt_sp = iga_find_positions_weights(p_sp, knotvector_sp)
basis_sp      = eval_basis_python(p_sp, knotvector_sp, qp_sp)

# Construct space matrices
Mcfs     = wgt_sp * rho * Cap * length
Mass_sp  = basis_sp[0] @ np.diag(Mcfs) @ basis_sp[0].T

Kcfs     = wgt_sp * Cond * 1.0/length
Stiff_sp = basis_sp[1] @ np.diag(Kcfs) @ basis_sp[1].T 

# Time discretization
p_t, nsteps  = 1, 9; nbctrlpts_t = p_t + nsteps
knotvector_t = create_knotvector(p_t, nsteps)
qp_t, wgt_t  = iga_find_positions_weights(p_t, knotvector_t)
basis_t      = eval_basis_python(p_t, knotvector_t, qp_t)

# Construct time matrices
Adv_time  = basis_t[0] @ np.diag(wgt_t) @ basis_t[1].T
Mass_time = basis_t[0] @ np.diag(wgt_t*Tspan) @ basis_t[0].T

# Get free nodes. To reuse algorithms we use the class thermomodel
all_ctrlpts = set(np.arange(nbctrlpts_sp*nbctrlpts_t))
all_dof = []
for i in range(1, nbctrlpts_sp - 1):
    for j in range(1, nbctrlpts_t):
        k = i + j*nbctrlpts_sp
        all_dof.append(k)
all_dof = set(all_dof); all_dod = all_ctrlpts.difference(all_dof)
all_dof = list(all_dof); all_dod = list(all_dod)

spe_dod = []
for i in [nbctrlpts_sp-1]:
    for j in range(nbctrlpts_t):
        k = i + j*nbctrlpts_sp
        spe_dod.append(k)

# Define space-time temperature
u = np.zeros(len(all_ctrlpts))
u[spe_dod] = 1.0

# Solve problem
A = np.kron(Adv_time, Mass_sp) + np.kron(Mass_time, Stiff_sp)
Ann = A[all_dof, :][:, all_dof]
b = - A @ u; bn = b[all_dof]
u_n = np.linalg.solve(Ann, bn)
u[all_dof] = u_n

# Interpolate solution
qp_sp, qp_t = np.linspace(0, 1, 101), np.linspace(0, 1, 101)
basis_sp  = eval_basis_python(p_sp, knotvector_sp, qp_sp) 
basis_t   = eval_basis_python(p_t, knotvector_t, qp_t)
spt_basis = sp.kron(basis_t[0], basis_sp[0])
u_interp  = spt_basis.T @ u
UU = np.reshape(u_interp, (len(qp_sp), len(qp_t)), order='F')
print(UU.min())

# Plot
XX, TT  = np.meshgrid(qp_sp*length, qp_t*Tspan)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,4))
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,4))
levels  = np.array([-0.2]); levels = np.append(levels, np.linspace(0, 1, 9))
norm    = mpl.colors.BoundaryNorm(levels, len(levels))
colors  = list(plt.cm.Greys(np.linspace(0, 1, len(levels)-1))); colors[0] = "red"
cmap = mpl.colors.ListedColormap(colors,"", len(colors))
im   = ax.contourf(XX, TT, UU.T, norm=norm, cmap=cmap)
cbar = plt.colorbar(im)
cbar.set_label('Temperature (K)')

ax.grid(False)
ax.set_ylabel('Time (s)')
ax.set_xlabel('Position (m)')
fig.tight_layout()
fig.savefig(folder + 'space_time.png')