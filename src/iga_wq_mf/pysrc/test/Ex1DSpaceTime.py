from lib.__init__ import *
from lib.lib_base import createKnotVector, evalDersBasisPy
from lib.lib_quadrules import GaussQuadrature

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/d1heat/'
if not os.path.isdir(folder): os.mkdir(folder)

# Material and geometry
rho, Cap, Cond = 1.0, 1.0, 1.0
length 		   = 1.0
Tspan          = 0.02

# Space discretization
degree_sp, nbel = 2, 10
knotvector_sp   = createKnotVector(degree_sp, nbel)
spDscrt         = GaussQuadrature(degree_sp, knotvector_sp)
spDscrt.getQuadratureRulesInfo()
spDenseBasis, spDenseWeights = spDscrt.getDenseQuadRules()

# Construct space matrices
Mcfs     = rho * Cap * length * np.ones(spDscrt.nbqp)
Mass_sp  = spDenseWeights[0] @ np.diag(Mcfs) @ spDenseBasis[0].T

Kcfs     = Cond/length * np.ones(spDscrt.nbqp)
Stiff_sp = spDenseWeights[-1] @ np.diag(Kcfs) @ spDenseBasis[1].T 

# Time discretization
degree_t, nsteps  = 2, 9
knotvector_t      = createKnotVector(degree_t, nsteps)
tmDscrt           = GaussQuadrature(degree_t, knotvector_t)
tmDscrt.getQuadratureRulesInfo()
tmDenseBasis, tmDenseWeights = tmDscrt.getDenseQuadRules()

# Construct time matrices
Adv_time  = tmDenseWeights[0] @ np.diag(np.ones(tmDscrt.nbqp)) @ tmDenseBasis[1].T
Mass_time = tmDenseWeights[0] @ np.diag(Tspan*np.ones(tmDscrt.nbqp)) @ tmDenseBasis[0].T

# Get free nodes. To reuse algorithms we use the class thermomodel
all_ctrlpts = set(np.arange(spDscrt.nbctrlpts*tmDscrt.nbctrlpts))
all_dof = []
for i in range(1, spDscrt.nbctrlpts - 1):
    for j in range(1, tmDscrt.nbctrlpts):
        k = i + j*spDscrt.nbctrlpts
        all_dof.append(k)
all_dof = set(all_dof); all_dod = all_ctrlpts.difference(all_dof)
all_dof = list(all_dof); all_dod = list(all_dod)

spe_dod = []
for i in [spDscrt.nbctrlpts-1]:
    for j in range(tmDscrt.nbctrlpts):
        k = i + j*spDscrt.nbctrlpts
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

# ------------------
# Post-treatement
# ------------------
# Interpolate solution
qp_sp, qp_t = np.linspace(0, 1, 101), np.linspace(0, 1, 101)
basis_sp  = evalDersBasisPy(degree_sp, knotvector_sp, qp_sp) 
basis_t   = evalDersBasisPy(degree_t, knotvector_t, qp_t)
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