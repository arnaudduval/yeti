# Python libraries
import numpy as np
from matplotlib import pyplot as plt

# My libraries
from lib.base_functions import (eval_basis_python,
                                iga_find_positions_weights,
                                create_knotvector
)
from lib.D1transient import *

def conductivity(T):
    # K = 1 + 2*np.exp(-abs(T))
    K = 0.1*np.ones(np.shape(T))
    return K

def capacity(T):
    C = 1 + 0.2*np.exp(-0.1*abs(T))
    # C = 1*np.ones(np.shape(T))
    return C

def source(qp):
    f = 1*np.sin(np.pi*qp)
    return f

# Define geometry
alpha, JJ = 0.5, 1
degree, nbel = 5, 8
nb_ctrlpts = degree + nbel
ctrlpts = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)
qp_cgg, W = iga_find_positions_weights(degree, knotvector)
DB = eval_basis_python(degree, knotvector, qp_cgg)
Fcoefs = source(qp_cgg)
properties = [JJ, conductivity, capacity, alpha]

# Define boundary and initial conditions
N = 100
n = 30
t = np.linspace(0, 1, N)
dof = np.arange(1, nb_ctrlpts-1, dtype=int)
Ffinal = compute_source_1D(JJ, DB, W, Fcoefs)
Ffinal = np.atleast_2d(Ffinal).reshape(-1, 1)
Fext = np.zeros((len(Ffinal), len(t)+n))
Fext[:,:len(t)] = np.kron(Ffinal, t)
for i in range(len(t), len(t)+n):
    Fext[:,i] = Fext[:,len(t)-1]

# Solve problem
time_list = np.zeros(N+n)
time_list[:N] = t
time_list[N:] = [1 + 0.2*(i+1) for i in range(n)]
TT = solve_transient_1D(properties, DB=DB, W =W, Fext=Fext, time_list=time_list, dof=dof)

# Post-treatement
xi = np.linspace(0, 1, 101)
B0, B1 = eval_basis_python(degree, knotvector, xi)
TT_it = B0.T*TT

# Plot figure
names = ['Temperature field']
XX, TT = np.meshgrid(xi*JJ, time_list)
fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(14,4))
for ax, variable, name in zip([ax1], [TT_it], names):
    # Plot
    im = ax.contourf(XX, TT, variable.T, 100)

    # Properties
    plt.colorbar(im)
    ax.set_title(name, fontsize=14)
    ax.set_ylabel('Time (s)', fontsize=12)
    ax.set_xlabel('Position (m)', fontsize=12)
    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)

plt.tight_layout()
plt.savefig('Transient.png')

