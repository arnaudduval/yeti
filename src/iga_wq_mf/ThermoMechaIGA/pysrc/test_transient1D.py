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
    K = 1 + 20*np.exp(-T)
    return K

def capacity(T):
    C = 2 + 200*np.exp(-T)
    return C

def source(qp):
    f = 100*np.sin(np.pi*qp)
    return f

# Define geometry
alpha, JJ = 0.5, 1
degree, nbel = 5, 32
nb_ctrlpts = degree + nbel
ctrlpts = np.linspace(0, 1, nb_ctrlpts)
knotvector = create_knotvector(degree, nbel)
qp_cgg, W = iga_find_positions_weights(degree, knotvector)
DB = eval_basis_python(degree, knotvector, qp_cgg)
Fcoefs = source(qp_cgg)
properties = [JJ, conductivity, capacity, alpha]

# Define boundary and initial conditions
N = 500
t = np.linspace(0, 1, N)
dof = np.arange(1, nb_ctrlpts, dtype=int)
Ffinal = compute_source_1D(JJ, DB, W, Fcoefs)
Ffinal = np.atleast_2d(Ffinal).reshape(-1, 1)
Fext = np.kron(Ffinal, t)

# Solve problem
TT = solve_transient_1D(properties, DB=DB, W =W, Fext=Fext, time_list=t, dof=dof)

# Post-treatement
xi = np.linspace(0, 1, 101)
B0, B1 = eval_basis_python(degree, knotvector, xi)
TT_it = B0.T*TT

# Plot figure
names = ['Temperature field']
XX, TT = np.meshgrid(xi*JJ, np.arange(N))
fig, [ax1] = plt.subplots(nrows=1, ncols=1, figsize=(14,4))
for ax, variable, name in zip([ax1], [TT_it], names):
    # Plot
    ax.contourf(XX, TT, variable.T, 20)

    # Properties
    ax.set_title(name, fontsize=14)
    ax.set_ylabel('Step', fontsize=12)
    ax.set_xlabel('Position (m)', fontsize=12)
    ax.tick_params(axis='x', labelsize=11)
    ax.tick_params(axis='y', labelsize=11)

plt.tight_layout()
plt.savefig('Transient.png')

