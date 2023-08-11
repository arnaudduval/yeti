"""
    This file contains FEM functions to 1D heat transfer problem
    We consider:
        - Isotropic material
        - There is no heat source, f = 0
        - Is a linear problem
"""

import os, numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import interpolate

full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

def compute_thermal_Fint_1D(C, K, T, dT):
	"Returns the internal heat force in transient heat"
	# Compute internal heat force 
	Fint = C @ dT + K @ T
	return Fint

def compute_tangent_thermal_matrix_1D(C, K, theta=0.5, dt=0.1):
	""" Computes tangent matrix in transient heat
		S = C + theta dt K
		K = int_Omega dB/dx Kprop dB/dx dx = int_[0, 1] J^-1 dB/dxi Kprop J^-1 dB/dxi detJ dxi.
		But in 1D: detJ times J^-1 get cancelled.
		C = int_Omega B Cprop B dx = int [0, 1] B Cprop det J B dxi
	"""
	# Compute tangent matrix 
	M = C + theta*dt*K
	return M

def solve_transient_heat_1D(C, K, time_list=None, dof=None, dod=None, Tinout=None, theta=0.5, threshold=1e-12, nbIterNL=20):
    " Solves transient heat problem in 1D. "

    ddGG = np.zeros(len(dod)) # d Temperature/ d time
    VVn0 = np.zeros(len(dof)+len(dod))

    # Compute initial velocity from boundry conditions (for i = 0)
    if np.shape(Tinout)[1] == 2:
        delta_t = time_list[1] - time_list[0]
        ddGG = (Tinout[dod, 1] - Tinout[dod, 0])/delta_t
    elif np.shape(Tinout)[1] >= 2:
        delta_t1 = time_list[1] - time_list[0]
        delta_t2 = time_list[2] - time_list[0]
        factor = delta_t2/delta_t1
        ddGG = 1.0/(delta_t1*(factor-factor**2))*(Tinout[dod, 2] - (factor**2)*Tinout[dod, 1] - (1 - factor**2)*Tinout[dod, 0])
    else: raise Warning('We need more than 2 steps')
    VVn0[dod] = ddGG

    for i in range(1, np.shape(Tinout)[1]):
        # Get delta time
        delta_t = time_list[i] - time_list[i-1]

        # Get values of last step
        TTn0 = Tinout[:, i-1]; TTn1 = np.copy(TTn0)

        # Get values of new step
        TTn1 = TTn0 + delta_t*(1-theta)*VVn0; TTn1[dod] = Tinout[dod, i]
        TTn1i0 = np.copy(TTn1); VVn1 = np.zeros(len(TTn1))
        VVn1[dod] = 1.0/theta*(1.0/delta_t*(Tinout[dod, i] - Tinout[dod, i-1]) - (1-theta)*VVn0[dod])

        for j in range(nbIterNL): # Newton-Raphson

            # Compute residue
            Fint = compute_thermal_Fint_1D(C, K, TTn1, VVn1)
            dF = - Fint[dof]
            prod1 = np.dot(dF, dF)
            relerror = np.sqrt(prod1)
            if relerror <= threshold: break

            # Compute tangent matrix
            MM = compute_tangent_thermal_matrix_1D(C, K, theta=theta, dt=delta_t)[np.ix_(dof, dof)]

            # Compute delta dT 
            ddVV = np.linalg.solve(MM, dF)

            # Update values
            VVn1[dof] += ddVV
            TTn1[dof] = TTn1i0[dof] + theta*delta_t*VVn1[dof]

        # Update values in output
        Tinout[:, i] = np.copy(TTn1)
        VVn0 = np.copy(VVn1)

    return 

# Elementary matrices
Cel = np.array([[2., 1.], [1., 2.]])
Kel = np.array([[1., -1.], [-1., 1.]])

# Properties
length       = 1.0
conductivity = 1.0
capacity     = 1.0
density      = 1.0

# Discretization
theta = 0.5
nbel  = 10; dh = 1/nbel
T     = 0.02
N     = 9
time_list = np.linspace(0, T, N)

# Assembly
C = np.zeros((nbel+1, nbel+1))
K = np.zeros((nbel+1, nbel+1))
for i in range(nbel):
    C[i:i+2, i:i+2] += Cel*capacity*density*dh/6.
    K[i:i+2, i:i+2] += Kel*conductivity/dh
    
temperature = np.zeros((nbel+1, len(time_list)))
temperature[0, :] = 0.0
temperature[-1,:] = 1.0
dod = [0, -1]; dof = np.arange(1, nbel, dtype=int)
solve_transient_heat_1D(C, K, time_list=time_list, dof=dof, dod=dod, Tinout=temperature, theta=1)

# ------------------
# Post-treatement
# ------------------
f = interpolate.interp2d(np.linspace(0, 1, nbel+1), time_list, temperature.T, kind='linear')
xnew = np.linspace(0, 1, 101)
ynew = np.linspace(0, T, 101)
znew = f(xnew, ynew)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,4))
ax.grid(False)
levels = np.array([-0.2]); levels = np.append(levels, np.linspace(0, 1, 9))
norm = mpl.colors.BoundaryNorm(levels, len(levels))
colors = list(plt.cm.Greys(np.linspace(0, 1, len(levels)-1)))
colors[0] = "red"
cmap = mpl.colors.ListedColormap(colors,"", len(colors))
im = ax.contourf(xnew, ynew, znew, norm=norm, cmap=cmap)
cbar = fig.colorbar(im)
cbar.set_label('Temperature (K)')

ax.set_ylabel('Time (s)')
ax.set_xlabel('Position (m)')
fig.tight_layout()
fig.savefig(folder + 'TransientHeat_1D.png')