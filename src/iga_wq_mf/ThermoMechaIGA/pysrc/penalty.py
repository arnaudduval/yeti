"""
Analysis of penalty method
"""

# Python libraries
from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt

from lib.base_functions import *

def array2csr(indi, indj, B, W, isfortran=True):

    # Get information 
    nr = len(indi) -1
    nc = np.max(indj)
    indit, indjt = deepcopy(indi), deepcopy(indj)
    if isfortran: indit -=1; indjt -= 1

    # Create basis and weights from fortran
    B0 = sp.csr_matrix((B[:, 0], indjt, indit), shape=(nr, nc))
    B1 = sp.csr_matrix((B[:, 1], indjt, indit), shape=(nr, nc))
    W00 = sp.csr_matrix((W[:, 0], indjt, indit), shape=(nr, nc))
    W01 = sp.csr_matrix((W[:, 1], indjt, indit), shape=(nr, nc))
    W10 = sp.csr_matrix((W[:, 2], indjt, indit), shape=(nr, nc))
    W11 = sp.csr_matrix((W[:, 3], indjt, indit), shape=(nr, nc))

    return [B0, B1, W00, W01, W10, W11]

# ========================
# INPUTS
# ========================
"""
We want to prove accuracy of penalty method. For that, 2 cases are going to be solve 
"""
# Create basis and weights
degree = 4
nbel = 3
nb_ctrlpts = degree + nbel
knotvector = create_knotvector(degree, nbel)
_, qpos, B, W, indi, indj = wq_find_basis_weights_fortran(degree, knotvector)
dof1D = np.arange(1, nb_ctrlpts-1)

# Eigen decomposition
data = [B[:, 0], B[:, 1], W[:, 0], W[:, -1]]
DD, UU = eigen_decomposition(indi, indj, data, t_robin=[1, 1])

# Get some matrices
B0, B1, W00, _, _, W11 = array2csr(indi, indj, B, W)
M = W00 @ B0.T
K = W11 @ B1.T

# # ========================
# # CASE 1-D
# # ========================
# def source1D(P):
#     f = np.sin(np.pi*P)
#     return f

# # Define quadrature points:
# QP1D = qpos

# # Define source vector
# F = W00 @ source1D(QP1D)

# # Solve with Dirichlet lift
# K2solve = K[dof1D, :][:, dof1D].todense()
# F2solve = F[dof1D]
# T1 = np.zeros(np.shape(F))
# T1[dof1D] = np.linalg.solve(K2solve, F2solve)
# print(T1)

# # Solve with penalty
# T2 = UU @ ((UU.T @ F)/DD) 
# print(T2)

# ========================
# CASE 2-D
# ========================
def source2D(P):
    x = P[:, 0]
    y = P[:, 1]
    f = np.sin(np.pi*x)*np.sin(np.pi*y)
    return f

# Define quadrature points:
QP2D = np.zeros((len(qpos)**2,2))
for j, y in enumerate(qpos):
    for i, x in enumerate(qpos):
        k = i + j*len(qpos)
        QP2D[k, 0] = x
        QP2D[k, 1] = y

# Define source vector
F = sp.kron(W00, W00) @ source2D(QP2D)

dof2D = np.zeros(len(dof1D)**2, dtype=int)
for j, y in enumerate(dof1D):
    for i, x in enumerate(dof1D):
        k = i + j*len(dof1D)
        dof2D[k] = x + y*nb_ctrlpts

# Solve with Dirichlet lift
K2D = sp.kron(K, M) + sp.kron(M, K)
K2solve = K2D.todense()[np.ix_(dof2D, dof2D)]
F2solve = F[dof2D]
T1 = np.zeros(np.shape(F))
T1[dof2D] = np.linalg.solve(K2solve, F2solve)
Tinterp1 = sp.kron(B0, B0).T @ T1

# Solve with penalty
Deigen = np.kron(DD, np.ones(len(DD))) + np.kron(np.ones(len(DD)), DD) 
Ueigen = np.kron(UU, UU)
T2 = Ueigen @ ((Ueigen.T @ F)/Deigen) 
Tinterp2 = sp.kron(B0, B0).T @ T2

# Plots
XX, YY = np.meshgrid(qpos, qpos)
Tf1 = np.reshape(Tinterp1, (len(qpos), len(qpos)))
Tf2 = np.reshape(Tinterp2, (len(qpos), len(qpos)))
Tf3 = abs(Tf1 - Tf2)
fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(16, 5))
ax1.contourf(XX, YY, Tf1, 200)
ax1.set_title('max = %.3e, min = %.3e' %(Tf1.max(), Tf1.min()))
ax2.contourf(XX, YY, Tf2, 200)
ax2.set_title('max = %.3e, min = %.3e' %(Tf2.max(), Tf2.min()))
ax3.contourf(XX, YY, Tf3, 200)
ax3.set_title('max = %.3e, min = %.3e' %(Tf3.max(), Tf3.min()))
plt.tight_layout()
plt.savefig('Penalty.png')
