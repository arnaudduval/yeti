from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity, evalDersBasisPy
from pysrc.lib.lib_quadrules import GaussQuadrature

def GMRES(A, b, x0=None, n_iter=15, n_restarts=1):
    
    if x0 is None: x = np.zeros(len(b))
    else: x = np.copy(x0)
    n = x.shape[0]
    H = np.zeros((n_iter + 1, n_iter))
    V = np.zeros((n_iter + 1, n))
    
    for k in range(n_restarts):
        r = b - A @ x
        beta = np.linalg.norm(r)
        V[0] = r / beta

        for j in range(n_iter):
            w = A @ V[j]
            for i in range(j + 1):
                H[i, j] = np.dot(w, V[i])
                w -= H[i, j] * V[i]
            H[j + 1, j] = np.linalg.norm(w)
            V[j + 1] = w / H[j + 1, j]

        e1 = np.zeros(n_iter + 1); e1[0] = beta
        y = np.linalg.lstsq(H, e1)[0]
        xk = x + V[:-1].T @ y
        # ADD STOP CRITERIA
        x = np.copy(xk)

    return xk

# Time discretization
degree, nbel = 4, 128
knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel)
dscrt = GaussQuadrature(degree, knotvector, {})
dscrt.getQuadratureRulesInfo()
DenseBasis, DenseWeights = dscrt.getDenseQuadRules()

# Construct time matrices
Adv  = DenseWeights[0] @ np.diag(np.ones(dscrt.nbqp)) @ DenseBasis[1].T; Adv  = Adv[1:, 1:]
Mass = DenseWeights[0] @ np.diag(np.ones(dscrt.nbqp)) @ DenseBasis[0].T; Mass = Mass[1:, 1:]
Adv_sym = 0.5*(Adv + Adv.T); Adv_asym = Adv - Adv_sym

# Eigen decomposition
Meigval, Meigvec = np.linalg.eig(Mass)
# print(Meigvec.T@Adv_sym@Meigvec)
# print(Meigvec.T@Adv_asym@Meigvec)
# Just to show that U^T Adv U is not a quase skew matrix 

import random
Z = Meigvec.T@Adv@Meigvec
D = random.randint(1, 100)*np.diag(Meigval)
M = D + Z

b = np.random.random(len(Meigval))
for m in range(1, 8):
    x = GMRES(M, b, n_restarts=m, n_iter=15)
    print('error:%.3e'%(np.linalg.norm(b - M @ x)/np.linalg.norm(b)))

print('%===========%') 

for m in range(1, 8):
    x = GMRES(M, b, n_iter=m*15, n_restarts=1)
    print('error:%.3e'%(np.linalg.norm(b - M @ x)/np.linalg.norm(b)))

