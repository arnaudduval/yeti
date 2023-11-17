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

def PBiCGSTAB(Afun, Pfun, b, threshold=1e-5, nbIter=50):
    x = np.zeros(np.shape(b))
    r = b; normb = np.linalg.norm(r)
    resPCG = [1.0]
    if normb <=threshold: return
    rhat = r; p = r
    rsold = np.dot(r, rhat)

    for i in range(nbIter):
        ptilde = Pfun(p)
        Aptilde = Afun(ptilde)
        alpha = rsold/np.dot(Aptilde, rhat)
        s = r - alpha*Aptilde
        stilde = Pfun(s)
        Astilde = Afun(stilde)
        omega = np.dot(Astilde, s)/np.dot(Astilde, Astilde)
        x += alpha*ptilde + omega*stilde
        r = s - omega*Astilde

        resPCG.append(np.linalg.norm(r)/normb)
        if (resPCG[-1]<=threshold): break

        rsnew = np.dot(r, rhat)
        beta = (alpha/omega)*(rsnew/rsold)
        p = r + beta*(p - omega*Aptilde)
        rsold = rsnew

    return x, resPCG
    
# Time discretization
degree, nbel = 4, 512
knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel)
dscrt = GaussQuadrature(degree, knotvector, {})
dscrt.getQuadratureRulesInfo()
DenseBasis, DenseWeights = dscrt.getDenseQuadRules()

# Construct time matrices
Adv   = DenseWeights[1] @ np.diag(np.ones(dscrt.nbqp)) @ DenseBasis[-1].T; Adv  = Adv[1:, 1:]
Mass  = DenseWeights[0] @ np.diag(np.ones(dscrt.nbqp)) @ DenseBasis[0].T; Mass = Mass[1:, 1:]
Stiff = DenseWeights[-1] @ np.diag(np.ones(dscrt.nbqp)) @ DenseBasis[-1].T; Stiff = Stiff[1:, 1:]

import scipy.linalg as sclin
eigval, eigvec = sclin.eig(a=Stiff, b=Mass)
eigval = np.abs(eigval)
ones = np.ones(len(eigval))
D = np.kron(eigval, np.kron(ones, ones)) + np.kron(ones, np.kron(eigval, ones)) + np.kron(ones, np.kron(ones, eigval))
print(D.max(), D.min())

# PBiCGStab
nr = np.size(Adv, axis=0)
# invP = np.linalg.inv(Adv)
invP = np.linalg.inv(Mass)

for i in range(1):
    # coef = D[np.random.randint(0, len(D))]
    coef = 1000
    M = Adv + coef*Mass
    b = np.random.random(nr)

    Afun = lambda x: M @ x
    Pfun = lambda x: invP @ x
    Pfun = lambda x: 1/coef*invP @ x

    x, resPCG = PBiCGSTAB(Afun, Pfun, b)
    plt.semilogy(resPCG)
    plt.savefig("IterSolvers")


# # GMRES
# for m in range(1, 8):
#     x = GMRES(M, b, n_restarts=m, n_iter=15)
#     print('error:%.3e'%(np.linalg.norm(b - M @ x)/np.linalg.norm(b)))

# print('%===========%') 

# for m in range(1, 8):
#     x = GMRES(M, b, n_iter=m*15, n_restarts=1)
#     print('error:%.3e'%(np.linalg.norm(b - M @ x)/np.linalg.norm(b)))

