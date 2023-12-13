from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity, evalDersBasisPy
from pysrc.lib.lib_quadrules import GaussQuadrature

def GMRES(A, b, invP=None, x0=None, n_iter=15, n_restarts=1, tol=1e-8):
    
    if x0 is None: x = np.zeros(len(b))
    else: x = np.copy(x0)
    if invP is None: invP = np.eye(len(b))

    H = np.zeros((n_iter + 1, n_iter))
    V = np.zeros((n_iter + 1, len(b)))
    Z = np.zeros((n_iter + 1, len(b)))
    beta = np.zeros(n_restarts)
    
    for m in range(n_restarts):
        r = b - A @ x
        beta[m] = np.linalg.norm(r)
        if beta[m] <= tol*beta[0]: break
        V[0] = r / beta[m]
        e1 = np.zeros(n_iter + 1); e1[0] = beta[m]

        for k in range(n_iter):
            Z[k] = invP @ V[k]
            w = A @ Z[k]

            for j in range(k+1):
                H[j, k] = np.dot(w, V[j])
                w -= H[j, k] * V[j]
            H[j+1, j] = np.linalg.norm(w)
            if H[j+1, j] != 0: V[k+1] = w/H[j+1, j]
            y = np.linalg.lstsq(H[:k+2, :k+1], e1[:k+2])[0]
            rho = np.linalg.norm(H[:k+2, :k+1] @ y - e1[:k+2])
            if rho <= tol*beta[0]: break

        xk = x + Z[:k+1].T @ y
        x = np.copy(xk)
    print('outer: %d, inner: %d, total: %d' %(m+1, k+1, m*n_iter+k+1))

    return x

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
degree, nbel = 4, 128
knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel)
dscrt = GaussQuadrature(degree, knotvector, {})
dscrt.getQuadratureRulesInfo()
DenseBasis, DenseWeights = dscrt.getDenseQuadRules()

# Construct time matrices
Adv   = DenseWeights[1] @ np.diag(np.ones(dscrt.nbqp)) @ DenseBasis[-1].T; Adv  = Adv[1:, 1:]
Mass  = DenseWeights[0] @ np.diag(np.ones(dscrt.nbqp)) @ DenseBasis[0].T; Mass = Mass[1:, 1:]
Stiff = DenseWeights[-1] @ np.diag(np.ones(dscrt.nbqp)) @ DenseBasis[-1].T; Stiff = Stiff[1:, 1:]

# import scipy.linalg as sclin
# eigval, eigvec = sclin.eig(a=Stiff, b=Mass)
# eigval = np.abs(eigval)
# ones = np.ones(len(eigval))
# D = np.kron(eigval, np.kron(ones, ones)) + np.kron(ones, np.kron(eigval, ones)) + np.kron(ones, np.kron(ones, eigval))
# print(D.max(), D.min())

# # PBiCGStab
# nr = np.size(Adv, axis=0)
# # invP = np.linalg.inv(Adv)
# invP = np.linalg.inv(Mass)

# for i in range(1):
#     # coef = D[np.random.randint(0, len(D))]
#     coef = 1000
#     M = Adv + coef*Mass
#     b = np.random.random(nr)

#     Afun = lambda x: M @ x
#     Pfun = lambda x: invP @ x
#     Pfun = lambda x: 1/coef*invP @ x

#     x, resPCG = PBiCGSTAB(Afun, Pfun, b)
#     plt.semilogy(resPCG)
#     plt.savefig("IterSolvers")


# GMRES
invP = np.diag(1/np.diag(Mass))
b = np.random.random(np.size(Adv, axis=0))
for m in range(1, 5):
    x = GMRES(Mass, b, invP=invP, n_restarts=m, n_iter=20)
    print('error:%.3e'%(np.linalg.norm(b - Mass @ x)/np.linalg.norm(b)))

print('%===========%') 

for m in range(1, 5):
    x = GMRES(Mass, b, invP=invP, n_iter=m*20, n_restarts=1)
    print('error:%.3e'%(np.linalg.norm(b - Mass @ x)/np.linalg.norm(b)))

