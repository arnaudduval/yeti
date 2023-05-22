from lib.__init__ import *
from lib.lib_base import evalDersBasisPy, evalDersBasisFortran, createKnotVector, array2csr_matrix

degree = 1
nbel   = 10
knotvector = createKnotVector(degree, nbel)
knots = np.linspace(0, 1, 10)

pB0, pB1 = evalDersBasisPy(degree, knotvector, knots)
basis, indi, indj = evalDersBasisFortran(degree, knotvector, knots)
fB0 = array2csr_matrix(basis[:, 0], indi, indj, isfortran=True)
fB1 = array2csr_matrix(basis[:, 1], indi, indj, isfortran=True)
error = np.abs(pB1.todense() - fB1.todense())
print(np.max(error))

