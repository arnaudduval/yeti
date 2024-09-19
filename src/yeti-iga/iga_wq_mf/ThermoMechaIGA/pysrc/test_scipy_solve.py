from numpy.random import default_rng
import numpy as np
from scipy.sparse.linalg import splu
import scipy.sparse as sparse
import time
from scipy.sparse.linalg import use_solver

def add_to_diag(res, vector):
    diag = sparse.eye(len(vector), format="coo") ##make variance
    diag.setdiag(vector) ##make variance
    sparse_covariance = res + diag  ##add variance
    return sparse_covariance

def maked(n,m):
    rng1 = default_rng()
    rows = rng1.choice(n, size=m)
    rng2 = default_rng()
    cols = rng2.choice(n, size=m)
    data = np.random.rand(m)
    return rows,cols,data

use_solver(useUmfpack=False)

N = 10000
n = 2000000

matrix = sparse.coo_matrix((N,N))

r,c,d = maked(N,n)

res = sparse.coo_matrix((d,(r,c)), dtype = np.float64)
res = res @ res.T

vector = np.random.rand(N) * 1.0
res = add_to_diag(res,vector)

print("Number of Non-Zero Entries: ",res.count_nonzero())

st = time.time()
#########LU###########
A = res.tocsc()
print("Transferred after: ", time.time() - st)
LU = splu(A)

print(LU)
print("concluded after: ", time.time() - st)