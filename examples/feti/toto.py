from solver import pseudoLUstep
import numpy as np
import scipy.sparse as sp

# A = sp.coo_matrix((np.array([1,2,1]),(np.array([0,1,3]),(np.array([0,1,3])))), shape=(4,4)).tocsc()
A = sp.coo_matrix((np.array([1,2e-9,2,1e-9,1,3]),(np.array([0,1,2,3,4,3]),(np.array([0,1,2,3,4,1])))), shape=(5,5)).tocsc()


print(A.data)
print(A.indptr)
print(A.indices)

LU_ = sp.linalg.splu(A,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0.,
                            options=dict(SymmetricMode=True))



print(LU_.U)
print("\n")
print(LU_.L)
print("\n")

tol = 1e-8
Ukk = np.abs(LU_.U.diagonal())
ind = np.where(Ukk<tol*Ukk.max())[0]


print(Ukk)
print("\n")
print(ind)
print("\n")

dataU   = LU_.U.data
indptrU = LU_.U.indptr
indicesU= LU_.U.indices
dataL   = LU_.L.data
indptrL = LU_.L.indptr
indicesL= LU_.L.indices

print(dataU)
print(indptrU)
print(indicesU)
print(dataL)
print(indptrL)
print(indicesL)


y = []
for i in np.array([3]):
    dataU   = LU_.U.data
    indptrU = LU_.U.indptr
    indicesU= LU_.U.indices
    dataL   = LU_.L.data
    indptrL = LU_.L.indptr
    indicesL= LU_.L.indices

    print(i)        
    print(id(dataU))
    print(id(LU_.U.data))

    #v1  = LU.U[:i,i]
    v2t = LU_.U[i,i+1:]
    print(dataU)
    dataU[indptrU[i+1]-1] = 1.                 #  U[ i,i ] = 1.
    print(dataU)
    dataU[indptrU[i]:indptrU[i+1]-1] = 0.       #  U[:i,i ] = 0. début de la colonne
    print(dataU)
    dataU[np.where(indicesU == i)[0][1:]] = 0.  #  U[ i,i+1:] = 0. fin de la ligne 
    print(np.where(indicesU == i))
    print(dataU)
            
    dataL[indptrL[i]+1:indptrL[i+1]] = 0.       #  L[i+1:,i ] = 0. fin de la colonne 
    dataL[np.where(indicesL == i)[0][:-1]] = 0. #  L[ i,:i-1] = 0. début de la ligne
            
    #y1 = sp.linalg.spsolve_triangular(LU.U[:i,:i].tocsr(),v1.toarray(),lower=False)
    #tmp= sp.linalg.splu(LU.U[:i,:i],permc_spec='NATURAL')
    #y1 = tmp.solve(v1.toarray())
    #y.append(y1)
    print(LU_.U)
    y2 = sp.linalg.spsolve_triangular(LU_.U[i+1:,i+1:].T, v2t.T.toarray(),lower=True )
    y.append(y2)

# LU = pseudoLUstep(A,tol=1.e-8)

# print(LU.R)
# print(LU.R.shape)
# print("\n")
# print(LU._LU.U)
# print(LU._LU.U.shape)
# print("\n")
# print(LU._LU.L)
# print(LU._LU.L.shape)