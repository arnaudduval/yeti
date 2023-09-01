# Copyright 2019 Thibaut Hirschler

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the terms
# of the GNU Lesser General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Yeti. If not, see <https://www.gnu.org/licenses/>

import numpy as np
import scipy.sparse as sp

import scipy.linalg as la


class pseudoLU:
    '''
    Modification of the superLU factorization for singular Matrix
    - generate a pseudo-inverse callable through fct solve()
    - build of the null space of the matrix
    '''
    def __init__(self,A,tol=1.e-10):
        self._modifiedLU(A,tol=tol)

    def _modifiedLU(self,A,tol=1.e-10):
        #LU = sp.linalg.splu(A,permc_spec='NATURAL')

        LU = sp.linalg.splu(A,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0.,
                            options=dict(SymmetricMode=True))

        n  = LU.shape[0]
        Ltsave = sp.linalg.splu(LU.L.T.tocsc(),permc_spec='NATURAL')

        Ukk = np.abs(LU.U.diagonal())
        ind = np.where(Ukk<tol*Ukk.max())[0]

        y = []
        for i in ind:
            dataU   = LU.U.data
            indptrU = LU.U.indptr
            indicesU= LU.U.indices
            dataL   = LU.L.data
            indptrL = LU.L.indptr
            indicesL= LU.L.indices

            #v1  = LU.U[:i,i]
            v2t = LU.U[i,i+1:]
            dataU[indptrU[i+1]-1] = 1.                  #  U[ i,i ] = 1.
            dataU[indptrU[i]:indptrU[i+1]-1] = 0.       #  U[:i,i ] = 0.
            dataU[np.where(indicesU == i)[0][1:]] = 0.  #  U[ i,i+1:] = 0.

            dataL[indptrL[i]+1:indptrL[i+1]] = 0.       #  L[i+1:,i ] = 0.
            dataL[np.where(indicesL == i)[0][:-1]] = 0. #  L[ i,:i-1] = 0.

            #y1 = sp.linalg.spsolve_triangular(LU.U[:i,:i].tocsr(),v1.toarray(),lower=False)
            #tmp= sp.linalg.splu(LU.U[:i,:i],permc_spec='NATURAL')
            #y1 = tmp.solve(v1.toarray())
            #y.append(y1)
            y2 = sp.linalg.spsolve_triangular(LU.U[i+1:,i+1:].T, v2t.T.toarray(),lower=True )
            y.append(y2)

        LU.U.eliminate_zeros()
        LU.L.eliminate_zeros()

        Pc = sp.csc_matrix((np.ones(n,dtype=np.float64), (np.arange(n),LU.perm_c) ))
        Pr = sp.csc_matrix((np.ones(n,dtype=np.float64), (LU.perm_r,np.arange(n)) ))
        ng = len(y)
        dataR    = []
        indicesR = []
        indptrR  = [0]
        for i in np.arange(ng):
            yi = y[i].flatten()
            #dataR.append(-yi)
            dataR.append( 1.)
            dataR.append(-yi)
            #indicesR.append(np.arange(yi.size+1))
            indicesR.append(ind[i]+np.arange(yi.size+1))
            indptrR.extend([indptrR[i] + yi.size+1])
        if ng>0:
            dataR    = np.block(dataR)
            indicesR = np.block(indicesR)
            indptrR  = np.array(indptrR)
            Rleft    = sp.csc_matrix((dataR,indicesR,indptrR),shape=(n,ng))
            Rleft    = sp.csc_matrix(Ltsave.solve(Rleft.A))
            Rleft    = Rleft.T * Pr
            nleft    = sp.linalg.norm(Rleft,axis=1)
            for i in np.arange(ng):
                Rleft.data[Rleft.indptr[i]:Rleft.indptr[i+1]] /= nleft[i]
            #Rright   = Pc * sp.csc_matrix((dataR,indicesR,indptrR),shape=(n,ng))
            #nright   = sp.linalg.norm(Rright,axis=0)
            #for i in np.arange(ng):
            #    Rright.data[Rright.indptr[i]:Rright.indptr[i+1]] /= nright[i]
        else:
            Rleft   = sp.csc_matrix((0,n))
            #Rright   = sp.csc_matrix((n,0))


        self.perm_c = LU.perm_c.copy()
        self.perm_r = LU.perm_r.copy()
        self._Pc    = Pc.copy()
        self._Pr    = Pr.copy()
        self._U     = sp.linalg.splu(LU.U,permc_spec='NATURAL')
        self._L     = sp.linalg.splu(LU.L,permc_spec='NATURAL')
        self.U      = self._U.U   #self.U = LU.U.tocsr()
        self.L      = self._L.L   #self.L = LU.L.tocsr()
        #self.R      = Rright.copy()
        self.R      = Rleft.transpose()
        #self.R = sp.csc_matrix(la.null_space(A.toarray(),rcond=1.e-10))
        self.shape  = LU.shape
        self._nulleq= ind.copy()

        return None

    def _solve(self,b):
        # warning : spsolve_triangular is not efficient
        # routine solve is probably faster
        y1 = sp.linalg.spsolve_triangular(self.L,self._Pr*b)
        y2 = sp.linalg.spsolve_triangular(self.U,y1,lower=False)
        y2[self._nulleq] = 0.
        x  = self._Pc*y2
        return x

    def solve(self,b):
        y1 = self._L.solve(self._Pr*b)
        y1[self._nulleq] = 0.
        y2 = self._U.solve(y1)
        y2[self._nulleq] = 0.
        x  = self._Pc*y2
        return x





class pseudoDense:
    '''
    Use scipy dense library to build the pseudo inverse and the null space
    '''
    def __init__(self,A,tol=1.e-10):
        self._factorized(A,rcond=tol)

    def _factorized(self,A,rcond=1.e-10):
        B = la.pinvh(A.toarray())
        self._pinv = sp.csc_matrix(B)
        self._shape= B.shape
        self.R = sp.csc_matrix(la.null_space(A.toarray(),rcond=rcond))
        return None

    def solve(self,b):
        x = self._pinv.dot(b)
        return x



class pseudoLUstep:
    '''
    Successive superLU factorization for soliving singular linear system
    - generate a pseudo-inverse callable through fct solve()
    - build of the null space of the matrix
    '''
    def __init__(self,A,tol=1.e-08):
        self._stepLU(A,tol=tol)

    def _stepLU(self,A,tol=1.e-08):
        n    = A.shape[0]
        Ai   = A.copy()
        Ab   = A.copy()
        Ptot = sp.identity(n,format='csc')
        ngmax= 10 # -- maximal size of nullspace
        nullpiv = True

        i = 0
        while i<ngmax+1 and nullpiv:

            LU = sp.linalg.splu(Ai,permc_spec='MMD_AT_PLUS_A',diag_pivot_thresh=0.,
                                options=dict(SymmetricMode=True))

            if not np.all(LU.perm_r == LU.perm_c):
                print(LU.perm_r, LU.perm_c)
                print('warning: different left and right permutation')

            ind = np.where(np.abs(LU.U.diagonal())<tol*LU.U.diagonal().max())[0]
            if ind.size == 0:
                nullpiv = False
            else:
                ni = Ai.shape[0]
                Pr = sp.csc_matrix((np.ones(ni,dtype=np.float64), (LU.perm_r,np.arange(ni)) ))
                Pi = sp.bmat([[Pr,None],[None,sp.eye(n-ni)]],format='csc')

                isort = np.concatenate((np.setdiff1d(np.arange(ni),ind[0]),[ind[0]],
                                        np.arange(ni,n)))
                Ji = sp.csc_matrix((np.ones(n,dtype=np.float64),(np.arange(n),isort)))
                Ab = Ji * Pi * Ab * Pi.T * Ji.T
                Ptot = Ji * Pi * Ptot
                i += 1
                Ai = Ab[:-i,:-i]


        # build nullspace
        if i> 0:
            Rp = LU.solve(Ab[:-i,-i:].A)
            Id = np.identity(i)
            Rb = np.block([[-Rp],
                           [ Id]])
            Rb/= np.linalg.norm(Rb,axis=0)
            R  = sp.csc_matrix(Ptot.T * Rb)
        else:
            R  = sp.csc_matrix((n,0))

        self._LU = LU
        self._P  = Ptot
        self.R   = R
        return None

    def solve(self,b):
        n = self._LU.shape[0]
        nb= b.size
        if nb == n:
            x = self._LU.solve(b)
        else:
            bp = self._P * b
            ng = nb-n
            xp = self._LU.solve(bp[:-ng])
            x  = self._P.T * np.concatenate((xp,np.zeros(ng,dtype=np.float64)))
        return x
