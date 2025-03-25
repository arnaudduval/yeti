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

# Python module
import numpy as np
import scipy.sparse as sp
from ...solver import pseudoDense

class FETI:
    def __init__(self,interfaceData,dirichletData=None):
        self._set_lgrgemultiplierInfos(interfaceData)
        self._set_dirichletBCsInfos(dirichletData)
        self._build_assemblyMATRIX()

    def _set_lgrgemultiplierInfos(self,tabWeak):
        n     = len(tabWeak)
        self._nb_subdomain = n
        count = 0
        mult  = []
        for i in range(0,n-1):
            nbi = tabWeak[i].shape[1]
            for j in range(i+1,n):
                nbj = tabWeak[j].shape[1]

                tabi = np.repeat(tabWeak[i],nbj,axis=1)
                tabj = np.tile(  tabWeak[j],nbi)

                coords = np.all(np.isclose((tabi-tabj)[:-5,:],0.),axis=0)
                mtrslv = np.abs((tabi-tabj)[-5,:])==1.
                rotdisp= (tabi-tabj)[-4,:]==0.
                ndof   = (tabi-tabj)[-1,:]==0.
                test   = np.all(np.array([coords,mtrslv,rotdisp,ndof]),axis=0)
                if test.any():
                    ind = np.where(test)[0]
                    nn  = ind.size
                    count += nn
                    mult.append(np.array([i*np.ones(nn),tabi[-3,ind],
                                          j*np.ones(nn),tabj[-3,ind],tabi[-1,ind]],
                                         dtype=np.intp).transpose())
        self._interfaceInfos = np.concatenate(mult) #np.reshape(mult,(count,5))
        self._nbmult = count
        return None

    def _set_dirichletBCsInfos(self,tabBCs):
        if tabBCs == None:
            self._dirichletInfos = np.zeros(self._nb_subdomain,dtype=np.intp)
        else:
            self._dirichletInfos = np.array(tabBCs,dtype=np.intp)
        return None


    def _imap_lgrge2subdomain(self,ID):
        test  = np.where(np.isin(self._interfaceInfos[:,(0,2)],ID))
        isort = np.argsort(self._interfaceInfos[:,(1,3)][test])

        ndof  = np.intp(np.append(0,np.cumsum(self._interfaceInfos[:,-1])))
        num_interface = test[0]
        idof = []
        for lg in num_interface[isort]:
            idof.extend(np.arange(ndof[lg],ndof[lg+1]))
        return np.array(idof,dtype=np.intp)

    def _imap_dirichlet2subdomain(self,ID):
        offset = np.sum(self._interfaceInfos[:,-1])
        ndof   = np.intp(np.append(0,np.cumsum(self._dirichletInfos))) + offset
        return np.arange(ndof[ID],ndof[ID+1])

    def _build_assemblyMATRIX(self):
        listID = np.unique(self._interfaceInfos[:,(0,2)])
        self._dofl= np.sum(self._interfaceInfos[:,-1]) + np.sum(self._dirichletInfos)
        matrixA   = {}
        for ID in listID:
            rowL= self._imap_lgrge2subdomain(ID)
            rowD= self._imap_dirichlet2subdomain(ID)
            row = np.concatenate((rowL,rowD))
            col = np.arange(0,row.size)
            data= np.ones(row.size,dtype=np.float64)
            As  = sp.coo_matrix((data,(row,col)),shape=(self._dofl,row.size))
            matrixA.update({'%i'%ID:As.tocsc()})
        self._assemblyA = matrixA
        return None

    def set_coarseMATRIX(self,localG):
        Gtab = []
        ID   = 0
        for Gs in localG:
            Gtab.append(self._assemblyA['%i'%ID] * Gs)
            ID += 1
        self._matrixG = -sp.bmat([Gtab]).tocsc()
        self._matrixG.sort_indices()
        return None

    def set_coarsePrecondQ(self,matrixQ):
        if isinstance(matrixQ,sp.spmatrix):
            self._matrixQ = matrixQ.copy()
        else:
            self._matrixQ = matrixQ
        return None

    def solve_coarsePB(self):
        if not hasattr(self,'_matrixQ'):
            self._LUfeti = sp.linalg.splu((self._matrixG.T * self._matrixG).tocsc())
        else:
            #temp = []
            #for i in np.arange(0,self._matrixG.shape[1]):
            #    temp.append(self._matrixQ.dot( self._matrixG.getcol(i).toarray()[:,0] ))
            #temp = sp.bmat(temp).T
            temp = self._matrixQ * self._matrixG
            temp.sorted_indices()
            self._LUfeti_wQ = sp.linalg.splu((self._matrixG.T * temp ).tocsc())
        return None

    def _projectVect(self,vect):
        if not hasattr(self,'_LUfeti'):
            self.solve_coarsePB()
        y = np.zeros(self._dofl,dtype=np.float64)
        y[:] = vect[:] - self._matrixG * self._LUfeti.solve( self._matrixG.T * vect[:] )
        return y
    def _projectVect_wQ(self,vect):
        if not hasattr(self,'_LUfeti_wQ'):
            self.solve_coarsePB()
        y = np.zeros(self._dofl,dtype=np.float64)
        y[:] = vect[:] - self._matrixQ.dot(
            self._matrixG.dot(self._LUfeti_wQ.solve( self._matrixG.T * vect[:] )) )
        return y
    def _projectVect_wQ_t(self,vect):
        if not hasattr(self,'_LUfeti_wQ'):
            self.solve_coarsePB()
        temp = self._matrixG.T * ( self._matrixQ.T * vect )
        y = np.zeros(self._dofl,dtype=np.float64)
        y[:] = vect[:] - self._matrixG * self._LUfeti_wQ.solve( temp[:] )
        return y

    def set_projectorP(self):
        if not hasattr(self,'_matrixQ'):
            self.projectorP = sp.linalg.LinearOperator((self._dofl,self._dofl),
                                                       matvec=self._projectVect)
            self.projectorPt= self.projectorP
        else:
            self.projectorP = sp.linalg.LinearOperator((self._dofl,self._dofl),
                                                       matvec =self._projectVect_wQ)
            self.projectorPt= sp.linalg.LinearOperator((self._dofl,self._dofl),
                                                       matvec=self._projectVect_wQ_t)
        return None

    def set_RHSt(self,localt):
        t = np.zeros(self._dofl,dtype=np.float64)
        ID = 0
        for ts in localt:
            t[:] += self._assemblyA['%i'%ID].dot(ts)
            ID   += 1
        self._vectT = t
        return None

    def set_RHSe(self,locale):
        self._vectE = -np.concatenate(locale)
        return None

    def set_globalcouplingMATRIX(self,localC):
        Ctab = []
        ID   = 0
        for Cs in localC:
            Ctab.append(self._assemblyA['%i'%ID] * Cs)
            ID += 1
        self._matrixC = sp.bmat([Ctab]).tocsc()
        self._matrixC.sort_indices()
        return None

    def set_globalinvdiaKMATRIX(self,localinvDiagK):
        invkii = np.concatenate(localinvDiagK)
        self._invDiagK = sp.dia_matrix((invkii,[0]),(invkii.size,invkii.size)).tocsc()
        return None

    def _build_invglobalC(self,scaled=False):
        if not hasattr(self,'_matrixC'):
            print('Error: global coupling matrix does not exist')
            return None
        if scaled is True:
            if not hasattr(self,'_invDiagK'):
                print('Error: inverse of diagonalized stiffness matrix does not exist')
                return None
            else:
                self._invC = sp.linalg.splu(self._matrixC * self._invDiagK * self._matrixC.T)
        else:
            self._invC = sp.linalg.splu(self._matrixC * self._matrixC.T)
            return None

    def evaluatedispcorrection(self,residual,scaled=False):
        if not hasattr(self,'_invC'):
            self._build_invglobalC(scaled=scaled)
        deltaU = self._invC.solve(residual)
        return deltaU

