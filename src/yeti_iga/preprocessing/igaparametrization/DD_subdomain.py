# Copyright 2019-2020 Thibaut Hirschler

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
import scipy.linalg as sla

# IGA module
#from . import IGAparametrization
from ...solver import pseudoLU,pseudoDense,pseudoLUstep
from ...preprocessing.igaparametrization import IGAparametrization
from ...stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from ...coupling.cplgmatrix   import cplg_matrix,cplg_dirichlet,tabweakinfos
#from coupling.cplginfos    import tabweakinfos
from . import IGA_manipulation as manip



class IGAsubdomain:

    def __init__(self,igapara,ID):

        if not isinstance(igapara,IGAparametrization):
            raise ValueError('Invalid input > should be an IGAparametrization')
        self._modeleIGA = igapara  # IGA model
        self._ID        = ID       # num of subdomain
        self._maplgrge  = np.array([],dtype=np.intp) # map global lagrange to sub-domain
        self.set_grpByType()
        self.set_dofInfos()


    def set_grpByType(self):
        self._list_patch = np.where(
            np.isin(self._modeleIGA._ELT_TYPE,np.array(['U1','U3','U30'])))[0] + 1
        self._list_curve = np.where(np.isin(self._modeleIGA._ELT_TYPE,np.array(['U00'])))[0] + 1
        self._list_lgrge = np.where(np.isin(self._modeleIGA._ELT_TYPE,np.array([ 'U4'])))[0] + 1
        return None

    def set_dofInfos(self):
        ndof = self._modeleIGA._nb_dof_free
        idof = self._modeleIGA._ind_dof_free[:ndof]-1

        idof_internal = np.array([],dtype=np.intp)
        mcrd = self._modeleIGA._mcrd
        for p in self._list_patch-1:
            idof_internal = np.concatenate(
                (idof_internal,np.repeat(np.unique(self._modeleIGA._IEN[p])-1,mcrd)*mcrd
                 + np.tile(np.arange(0,mcrd),self._modeleIGA._indCPbyPatch[p].size)))
        self._idof_patch    = idof_internal.copy()
        self._idof_internal = np.intersect1d(idof_internal,idof)

        idof_lgrge = np.array([],dtype=np.intp)
        idof_lgrge_list = []
        for p in self._list_lgrge-1:
            temp = np.repeat(np.unique(self._modeleIGA._IEN[p])-1,mcrd)*mcrd \
                   + np.tile(np.arange(0,mcrd),self._modeleIGA._indCPbyPatch[p].size)
            idof_lgrge = np.concatenate((idof_lgrge,temp))
            idof_lgrge_list.append(np.intersect1d(temp,idof))

        self._idof_lgrge_tot = np.intersect1d(idof_lgrge,idof)
        self._idof_lgrge     = idof_lgrge_list

        return None

    def get_lgrge_ndof(self):
        ndof= np.array([i.size for i in self._idof_lgrge])
        tab = np.array([self._list_lgrge,ndof],dtype=np.intp)
        return tab

    def get_weakInfos(self,nb_pts=5):
        tab = tabweakinfos( *self._modeleIGA.get_inputs4cplginfos(nb_pts=nb_pts) )
        ndof= self.get_lgrge_ndof()
        return np.append(tab,ndof[1][np.newaxis,:],axis=0)

    def get_bcsInfos(self):
        idof_fixed = np.setdiff1d(self._idof_patch,self._idof_internal)
        return idof_fixed.size

    def set_mapLgrge(self,interfaceInfos):
        test  = np.where(np.isin(interfaceInfos[:,(0,2)],self._ID))
        isort = np.argsort(interfaceInfos[:,(1,3)][test])

        ndof  = np.intp(np.append(0,np.cumsum(interfaceInfos[:,-1])))
        num_interface = test[0]
        idof = []
        for lg in num_interface[isort]:
            idof.extend(np.arange(ndof[lg],ndof[lg+1]))
        self._maplgrge = np.array(idof)
        return None

    def set_stiffnessMATRIX(self,strongDirichlet=True):
        data,row,col,Fb = build_stiffmatrix( *self._modeleIGA.get_inputs4system_elemStorage() )
        Kside = sp.coo_matrix((data,(row,col)),
                              shape=(self._modeleIGA._nb_dof_tot,self._modeleIGA._nb_dof_tot),
                              dtype='float64').tocsc()
        Ktot = Kside + Kside.transpose()
        if strongDirichlet:
            self._K2solve = Ktot[self._idof_internal,:][:,self._idof_internal]
            self._f2solve =   Fb[self._idof_internal]
        else:
            self._K2solve = Ktot[self._idof_patch,:][:,self._idof_patch]
            self._f2solve =   Fb[self._idof_patch]
        return None

    def _build_weakDirichletMATRIX(self):
        Cdata,Crow,Ccol,Ub = cplg_dirichlet( *self._modeleIGA.get_inputs4dirichletmatrix() )
        Cside = sp.coo_matrix((Cdata,(Crow,Ccol)),
                               shape=(self._modeleIGA._nb_dof_bloq,self._modeleIGA._nb_dof_tot),
                               dtype='float64').tocsr()
        mask = np.isin(self._modeleIGA._ind_dof_bloq[:self._modeleIGA._nb_dof_bloq]-1,
                       self._idof_patch)
        Cbound = Cside[np.where(mask)[0],:][:,self._idof_patch]
        return Cbound


    def set_couplingMATRIX(self,strongDirichlet=True):
        Cdata,Crow,Ccol = cplg_matrix( *self._modeleIGA.get_inputs4cplgmatrix() )
        Cside = sp.coo_matrix((Cdata,(Crow,Ccol)),
                              shape=(self._modeleIGA._nb_dof_tot,self._modeleIGA._nb_dof_tot),
                              dtype='float64').tocsc()
        if strongDirichlet:
            self._C2solve = Cside[self._idof_internal,:][:,self._idof_lgrge_tot].transpose()
        else:
            Cbound = self._build_weakDirichletMATRIX()
            self._C2solve = sp.bmat(
                [[Cside[self._idof_patch,:][:,self._idof_lgrge_tot].transpose()],[Cbound]]).tocsr()
        self._C2solve.eliminate_zeros()
        if self._C2solve.format == 'csr':
            self._itracedisp = np.unique( self._C2solve.indices )
        else:
            self._itracedisp = np.unique( self._C2solve.tocsr().indices )
        return None

    def set_factorizationMATRIX(self,tol=1.e-08):
        #self._LU = sp.linalg.splu(self._K2solve)
        #self._LU = pseudoLU(self._K2solve,tol=1.e-5)
        #self._LU = pseudoDense(self._K2solve)
        if tol == 1.:
            self._LU = pseudoDense(self._K2solve)
        else:
            self._LU = pseudoLUstep(self._K2solve,tol=tol)
        return None


    def set_pseudoinverseMATRIX(self,rho=1.):
        Rtot   = self.get_unconstrainedRigidBodyModes()[self._idof_patch,:]

        Rgamma = sp.lil_matrix(Rtot.shape)
        Rgamma[self._itracedisp,:] = Rtot[self._itracedisp,:]
        Rgamma = Rgamma.tocsc()
        Bgamma = Rgamma.dot(Rgamma.T)

        self._LUstar = sp.linalg.splu(self._K2solve + rho*Bgamma)
        return None

    def set_pseudoinverseMATRIX_fixednodes(self,rho=None):
        vertex = []
        for patch in self._list_patch:
            indcp = manip.get_vertexCPindice(self._modeleIGA._Nkv,self._modeleIGA._Jpqr,
                                             self._modeleIGA._dim,num_patch=patch-1)
            vertex.extend(list(self._modeleIGA._indCPbyPatch[patch-1][indcp]-1))

        print('vertex',vertex)
        mcrd = self._modeleIGA._mcrd
        dofvertex = np.repeat(vertex,mcrd) + np.tile(np.arange(0,mcrd),len(vertex))
        Rtot = self.get_unconstrainedRigidBodyModes()

        RI = sp.lil_matrix((self._idof_patch.size,Rtot.shape[1]))
        RI[np.where(np.isin(self._idof_patch,dofvertex))[0],:] = Rtot[dofvertex,:]
        RI = RI.tocsc()
        TI = RI.T * RI
        invRI = sp.csc_matrix(sla.inv(TI.toarray()))
        BI = RI * invRI * RI.T
        if rho is None:
            rho = self._K2solve.diagonal().max()
        self._LUstar = sp.linalg.splu(self._K2solve + rho*BI)
        return None

    def compute_condensedRHSvect(self):
        if   hasattr(self,'_LU'):
            lu = self._LU
        elif hasattr(self,'_LUstar'):
            lu = self._LUstar
        else:
            return np.zeros(self._idof_lgrge_tot.size)
        ts = self._C2solve.dot(lu.solve(self._f2solve))
        return ts

    def compute_rigidbodyRHSvect(self):
        es = self._rigidbodyR.T.dot(self._f2solve)
        return es

    def evaluatedualshur(self,lmbda):
        if   hasattr(self,'_LU'):
            lu = self._LU
        elif hasattr(self,'_LUstar'):
            lu = self._LUstar
        else:
            return np.zeros(self._idof_lgrge_tot.size)
        ub = self._C2solve * lu.solve(self._C2solve.T * lmbda)
        return ub

    # For preconditionner
    def set_factorizationInternalMATRIX(self):
        ind = np.setdiff1d(np.arange(0,self._idof_internal.size),self._itracedisp)
        self._LUinternal = sp.linalg.splu(self._K2solve[ind,:][:,ind])
        self._Kbb = self._K2solve[self._itracedisp,:][:,self._itracedisp]
        self._Kib = self._K2solve[ind,:][:,self._itracedisp]
        self._Kbi = self._Kib.transpose()
        return None

    def evaluateprimalshur(self,ub):
        if not hasattr(self,'_LUinternal'):
            self.set_factorizationInternalMATRIX()
        lb = np.zeros(self._itracedisp.size)
        if ub.size == lb.size:
            y1 = self._Kib.dot(ub)
            y2 = self._LUinternal.solve(y1)
            y3 = self._Kbi.dot(y2)
            lb[:] = self._Kbb.dot(ub) - y3
        else:
            print('input has wrong size')
        return lb

    def evaluatedirichletprecond(self,rk):
        if not hasattr(self,'_factorR'):
            self.set_factorizationCPLGMATRIX()
        y1 = sla.solve_triangular(self._factorR,rk,trans=1)
        y2 = sla.solve_triangular(self._factorR,y1)
        y3 = self._C2solve.transpose().dot(y2)
        y4 = self.evaluateprimalshur(y3[self._itracedisp])
        y4tot = np.zeros(self._idof_internal.size)
        y4tot[self._itracedisp] = y4[:]
        y5 = self._C2solve.dot(y4tot)
        y6 = sla.solve_triangular(self._factorR,y5,trans=1)
        zk = sla.solve_triangular(self._factorR,y6)
        return zk


    def set_invdiaKMATRIX(self):
        diaK = 1./self._K2solve.diagonal()
        self._invDiagK = sp.dia_matrix((diaK,[0]),self._K2solve.shape).tocsc()
        return None

    def set_invcouplingMATRIX(self,scaled=False):
        if scaled is True:
            if not hasattr(self,'_invDiagK'):
                self.set_invdiaKMATRIX()
            self._invCCt   = sp.linalg.splu(
                self._C2solve * self._invDiagK * self._C2solve.transpose())
        else:

            #self._invCCt = sp.linalg.splu((self._C2solve * self._C2solve.transpose()).tocsc())
            self._invCCt = pseudoDense(self._C2solve * self._C2solve.transpose())
        return None

    def evaluateprimalshur_winvC(self,ub,scaled=False):
        if not hasattr(self,'_LUinternal'):
            self.set_factorizationInternalMATRIX()
        if not hasattr(self,'_invCCt'):
            self.set_invcouplingMATRIX(scaled=scaled)
        lb = np.zeros(ub.size)
        y0 = self._C2solve.T.dot( self._invCCt.solve(ub) )
        if scaled is True:
            y1 = self._invDiagK.dot(y0)[self._itracedisp]
        else:
            y1 = y0[self._itracedisp]
        y2 = self._Kib.dot(y1)
        y3 = self._LUinternal.solve(y2)
        y4 = self._Kbb.dot(y1) - self._Kbi.dot(y3)
        y4tot = np.zeros(self._C2solve.shape[1])
        y4tot[self._itracedisp] = y4[:]
        if scaled is True:
            lb[:] = self._invCCt.solve(self._C2solve.dot( self._invDiagK.dot(y4tot) ))
        else:
            lb[:] = self._invCCt.solve(self._C2solve.dot(y4tot))
        return lb

    def evaluatelumpedprimalshur(self,ub):
        if not hasattr(self,'_Kbb'):
            self._Kbb = self._K2solve[self._itracedisp,:][:,self._itracedisp]
        lb = np.zeros(ub.size)
        y0 = self._C2solve.T.dot( ub )
        y1 = self._Kbb.dot(y0[self._itracedisp])
        y1tot = np.zeros(self._C2solve.shape[1])
        y1tot[self._itracedisp] = y1[:]
        lb[:] = self._C2solve.dot(y1tot)
        return lb

    def set_factorizationCPLGMATRIX(self):
        P = self._C2solve[:,self._itracedisp]
        self._factorR = np.linalg.qr(P.transpose().toarray(),mode='r')
        return None

    def get_unconstrainedRigidBodyModes(self):
        nb_cp  = self._modeleIGA._nb_cp
        mcrd   = self._modeleIGA._mcrd
        COORDS = self._modeleIGA._COORDS.T
        Rtrl   = np.tile(np.identity(mcrd),(nb_cp,1))
        Rrotz  = np.cross(np.array([0,0,1]),COORDS)
        if mcrd==3:
            Rrotx  = np.cross(np.array([1,0,0]),COORDS)
            Rroty  = np.cross(np.array([0,1,0]),COORDS)

            #patch2interpol = np.where(self._modeleIGA._ELT_TYPE == 'U30')[0]
            #activePatch = np.zeros(self._modeleIGA._nb_patch,dtype=np.intp)
            #activePatch[patch2interpol] = 1
            #pts = getgrevabscphysicalcoords(
            #    *self._modeleIGA.get_inputs4grevphyscoords( activePatch=activePatch) )


            Rrot   = np.reshape(np.concatenate((Rrotx,Rroty,Rrotz),axis=1),(nb_cp*3,3))
        else:

            Rrot = np.vstack(Rrotz[:,:2].flatten())
        rigidbodyRtot = np.concatenate((Rtrl,Rrot),axis=1)
        return rigidbodyRtot

    def build_rigidBodyModes(self):
        if hasattr(self,'_LU'):
            self._rigidbodyR = self._LU.R
            #from scipy.linalg import null_space
            #self._rigidbodyR = sp.csc_matrix(null_space(self._K2solve.A))
        else:
            self._rigidbodyR = sp.csc_matrix(
                self.get_unconstrainedRigidBodyModes()[self._idof_patch])
        return None

    def set_admissibleconstMATRIX(self):
        if not hasattr(self,'_rigidbodyR'):
            self.build_rigidBodyModes()
        if not hasattr(self,'_C2solve'):
            self.set_couplingMATRIX()
        self._localG = self._C2solve * self._rigidbodyR
        return None

    def evaluate_displacement(self,lmbda,alpha):
        if   hasattr(self,'_LU'):
            lu = self._LU
        elif hasattr(self,'_LUstar'):
            lu = self._LUstar
        disp = lu.solve(self._f2solve - self._C2solve.T * lmbda) + self._rigidbodyR.dot(alpha)
        return disp

