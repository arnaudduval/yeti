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

import numpy as np
import scipy.sparse as sp


# resolution
def PCPGortho(A,b,x0=None,tol=1.e-5,maxiter=50,M=None,P=None,Pt=None,savetxt=False):
    '''
    Preconditioned Conjugate Projected Gradient
    -- with orthogonalization of directions --
    '''
    # initialization
    if x0 is None:
        xk = np.zeros(b.size)
    else:
        xk = x0.copy()
    if P is None:
        P = sp.identity(b.size)
    if Pt is None:
        Pt = P
    if M is None:
        M = sp.identity(b.size)
    if savetxt is True:
        saveresidual   = []
        saveresidualrk = []
        
    rk = b - A.dot(xk)
    wk = Pt.dot(rk)
    cvg= False; norm0 = np.linalg.norm(wk)
    k  = 0
    pi = [] # save grad
    Api_piApi = [] # save dot product
    while k<maxiter and not cvg:
        
        zk = M.dot(wk)
        yk = P.dot(zk)
        pk = yk.copy()
        for i in np.arange(0,k):
            bik = yk.dot(Api_piApi[i])
            pk -= bik*pi[i]
        
        pi.append(pk.copy())
        Apk   = A.dot(pk)
        pkApk = pk.dot(Apk)
        
        Api_piApi.append(Apk/pkApk)
        
        alphak= yk.dot(wk)/pkApk
        
        xk += alphak*pk
        rk -= alphak*Apk
        wk  = Pt.dot(rk)

        k += 1
        normk = np.linalg.norm(wk)
        cvg   = normk<tol*norm0
        if savetxt is True:
            saveresidual.append(normk/norm0)
            saveresidualrk.append(np.linalg.norm(rk)/norm0)
            
    info = k
    if savetxt is True:
        np.savetxt('pcpg_cvrg.txt',np.array([saveresidual,saveresidualrk]).T,
                   delimiter=',')
    return xk,info
