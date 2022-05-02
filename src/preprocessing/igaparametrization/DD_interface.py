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

'''
Quelques fonctions pour ajouter des interfaces de couplage
'''

import numpy as np

from . import IGA_parametrization as igapara
from . import IGA_manipulation    as manip


def add1Dinterface(igapara,numpatch,COORDS,vectW,Uknots,degree,masterslave=1,disprot=0):
    '''
    Add one dimensional interface to a patch of an IGA model
    Entry:
    * igapara ...... an IGA model containing the patch
    * numpatch ..... number of patch where the interface curve is embedded
    * COORDS ....... Coordinates of the control points of the curve
    * vectW ........ Weights of the control points of the curve
    * Uknots ....... Knot vector of the curve
    * degree ....... Degree of the curve
    Option:
    * masterslave .. master=1, slave=0 (Default=1)
    * disprot ......   disp=0,   rot=1 (Default=0)
    Output
    * None --> properties of igapara are directly modified
    '''
    print(' Add interface to patch %2d -- info : %s, %s' % (numpatch,
                                                            'master' if masterslave else 'slave',
                                                            'rot' if disprot else 'disp'))
    
    n = igapara._nb_patch
    propsedge = np.array([n+1,numpatch,n+2,masterslave],dtype=np.float64)
    k0       = np.arange(degree,Uknots[0].size-(degree+1))
    test0    = Uknots[0][k0] != Uknots[0][k0+1]
    list0    = degree+1 + np.where(test0)[0]
    Nijkedge = np.array([list0,np.zeros_like(list0),np.zeros_like(list0)])
    nb_elem  = Uknots[0].size-(degree+1)*2 + 1
    ienxi    = np.tile(np.arange(0,degree+1), (nb_elem,1)) \
               + np.tile(np.vstack(np.arange(0,nb_elem)), (1,degree+1))
    IENedge  = np.arange(1,COORDS.shape[1]+1)[ienxi[:,::-1]]
    weightedge = []
    for num_elem in np.arange(0,nb_elem):
        weightedge.append(vectW[IENedge[num_elem]-1])

    ulgrge = [np.array([0.,1.])]
    plgrge = 0
    propslgrge= np.array([n+2,disprot],dtype=np.float64)
    k0       = np.arange(plgrge,ulgrge[0].size-(plgrge+1))
    test0    = ulgrge[0][k0] != ulgrge[0][k0+1]
    list0    = plgrge+1 + np.where(test0)[0]
    Nijklgrge= np.array([list0,np.zeros_like(list0),np.zeros_like(list0)])
    IENlgrge = np.array([[COORDS.shape[1]+1]])
    weightlgrge = np.array([1.])
    
    
    elt_type = np.array(['U00','U4'])
    nbpint   = np.array([6,1])
    tensor   = np.array(['NONE','NONE'])
    COORDStot= np.block([COORDS,np.zeros((3,1))])
    nb_cp    = COORDStot.shape[1]
    ien      = [IENedge, IENlgrge]
    matprops = np.zeros((3,2))
    props    = [propsedge,propslgrge]
    jprops   = np.array([4,2])
    dim      = np.array([1,1])
    nkv      = np.array([[Uknots[0].size,0,0],[ulgrge[0].size,0,0]]).T
    ukv      = [Uknots,ulgrge]
    jpqr     = np.array([[degree,0,0],[0,0,0]]).T
    nijk     = np.block([Nijkedge,Nijklgrge])
    weights  = weightedge
    weights.append(weightlgrge)
    elmtsbyPatch = np.array([ien[0].shape[0],ien[1].shape[0]])
    nb_patch = 2
    nnode    = np.array([degree+1,1])
    nb_elem  = np.sum(elmtsbyPatch)
    
    mechSet = []
    mechSet.append([elt_type,nbpint,tensor,3])
    mechSet.append([COORDStot,nb_cp])
    mechSet.append( ien )
    mechSet.append([matprops])
    mechSet.append([props,jprops])
    
    geoSet = [dim,nkv,ukv,jpqr,nijk,weights,elmtsbyPatch,nb_patch,nnode,nb_elem]
    
    igapara.add_multiplePatch(mechSet,geoSet)

    for i in [1,2,3]:
        manip.add_displacementBC(igapara,igapara._indCPbyPatch[-2], i, 0.) # crve1
    if disprot == 1:
        for i in [2,3]:
            manip.add_displacementBC(igapara,igapara._indCPbyPatch[-1], i, 0.) # lgrge rot
    
    return None


def autosetCPLG4embeddedpatches(surfaceIGA,embeddedIGA):
    '''
    Automatic definition of the interface between an embedded KL shell and a standard KL shell.
    - One face of the mapping involved for the embedded KL shell should match with the mid-surface 
      of the standard KL shell.
    Entry:
    * surfaceIGA ... an IGA model which contains a standard KL shell (Elt_type 'U3' )
    * embeddedIGA .. an IGA model which contains a embdedd KL shell  (Elt_type 'U30')
    Output:
    * None --> new interfaces (curve + lgrge) are directly added to both IGA model
    '''
    print('\nAdd interface between embedded KL and standard KL')
    # find face where the mapping and the surface match
    numSurf = np.where( surfaceIGA._ELT_TYPE == 'U3' )[0][0]
    cpSurf  = surfaceIGA._COORDS[:,surfaceIGA._indCPbyPatch[numSurf]-1]

    matchingface = 0
    numMap  = np.where( embeddedIGA._ELT_TYPE == 'U0' )[0][0]
    for face in np.arange(1,7):
        iface = manip.get_directionCP(embeddedIGA, face, numMap)-1
        cpMap = embeddedIGA._COORDS[:,iface]
        if cpSurf.shape == cpMap.shape:
            if np.all(np.isclose(cpSurf,cpMap)):
                matchingface = face
                break
    if matchingface == 0:
        print(' Warning: no matching face found')
        return None

    # find edge of the embedded surface which lies on the matching face
    numEmb  = np.where( embeddedIGA._ELT_TYPE == 'U30' )[0][0]
    if matchingface == 1:
        ei = 0; val = 0.
    if matchingface == 2:
        ei = 0; val = 1.
    if matchingface == 3:
        ei = 1; val = 0.
    if matchingface == 4:
        ei = 1; val = 1.
    if matchingface == 5:
        ei = 2; val = 0.
    if matchingface == 6:
        ei = 2; val = 1.  
    for edge in np.arange(4):
        temp = manip.get_boundCPindice_wEdges(embeddedIGA._Nkv,embeddedIGA._Jpqr,embeddedIGA._dim,
                                              101+edge, num_patch=numEmb)
        iedge = embeddedIGA._indCPbyPatch[numEmb][temp]-1
        if np.all(embeddedIGA._COORDS[ei,iedge] == val):
            findedge = 101+edge
            COORDSedge = np.zeros((3,iedge.size))
            COORDSedge[:2,:] = np.delete(embeddedIGA._COORDS[:,iedge],ei,axis=0)
            vweightedge = embeddedIGA._vectWeight[iedge]
            break
            
    # exctract infos 
    print(' Standard KL')
    if findedge == 101 or findedge == 102:
        side = 0
    if findedge == 103 or findedge == 104:
        side = 1
    uedge = [embeddedIGA._Ukv[numEmb][side]]
    pedge = embeddedIGA._Jpqr[side,numEmb]
    
    add1Dinterface(surfaceIGA,numSurf+1,COORDSedge,vweightedge,uedge,pedge,masterslave=1)
    add1Dinterface(surfaceIGA,numSurf+1,COORDSedge,vweightedge,uedge,pedge,masterslave=1,disprot=1)
    #for i in [1,2,3]:
    #    manip.add_displacementBC(surfaceIGA,surfaceIGA._indCPbyPatch[-2], i, 0.) # crve1
    #    manip.add_displacementBC(surfaceIGA,surfaceIGA._indCPbyPatch[-4], i, 0.) # crve2
    #for i in [2,3]:
    #    manip.add_displacementBC(surfaceIGA,surfaceIGA._indCPbyPatch[-1], i, 0.) # lgrge rot
        
    print(' Embedded KL')
    uedge2 = [np.array([0.,0.,1.,1.])]
    if findedge == 101:
        COORDSedge2 = np.array([[0.,0.,0.],[1.,0.,0.]]).T
    if findedge == 102:
        COORDSedge2 = np.array([[0.,1.,0.],[1.,1.,0.]]).T
    if findedge == 103:
        COORDSedge2 = np.array([[0.,0.,0.],[0.,1.,0.]]).T
    if findedge == 104:
        COORDSedge2 = np.array([[1.,0.,0.],[1.,1.,0.]]).T
    vweightedge2 = np.array([1.,1.])
    
    add1Dinterface(embeddedIGA,numEmb+1,COORDSedge2,vweightedge2,uedge2,1,masterslave=0)
    add1Dinterface(embeddedIGA,numEmb+1,COORDSedge2,vweightedge2,uedge2,1,masterslave=0,disprot=1)
    #for i in [1,2,3]:
    #    manip.add_displacementBC(embeddedIGA,embeddedIGA._indCPbyPatch[-2], i, 0.) # crve1
    #    manip.add_displacementBC(embeddedIGA,embeddedIGA._indCPbyPatch[-4], i, 0.) # crve2
    #for i in [2,3]:
    #    manip.add_displacementBC(embeddedIGA,embeddedIGA._indCPbyPatch[-1], i, 0.) # lgrge rot
        
    return None
    




def automaticDomainDecomposition(igaPara):
    '''
    Automatic decomposition of a IGA model. Split the patches that are weakly coupled.
    Entry:
    * igaPara ...... an IGA model which possibly contains multiple weakly coupled patches 
    Output:
    * multipleIGA .. a list of independent IGA models
    '''
    
    print('\nAutomatic domain decomposition')
    
    # interfaces
    tabWeak = manip.get_tabCLG(igaPara)
    
    # test couplage fort
    tabStrong = manip.get_patchConnectionInfos(
        igaPara._Nkv,igaPara._Jpqr,igaPara._dim,igaPara._indCPbyPatch,igaPara._nb_patch)
    
    linkedPatches = []
    if tabStrong.size == 0:
        alonePatches = np.arange(1,igaPara._nb_patch+1)
    else:
        for cplg in tabStrong[:,(0,2)]+1:
            index   = []
            for i in np.arange(0,len(linkedPatches)):
                grp     = linkedPatches[i]
                isingrp = np.intersect1d(cplg,grp)
                if isingrp.size > 0:
                    index.append(i)
                
            if np.size(index)==0:
                linkedPatches.append(cplg)
            else:
                merge = list(cplg)
                for j in index:
                    merge.extend(linkedPatches[j])
                for j in index:
                    del linkedPatches[j]
                linkedPatches.append(np.unique(merge))
        alonePatches = np.setxor1d(np.arange(1,igaPara._nb_patch+1),np.concatenate(linkedPatches))
        
    for patch in alonePatches:
        linkedPatches.append(np.array([patch],dtype=np.intp))
        
    # sous-domaines
    # - nettoyage
    domain = []
    for grp in linkedPatches:
        mask = np.isin(igaPara._ELT_TYPE[grp-1],np.array(['U1','U3','U30']))
        if mask.any():
            domain.append(grp[mask])
            
    # - completer si couplage faible ou si composition de NURBS
    for d in range(0,len(domain)):
        toadd = []
        for patch in domain[d]:
            ind = np.where(tabWeak[:,1:3] == patch)
            toadd.extend(tabWeak[ind[0],0])
            toadd.extend(tabWeak[ind[0],ind[1]+3])

            if igaPara._ELT_TYPE[patch-1] == 'U30':
                toadd.extend( [int(igaPara._PROPS[patch-1][1])] )
            
        domain[d] = np.concatenate((domain[d],np.unique(toadd).astype(np.intp)))
        domain[d].sort()
        
    deletePatches = np.setxor1d(np.arange(0,igaPara._nb_patch)+1,np.concatenate(domain))
    if deletePatches.size>0:
        print('\n**Warning: the following patches are ommited',deletePatches)
        
    # extract modele IGA
    multipleIGA = []
    for subdomain in domain:
        mechSet = igaPara.get_mechanicalSettings_somePatch(subdomain-1,withBC=True,
                                                           updatePROPS=True)
        geoSet  = igaPara.get_geometricSettings_somePatch( subdomain-1)
        multipleIGA.append(igapara.IGAparametrization(
            mechanicalSettings=mechSet,geometricSettings=geoSet))
        
    return multipleIGA
