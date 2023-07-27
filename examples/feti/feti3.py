

# Python module
import numpy as np
import scipy.sparse as sp
import sys
import time

#IGA module
from solver import PCPGortho
from preprocessing.igaparametrization import IGAparametrization
from preprocessing.igaparametrization import IGAsubdomain,FETI
from preprocessing.igaparametrization import IGAmanip    as manip
from preprocessing.igaparametrization import DDinterface as ddmanip
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# Selection of .INP and .NB file
# ------------------------------
DIRECTORY='inputFiles/'
CASE=[]
CASE.append('plateStiff') # --------------  0
CASE.append('beamEmbded') # --------------  1
CASE.append('cylindreTroue_embded') # ----  2
CASE.append('scordelis-Lo2cplg') # -------  3
CASE.append('Tbeam2cplg') # --------------  4
CASE.append('plate2cplg') # --------------  5
CASE.append('stifftube') # ---------------  6
CASE.append('scordelis-Lo_3patches') # ---  7
CASE.append('curvedwallstrangeStiff') # --  8
CASE.append('wing2') # -------------------  9
CASE.append('twoplatesDD') # ------------- 10
CASE.append('threeplatesDD') # ----------- 11
CASE.append('scordelis-Lo_2mixpatches') #- 12
CASE.append('stiffenedmultiPanel_DD') # -- 13
CASE.append('stiffroof') # --------------- 14
CASE.append('wingcplg') # ---------------- 15
CASE.append('wingcplg2') # --------------- 16
CASE.append('supportedroof') # ----------- 17
CASE.append('beamKL2cplg') # ------------- 18
CASE.append('fourplatesDD') # ------------ 19
CASE.append('fourplatesDDv2') # ---------- 20
#CASE.append('fourembeddedplatesDD') # ---------- 20
CASE.append('fiveplatesDD') # ------------ 21
CASE.append('threeembdedDD') # ----------- 22
CASE.append('fourembdedDD') # ------------ 23
CASE.append('nineembdedDD') # ------------ 24
CASE.append('stiffenedcurvedPanel_DD') # - 25

EXEMPLE_NB = 11

FILENAME = DIRECTORY+CASE[EXEMPLE_NB]


# Decomposition
# -------------

completeIGA = IGAparametrization(filename=FILENAME)

nb_deg = np.zeros((3,completeIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,completeIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}
if EXEMPLE_NB == 3:
    p = 2
    r = 4
    # cylinder
    nb_deg[:2,0] = [p-1,p]
    nb_deg[:2,1] = [p-1,p]
    nb_ref[:2,0] = [r,r]
    nb_ref[:2,1] = [r,r]
    # linked curves
    nb_ref[0,np.array([2,3,4,5])] = r+2
    # lgrge
    nb_deg[0,6] = p
    nb_deg[0,7] = p-1
    nb_ref[0,np.array([6,7])] = r
if EXEMPLE_NB == 4:
    p = 3
    nb_deg[:2,:2] = p
    r = 2
    nb_ref[:,0]  = np.array([3+r,2+r,0])
    additional_knots["patches"] = np.array([1])
    additional_knots["1"] = np.linspace(0,1,10*2**r)[1:-1]
    additional_knots["2"] = np.linspace(0,1, 4*2**r)[1:-1]
    nb_ref[0,2:6]= 5+r
    nb_deg[0,6:] = np.array([p,p-1])
    nb_ref[0,6:] = 3+r
if EXEMPLE_NB == 5:
    p = 3
    r = 2
    nb_deg[:,0] = np.array([p-1,p,0])
    nb_deg[:,1] = np.array([p,p-1,0])

    nb_ref[:,0]  = np.array([r,r,0])
    nb_ref[:,1]  = np.array([r+1,r+1,0])

    nb_ref[0,2:6]= r+2

    nb_deg[0,6:] = np.array([p,p-1])
    nb_ref[0,6:] = r+1
if EXEMPLE_NB == 7:
    p = 1
    r = 4
    # cylinder
    nb_deg[ 0,1:4] = p
    nb_deg[ 1,1:4] = p
    nb_ref[ 0,1:4] = r
    nb_ref[ 1,1:4] = r
    additional_knots = {"patches":np.array([3]),"1":np.array([0.45]),"2":np.array([0.55]),
                        "3":np.array([])}
    nb_deg[:2,3] = p+1
    nb_ref[:2,3] = r-1
    nb_ref[:2,2] = r+1
    # linked curves
    nb_ref[0,np.array([5,6,7,8,11,12,13,14,17,18,19,20])-1] = r+1
    # lgrge
    nb_deg[0,np.array([ 9,15,21])-1] = np.max(p+1,0)
    nb_deg[0,np.array([10,16,22])-1] = np.max(p+0,0)
    nb_ref[0,np.array([ 9,15,21,10,16,22])-1] = r
if EXEMPLE_NB == 8:
    p = 1
    r = 5
    # wall
    nb_deg[:2,1] = p
    nb_ref[:2,1] = r
    # stiff
    nb_deg[ 0,2] = p
    nb_deg[ 1,2] = p+1
    nb_ref[ 0,2] = np.maximum(0,r-2)
    nb_ref[ 1,2] = r+1
    # curves
    nb_deg[0,np.array([3,5,9,11])] = p+1
    nb_ref[0,np.array([3,4,5,6,9,10,11,12])] = r+1
    # lgrge
    nb_deg[0,np.array([ 7, 8])] = p+1
    nb_ref[0,np.array([ 7, 8])] = r
    nb_deg[0,np.array([13,14])] = p
    nb_ref[0,np.array([13,14])] = r
if EXEMPLE_NB == 10:
    p = 2
    r = 2
    # domains
    nb_deg[:2,:2] = p
    nb_ref[:2, 0] = r+1
    nb_ref[:2, 1] = r
    # curves
    nb_ref[0,(2,3)] = r+2
    # lgrge
    nb_deg[0,4] = p
    nb_ref[0,4] = r
if EXEMPLE_NB == 11:
    p = 2
    r = 3
    # domains
    nb_deg[:2,:3] = p
    nb_ref[:2, 0] = r
    nb_ref[:2, 1] = r+1
    nb_ref[:2, 2] = r
    # curves
    nb_ref[0,(3,4,5,6)] = r+2
    # lgrge
    nb_deg[0,(7,8)] = p
    nb_ref[0,(7,8)] = r
if EXEMPLE_NB == 12:
    p = 3
    r = 4
    # cylinder
    nb_deg[:2,0] = [p-1,p]
    nb_deg[:2,2] = p-1
    nb_ref[:2,0] = [r,r]
    nb_ref[:2,2] = [r,r]
    # linked curves
    nb_ref[0,np.array([3,4,5,6])] = r+2
    # lgrge
    nb_deg[0,7] = p
    nb_deg[0,8] = p-1
    nb_ref[0,np.array([7,8])] = r
if EXEMPLE_NB == 13:
    p = 2
    r = 2
    # domains
    nb_deg[:2,:2] = p
    nb_ref[:2, 0] = r
    nb_ref[:2, 1] = [r+1,r-2]
    nb_deg[:2, 8] = p
    nb_ref[:2, 8] = [r+1,r-2]
    nb_deg[:2,15] = p
    nb_ref[:2,15] = [r+1,r-2]
    # curves
    nb_ref[0,( 2, 3, 4, 5)] = r+2
    nb_ref[0,( 9,10,11,12)] = r+2
    nb_ref[0,(16,17,18,19)] = r+2
    nb_ref[0,(22,23,24,25)] = r
    nb_ref[0,(28,29,30,31)] = r
    # lgrge
    nb_deg[0,6]  = p
    nb_deg[0,7]  = p-1
    nb_ref[0,(6,7)]   = r
    nb_deg[0,13] = p
    nb_deg[0,14] = p-1
    nb_ref[0,(13,14)] = r
    nb_deg[0,20] = p
    nb_deg[0,21] = p-1
    nb_ref[0,(20,21)] = r
    nb_deg[0,26] = p
    nb_deg[0,27] = p-1
    nb_ref[0,(26,27)] = r-2
    nb_deg[0,32] = p
    nb_deg[0,33] = p-1
    nb_ref[0,(32,33)] = r-2
if EXEMPLE_NB == 14:
    p = 2
    r = 6
    # plate
    nb_deg[:2,:2] = p
    nb_ref[:2,:2] = r
    # stiffeners and linked curves
    nb_deg[0,np.array([3,4,6,10,11,13])-1] = p
    nb_ref[0,np.array([3,4,6,10,11,13])-1] = r
    nb_deg[1,np.array([3,10])-1] = p
    nb_ref[1,np.array([3,10])-1] = r-2
    nb_ref[0,np.array([5,7,12,14])-1] = r+1
    # lgrge
    nb_deg[0,np.array([8,15])-1] = p
    nb_deg[0,np.array([9,16])-1] = np.maximum(p-1,0)
    nb_ref[0,np.array([8,9,15,16])-1] = r
if EXEMPLE_NB == 15:
    p = 2
    r = 5
    # domains
    nb_deg[:2, 0] = [np.maximum(p-2,0),p]
    nb_ref[:2, 0] = [np.maximum(r-3,0),r]
    nb_deg[:2, 1] = [np.maximum(p-2,0),p]
    nb_ref[:2, 1] = [np.maximum(r-3,0),r]
    # curves
    nb_ref[0,(2,3,4,5,6,7,8,9)] = r+1
    # lgrge
    nb_deg[0,(10,12)] = p
    nb_deg[0,(11,13)] = p-1
    nb_ref[0,(10,11,12,13)] = r
if EXEMPLE_NB == 16:
    p = 2
    r = 5
    # domains
    nb_deg[:2, 0] = [np.maximum(p-2,0),p]
    nb_ref[:2, 0] = [np.maximum(r-2,0),r]
    nb_deg[:2, 1] = [np.maximum(p-2,0),p]
    nb_ref[:2, 1] = [np.maximum(r-2,0),r]
    nb_deg[:2, 2] = [np.maximum(p-2,0),p]
    nb_ref[:2, 2] = [np.maximum(r-3,0),r]
    nb_deg[:2, 3] = [np.maximum(p-2,0),p]
    nb_ref[:2, 3] = [np.maximum(r-3,0),r]
    # curves
    nb_ref[0,( 4, 5, 6, 7, 8, 9,10,11)] = r+1
    nb_ref[0,(16,17,18,19,20,21,22,23)] = r+1
    nb_ref[0,(28,29,30,31,32,33,34,35)] = r+1
    # lgrge
    nb_deg[0,(12,14,24,26)] = p
    nb_deg[0,(13,15,25,27)] = p-1
    nb_ref[0,(12,13,14,15,24,25,26,27)] = r

    nb_deg[0,(36,38)] = p
    nb_deg[0,(37,39)] = p-1
    nb_ref[0,(36,37,38,39)] = r-1
if EXEMPLE_NB == 17:
    p = 2
    r = 3
    # roof and beam
    nb_deg[:,0] = np.array([p,p,0])
    nb_deg[:,1] = np.array([p,p,0])
    nb_ref[:,0]  = np.array([r+2,r+2,0])
    nb_ref[:,1]  = np.array([r  ,r+1,0])
    # curves
    nb_ref[0,2:6]= r+1
    # lgrge
    nb_deg[0,6:] = np.array([p,p-1])
    nb_ref[0,6:] = r
if EXEMPLE_NB == 18:
    p = 2
    r = 2
    # shell
    nb_deg[:,0] = np.array([p,p,0])
    nb_deg[:,1] = np.array([p,p,0])
    nb_ref[:,0]  = np.array([r+1,r+1,0])
    nb_ref[:,1]  = np.array([r,r,0])
    # curve
    nb_ref[0,2:6]= r+2
    # lgrge
    if True:
        nb_deg[0,6:] = np.array([p,p])+1
        nb_ref[0,6:] = r+1
    else:
        nb_deg[0,6:] = np.array([p,p-1])
        nb_ref[0,6:] = r
if EXEMPLE_NB == 19:
    p = 2
    r = 3
    # domains
    nb_deg[:2,:4] = p
    nb_ref[:2, 0] = r
    nb_ref[:2, 1] = r+1
    nb_ref[:2, 2] = r
    nb_ref[:2, 3] = r+1
    # curves
    nb_ref[0,(4,5,6,7,8,9)] = r+2
    # lgrge
    nb_deg[0,(10,11,12)] = p
    nb_ref[0,(10,11,12)] = r
if EXEMPLE_NB == 20:
    p = 2
    r = 3
    # domains
    nb_deg[:2,:4] = p
    nb_ref[:2, 0] = r+1
    nb_ref[:2, 1] =[r+1,r]
    nb_ref[:2, 2] = r
    nb_ref[:2, 3] = r
    # curves
    nb_ref[0,(4,5,6,7,8,9,10,11,12,13)] = r+2
    # lgrge
    nb_deg[0,(14,15,16,17,18)] = p
    nb_ref[0,(14,15,16,17,18)] = r
if EXEMPLE_NB == 21:
    p = 2
    r = 3
    # domains
    nb_deg[:2,:5] = p
    nb_ref[:2, 0] = r
    nb_ref[:2, 1] = r+1
    nb_ref[:2, 2] = r
    nb_ref[:2, 3] = r+1
    nb_ref[:2, 4] = r
    # curves
    nb_ref[0,(5,6,7,8,9,10,11,12)] = r+2
    # lgrge
    nb_deg[0,(13,14,15,16)] = p
    nb_ref[0,(13,14,15,16)] = r
if EXEMPLE_NB == 22:
    p = 2
    r = 3
    # domains
    nb_deg[1,:3] = 1
    completeIGA.refine(nb_ref,nb_deg,additional_knots)
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[0][3]-1] = 0.1
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[1][2]-1] = 0.1
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[1][3]-1] = 0.9
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[2][2]-1] = 0.9


    nb_deg[ 0,:3] = p
    nb_deg[ 1,:3] = p-1
    nb_ref[:2, 0] =[r-1,r+1]
    nb_ref[:2, 2] =[r-1,r+1]
    nb_ref[:2, 1] =[r,0]
    additional_knots["patches"] = np.array([1])
    additional_knots["2"] = np.linspace(0,1,2**(r+1)-1)[1:-1]
    nb_ref[:2, 2] =[r-1,r+1]
    # curves
    nb_ref[0,(3, 4, 5, 6)] = r+2
    nb_ref[0,(9,10,11,12)] = r+2
    # lgrge
    nb_deg[0,( 7, 8)] = p
    nb_ref[0,( 7, 8)] = r+1
    nb_deg[0,(13,14)] = p-1
    nb_ref[0,(13,14)] = r+1
if EXEMPLE_NB == 23:
    p = 2
    r = 2
    # domains
    nb_deg[:2,:4] = p
    nb_ref[:2, 0] = r
    #nb_ref[:2, 1] = r
    #nb_ref[:2, 2] = r
    additional_knots["patches"] = np.array([1,2])
    additional_knots["1"] = np.linspace(0,1,2**r+2)[1:-1]
    additional_knots["2"] = np.linspace(0,1,2**r+2)[1:-1]
    nb_ref[:2, 3] = r
    # curves
    nb_ref[0,( 4, 5, 6, 7, 8, 9,10,11)] = r+2
    nb_ref[0,(16,17,18,19,20,21,22,23)] = r+2
    # lgrge
    nb_deg[0,(12,13,14,15)] = p
    nb_ref[0,(12,13,14,15)] = r
    nb_deg[0,(24,25,26,27)] = p-1
    nb_ref[0,(24,25,26,27)] = r
if EXEMPLE_NB == 24:
    p = 2
    r = 2
    # domains
    nb_deg[:2,(0,2,4,6,8)] = p+1
    #nb_ref[:2,(0,2,4,6,8)] = r
    nb_deg[:2,(1,3,5,7)]   = p
    nb_ref[:2,(1,3,5,7)]   = r+1
    additional_knots["patches"] = np.array([0,2,4,6,8])
    additional_knots["1"] = np.linspace(0,1,2**r+2)[1:-1]
    additional_knots["2"] = np.linspace(0,1,2**r+2)[1:-1]
    # curves
    nb_ref[0,np.arange( 9,33)] = r+2
    nb_ref[0,np.arange(45,69)] = r+2
    # lgrge
    nb_deg[0,np.arange(33,45)] = p
    nb_ref[0,np.arange(33,45)] = r+1
    nb_deg[0,np.arange(69,81)] = p-1
    nb_ref[0,np.arange(69,81)] = r+1
if EXEMPLE_NB == 25:
    p = 2
    r = 3

    # shape modification
    nb_deg[0,(1,2,4, 8,9,11)] = 1
    completeIGA.refine(nb_ref,nb_deg,additional_knots)
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[1][list([1,4])]-1] = 0.5
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[2][1]-1] = 0.5
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[23]  -1] = 0.375
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[8][list([1,4])]-1] = 0.5
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[9][1]-1] = 0.5
    completeIGA._COORDS[0,completeIGA._indCPbyPatch[29]  -1] = 0.625

    # domains
    nb_deg[:,:] = 0
    nb_deg[:2, 0] = p
    nb_ref[:2, 0] = r
    nb_deg[:2, 1] = [p-1,p]
    nb_ref[:2, 1] = [r+1,r-2]
    nb_deg[:2, 8] = [p-1,p]
    nb_ref[:2, 8] = [r+1,r-2]
    nb_deg[:2,15] = p
    nb_ref[:2,15] = [r+1,r-2]
    # curves
    nb_ref[0,( 2, 3, 4, 5)] = r+2
    nb_ref[0,( 9,10,11,12)] = r+2
    nb_ref[0,(16,17,18,19)] = r+2
    nb_ref[0,(22,23,24,25)] = r+2
    nb_ref[0,(28,29,30,31)] = r+2
    # lgrge
    nb_deg[0,6]  = p
    nb_deg[0,7]  = p-1
    nb_ref[0,(6,7)]   = r #+1
    nb_deg[0,13] = p
    nb_deg[0,14] = p-1
    nb_ref[0,(13,14)] = r #+1
    nb_deg[0,20] = p
    nb_deg[0,21] = p-1
    nb_ref[0,(20,21)] = r #+1
    nb_deg[0,26] = p
    nb_deg[0,27] = p-1
    nb_ref[0,(26,27)] = r-2
    nb_deg[0,32] = p
    nb_deg[0,33] = p-1
    nb_ref[0,(32,33)] = r-2


completeIGA.refine(nb_ref,nb_deg,additional_knots)
completeIGA._NBPINT[ np.where(completeIGA._ELT_TYPE == 'U00') ] = 6


if EXEMPLE_NB==18:
    icps = manip.get_boundCPindice_wEdges(completeIGA._Nkv,completeIGA._Jpqr,completeIGA._dim, 1,
                                          num_patch=0, offset=1)
    manip.add_displacementBC(completeIGA, completeIGA._indCPbyPatch[0][icps], 3, 0.)


multipleIGA = ddmanip.automaticDomainDecomposition(completeIGA)

dd = []
count = 0
for modeleIGA in multipleIGA:
    dd.append( IGAsubdomain(modeleIGA,count) )
    count += 1


# FETI Analysis
# -------------

print('-'*4,'\nFETI analysis')

# --
# factorization
ti = time.time()
for domain in dd:
    domain.set_stiffnessMATRIX()
    domain.set_couplingMATRIX()
    domain.set_factorizationMATRIX(tol=1.0e-6)
    domain.set_admissibleconstMATRIX()
print(' Local Assembly and Factorization\n (duration : %.2f s).' % (time.time() - ti))


# --
# coarse problem
ti = time.time()
tabWeak = []
for domain in dd:
    interfaceInfos = domain.get_weakInfos()
    tabWeak.append(interfaceInfos)
feti = FETI(tabWeak)

localG = []
for domain in dd:
    localG.append(domain._localG)
feti.set_coarseMATRIX(localG)

localt = []
locale = []
for domain in dd:
    localt.append(domain.compute_condensedRHSvect())
    locale.append(domain.compute_rigidbodyRHSvect())
feti.set_RHSt(localt)
feti.set_RHSe(locale)

feti.solve_coarsePB()
feti.set_projectorP()
print(' Build and resolution of feti coarse problem\n (duration : %.2f s).' % (time.time() - ti))

# --
# shur complement and initialization
ti = time.time()
dofl = feti._dofl
def shurOp(d):
    global feti,dd
    y = np.zeros(feti._dofl,dtype=np.float64)
    for domain in dd:
        ds = feti._assemblyA['%i'%domain._ID].T * d
        ys = domain.evaluatedualshur(ds)
        y[:] += feti._assemblyA['%i'%domain._ID]* ys
    return y
Sd = sp.linalg.LinearOperator((dofl,dofl),matvec=shurOp)

lmbda0 = feti._matrixG.dot(feti._LUfeti.solve(feti._vectE))
print(' Definition of the dual Shur complement and initizalization')
print(' (duration : %.2f s).' % (time.time() - ti))


# --
# preconditionner
ti = time.time()
dofl = feti._dofl

precond_type = 0
if   precond_type is 0:
    # classical Dirichlet - Mortar case
    def precond_localdirichlet(jumpDisp):
        global feti,dd
        dL = np.zeros(feti._dofl,dtype=np.float64)
        for domain in dd:
            us = feti._assemblyA['%i'%domain._ID].T * jumpDisp
            ls = domain.evaluateprimalshur_winvC(us,scaled=False)
            dL[:] += feti._assemblyA['%i'%domain._ID]* ls
        return dL
    Sp = sp.linalg.LinearOperator((dofl,dofl),matvec=precond_localdirichlet)


elif precond_type is 1:
    # generalized Dirichlet
    localC     = []
    localinvdK = []
    for domain in dd:
        localC.append(domain._C2solve)
        domain.set_invdiaKMATRIX()
        localinvdK.append(domain._invDiagK.diagonal())
    feti.set_globalcouplingMATRIX(localC)
    feti.set_globalinvdiaKMATRIX(localinvdK)

    def precond_globaldirichlet_wScaled(jumpDisp):
        global feti,dd
        dL = np.zeros(feti._dofl,dtype=np.float64)
        corrU = feti.evaluatedispcorrection(jumpDisp,scaled=True)
        for domain in dd:
            us = domain._invDiagK*domain._C2solve.T * (feti._assemblyA['%i'%domain._ID].T * corrU)
            fb = domain.evaluateprimalshur(us[domain._itracedisp])
            fs = np.zeros_like(us)
            fs[domain._itracedisp] = fb[:]
            dL+= feti._assemblyA['%i'%domain._ID]* domain._C2solve * domain._invDiagK * fs
        return feti.evaluatedispcorrection(dL,scaled=True)
    SpGloScaled = sp.linalg.LinearOperator((dofl,dofl),matvec=precond_globaldirichlet_wScaled)


elif precond_type is 2:
    # generalized Dirichlet -- no scaling
    localC     = []
    localinvdK = []
    for domain in dd:
        localC.append(domain._C2solve)
    feti.set_globalcouplingMATRIX(localC)

    def precond_globaldirichlet(jumpDisp):
        global feti,dd
        dL = np.zeros(feti._dofl,dtype=np.float64)
        corrU = feti.evaluatedispcorrection(jumpDisp,scaled=False)
        for domain in dd:
            us = domain._C2solve.T * (feti._assemblyA['%i'%domain._ID].T * corrU)
            fb = domain.evaluateprimalshur(us[domain._itracedisp])
            fs = np.zeros_like(us)
            fs[domain._itracedisp] = fb[:]
            dL+= feti._assemblyA['%i'%domain._ID]* domain._C2solve * fs
        return feti.evaluatedispcorrection(dL,scaled=False)
    SpGlo = sp.linalg.LinearOperator((dofl,dofl),matvec=precond_globaldirichlet)


print(' Definition of the preconditionner\n (duration : %.2f s).' % (time.time() - ti))



# --
# Resolution
tcg = time.time()
if precond_type == None:
    lmbda,nbiter= PCPGortho(Sd, feti._vectT, x0=lmbda0, P=feti.projectorP, Pt=feti.projectorPt,
                            tol=1.e-08, maxiter=2000,savetxt=True)
if precond_type == 0:
    lmbda,nbiter= PCPGortho(Sd, feti._vectT, x0=lmbda0, P=feti.projectorP, Pt=feti.projectorPt,
                            M=Sp,
                            tol=1.e-08, maxiter=2000,savetxt=True)
if precond_type == 1:
    lmbda,nbiter= PCPGortho(Sd, feti._vectT, x0=lmbda0, P=feti.projectorP, Pt=feti.projectorPt,
                            M=SpGloScaled,
                            tol=1.e-08, maxiter=2000,savetxt=True)
if precond_type == 2:
    lmbda,nbiter= PCPGortho(Sd, feti._vectT, x0=lmbda0, P=feti.projectorP, Pt=feti.projectorPt,
                            M=SpGlo,
                            tol=1.e-08, maxiter=2000,savetxt=True)
print(' Resolution of feti iterface problem')
print(' (duration : %.2f s, nb iter : %i).\n\n' % (time.time() - tcg,nbiter))



# --
# Postprocessing
print('FETI Results')
alpha = feti._LUfeti.solve( feti._matrixG.T * (feti._vectT - shurOp(lmbda) ) )


count = 0
for domain in dd:
    i = domain._ID
    nbRG   = domain._rigidbodyR.shape[1]
    alphas = alpha[count:count+nbRG]
    count += nbRG
    lmbdas = feti._assemblyA['%i'%i].T * lmbda
    utot   = np.zeros(domain._modeleIGA._nb_dof_tot)
    utot[domain._idof_internal] = domain.evaluate_displacement(lmbdas,alphas)
    idof   = domain._modeleIGA._ind_dof_free[:domain._modeleIGA._nb_dof_free]-1
    SOL,u  = rsol.reconstruction(*domain._modeleIGA.get_inputs4solution(utot[idof]))
    pp.generatevtu(*domain._modeleIGA.get_inputs4postprocVTU(
        'FETI%i'%i,SOL.transpose(),nb_ref=3*np.array([1,1,1]),
        Flag=np.array([True,True,False])))




if True:
    # dof infos
    print('\nInfos:')
    print(' ref',r)
    print(' DOF',completeIGA._nb_dof_free)
    print(' lbd',dofl)



if False:
    # energy
    energy = 0.
    count  = 0
    for domain in dd:
        i = domain._ID
        nbRG   = domain._rigidbodyR.shape[1]
        alphas = alpha[count:count+nbRG]
        count += nbRG
        lmbdas = feti._assemblyA['%i'%i].T * lmbda
        ufree  = domain.evaluate_displacement(lmbdas,alphas)

        energy+= ufree.dot( domain._K2solve.dot( ufree ) )

    energyExact = 1
    if EXEMPLE_NB in [10,19,20]:
        energyExact = 0.11136413923449007
    #print '\nError measure:',np.abs(energy - energyExact)/energyExact
    print('\nEnergy:',energy)


if EXEMPLE_NB in [22] and False:
    # displacement
    print('\nDisp')
    count  = 0
    xi = np.array([1.,0.5,0.])
    for domain in [dd[-1]]:
        i = domain._ID
        nbRG   = domain._rigidbodyR.shape[1]
        alphas = alpha[count:count+nbRG]
        count += nbRG
        lmbdas = feti._assemblyA['%i'%i].T * lmbda
        utot   = np.zeros(domain._modeleIGA._nb_dof_tot)
        utot[domain._idof_internal] = domain.evaluate_displacement(lmbdas,alphas)
        idof   = domain._modeleIGA._ind_dof_free[:domain._modeleIGA._nb_dof_free]-1
        SOL,u  = rsol.reconstruction(*domain._modeleIGA.get_inputs4solution(utot[idof]))
        print(' domain %i:'%i, pp.evaldisp(*domain._modeleIGA.get_inputs4evaldisp(SOL.T,xi) )[0])

