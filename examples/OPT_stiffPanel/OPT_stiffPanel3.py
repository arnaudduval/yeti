

# Python module
import numpy as np
import sys
import time

import scipy.sparse as sp
from scipy.sparse import linalg as sla

#IGA module
from preprocessing.igaparametrization import IGAparametrization, IGAmanip as manip
#from stiffmtrx_embedded import sys_linmat_lindef_static as build_stiffmatrix_emb
#import coupling.penlty_embeddedstruct_rot as cplg
#import bandstorage as bs
import postprocessing.postproc as pp
#import entities2param as tt
import reconstructionSOL as rsol

# Selection of .INP and .NB file
# ------------------------------
DIRECTORY='inputFiles/'
CASE=[]
CASE.append('stiffenedPanelEmbded') # ------ 0

EXEMPLE_NB = 0

FILENAME = DIRECTORY+CASE[EXEMPLE_NB]


# Creation of the IGA object
# --------------------------
modeleIGA = IGAparametrization(filename=FILENAME)

nb_degOPT = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_refOPT = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

nb_degOPT[0,(0,1)] = 1

nb_degOPT[0,(2,4,6)] = 1
nb_refOPT[0,2] = 0

modeleIGA.refine(nb_refOPT,nb_degOPT)


if True:
    # shape modification --> curved panel
    maskcp = modeleIGA._COORDS[0,modeleIGA._indCPbyPatch[0]-1] == 5.
    #modeleIGA._COORDS[2,modeleIGA._indCPbyPatch[0][maskcp]-1] += 0.5
    maskcp = modeleIGA._COORDS[0,modeleIGA._indCPbyPatch[1]-1] == 5.
    #modeleIGA._COORDS[2,modeleIGA._indCPbyPatch[1][maskcp]-1] += 0.5

    maskcp = modeleIGA._COORDS[1,modeleIGA._indCPbyPatch[2]-1] == 0.5
    modeleIGA._COORDS[0,modeleIGA._indCPbyPatch[2][maskcp]-1] += 0.25
    maskcp = modeleIGA._COORDS[1,modeleIGA._indCPbyPatch[4]-1] == 0.5
    modeleIGA._COORDS[0,modeleIGA._indCPbyPatch[4][maskcp]-1] += 0.25

# --
# Parametrization

nb_var = 1
def translationStiffener(coords0,igapara,var):

    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[0,igapara._indCPbyPatch[2]-1] += var[0]
    igapara._COORDS[0,igapara._indCPbyPatch[4]-1] += var[0]

    return None


cptop = manip.get_directionCP(modeleIGA, 4,2,0) - 1
#nb_var = cptop.size
def heightStiffener(coords0,igapara,var):

    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[2,cptop] = var[:]

    return None

#nb_var = 1 + cptop.size
def combine(coords0,igapara,var):
    igapara._COORDS[:,:] = coords0[:,:]
    igapara._COORDS[0,igapara._indCPbyPatch[2]-1] = var[0]
    igapara._COORDS[0,igapara._indCPbyPatch[4]-1] = var[0]
    igapara._COORDS[2,cptop] = var[1:]
    return None

# --
# Definition pb optim

from preprocessing.igaparametrization import OPTmodelling

nb_deg = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,modeleIGA._nb_patch),dtype=np.intp)

p = 3
r = 5
# domains
nb_deg[:2,0] = [p,p]
nb_deg[:2,1] = [p,p]
nb_ref[:2,1] = [r,r]
nb_deg[:2,2] = [p,p]
nb_ref[:2,2] = [r+1,r-2]
# curves
nb_deg[0,(3,4,5,6)] = 1
nb_ref[0,(3,4,5,6)] = r+2
# lgrge
nb_deg[0,(7,8)] = [p,p-1]
nb_ref[0,(7,8)] = r


optPB = OPTmodelling(modeleIGA, nb_var, translationStiffener,
                     nb_degreeElevationByDirection = nb_deg - nb_degOPT,
                     nb_refinementByDirection      = nb_ref - nb_refOPT)

# Modification for analysis
optPB._fineParametrization._NBPINT[ np.where(optPB._fineParametrization._ELT_TYPE == 'U00') ] = 8



# Initialisation:
listpatch   = np.zeros(optPB._coarseParametrization._nb_patch, dtype=np.intp)
listpatch[2]= 1

x0 = np.block([0.3,np.ones(nb_var-1)*0.5])
c0 = optPB.compute_compliance_discrete(x0)
V0 = optPB.compute_volume(x0, listpatch)
i  = 0
saveX = []

import nlopt

def comp(xC,gradC):
    ci = optPB.compute_compliance_discrete(xC)/c0
    if gradC.size>0:
        global i,saveX
        i += 1
        print('\n--')
        print('Iter %3i' % i)
        gradC[:] = optPB.compute_gradCompliance_cplgOnly_AN(xC)/c0
        print(gradC)
        print("===================")
        gradC[:]+= optPB.compute_gradCompliance_AN(xC)/c0
        print(gradC)

        SOL,u = rsol.reconstruction(
            **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
        pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
            'stiff%0.2d'%i,  SOL.transpose(),
            nb_ref=np.array([3,3,0]), Flag=np.array([True,False,False])))

        saveX.append(xC.copy())

    return ci

def vol(xV,gradV):
    global listpatch
    if gradV.size>0:
        gradV[:] = optPB.compute_gradVolume_AN(xV,listpatch)/V0
    return optPB.compute_volume(xV,listpatch)/V0 - 1.




minimize = nlopt.opt( nlopt.LD_SLSQP, nb_var )

minimize.set_min_objective( comp )
#minimize.add_inequality_constraint(  vol, 1e-5 )

minimize.set_ftol_rel(1.0e-10)
minimize.set_xtol_rel(1.0e-10)
minimize.set_maxeval(300)

minimize.set_lower_bounds(-0.20*np.ones(nb_var))
minimize.set_upper_bounds( 0.80*np.ones(nb_var))

x = minimize.optimize( x0 )


