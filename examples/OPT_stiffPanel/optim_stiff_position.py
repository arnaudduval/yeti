

# Python module
import sys
import time
import numpy as np
import nlopt

# IGA module
from preprocessing.igaparametrization import IGAparametrization, \
    IGAmanip as manip
import postprocessing.postproc as pp
import reconstructionSOL as rsol

from preprocessing.igaparametrization import OPTmodelling


# Selection of .INP and .NB file
# ------------------------------
DIRECTORY = 'inputFiles/'
CASE = []
CASE.append('stiffroof')    # ----------- 0

EXEMPLE_NB = 0

FILENAME = DIRECTORY+CASE[EXEMPLE_NB]


# Creation of the IGA object
# --------------------------
modeleIGA = IGAparametrization(filename=FILENAME)

nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

# Plate
nb_deg[:2, 0:2] = 2
nb_ref[:2, 0:2] = 2

# Stiffeners and linked curves
nb_deg[0, np.array([3, 4, 6, 10, 11, 13])-1] = 1
# nb_ref[0,np.array([3,4,6,10,11,13])-1] = 1

modeleIGA.refine(nb_ref, nb_deg)


# --
# Parametrization

# Get CP of top and bottom of hull
cp5 = manip.get_directionCP(modeleIGA, 5, 0, 0)
cp6 = manip.get_directionCP(modeleIGA, 6, 0, 0)

# Get the control points of stiffener located at v=0
cpS = manip.get_boundCPindice_wEdges(modeleIGA._Nkv, modeleIGA._Jpqr,
                                     modeleIGA._dim, 3, num_patch=2)
nb_varS = cpS.size*2


def movestiffener(coords0, igapara, var):
    """
    Shape modification function : move stiffeners control points
    """

    print("var in movestiffener : ", var)

    igapara._COORDS[:, :] = coords0[:, :]

    # Move 1st stiffener
    igapara._COORDS[0, igapara._indCPbyPatch[2][cpS[:]]-1] =\
        0.50*var[:int(nb_varS/2)]
    igapara._COORDS[0, igapara._indCPbyPatch[2][cpS[:]+cpS.size]-1] =\
        0.50*var[:int(nb_varS/2)]

    # Move 2nd stiffener
    igapara._COORDS[0, igapara._indCPbyPatch[9][cpS[:]]-1] =\
        1. - 0.50*var[int(nb_varS/2):]
    igapara._COORDS[0, igapara._indCPbyPatch[9][cpS[:]+cpS.size]-1] =\
        1. - 0.50*var[int(nb_varS/2):]

    return None


movestiffener(modeleIGA._COORDS.copy(), modeleIGA,
              np.array([0., 0.7, 0., 0., 0.7, 0.]))


# --
# Define optimization problem

nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

# plate
nb_deg[:2, :2] = 0
nb_ref[:2, :2] = 3

# stiffeners and linked curves
nb_deg[0, np.array([3, 4, 6, 10, 11, 13])-1] = 1
nb_ref[0, np.array([3, 4, 6, 10, 11, 13])-1] = 6

nb_deg[1, np.array([3, 10])-1] = 2
nb_ref[1, np.array([3, 10])-1] = 2

nb_ref[0, np.array([5, 7, 12, 14])-1] = 6

# Lagrange multipliers
nb_deg[0, np.array([8, 15])-1] = 2
nb_deg[0, np.array([9, 16])-1] = 1
nb_ref[0, np.array([8, 9, 15, 16])-1] = 5


optPB = OPTmodelling(modeleIGA, nb_varS, movestiffener,
                     nb_degreeElevationByDirection=nb_deg,
                     nb_refinementByDirection=nb_ref)


# Modification for analysis (number of integration points along interfaces)
optPB._fineParametrization._NBPINT[
    np.where(optPB._fineParametrization._ELT_TYPE == 'U00')] = 6

# --
# OPTIMIZATION

# Set initial value of design variables
# x0 = np.random.random(nb_varS)*0.9+0.10 #np.ones(nb_varS)*0.4
# x0 = np.concatenate((np.ones(nb_varP)*0.2,np.ones(nb_varS)*0.4))

# x0 = np.ones(nb_varS)*0.0   # 25
x0 = np.array([0., 0.7, 0., 0., 0.7, 0.])

# Set list of patches to be taken into account in volume comptation
listpatch = np.zeros(optPB._coarseParametrization._nb_patch, dtype=np.intp)
listpatch[1] = 1

V0 = 1.10 * optPB.compute_volume(x0*0, listpatch)
c0 = optPB.compute_compliance_discrete(x0)

i = 0
iplt = 0
compprev = 1.
optimHistory = np.zeros((400, 3))


def comp(xC, gradC):
    """
    Compliance and gradient of compliance computation
    (relative to initial value)
    """
    print("xC in comp : ", xC)
    ci = optPB.compute_compliance_discrete(xC)/c0
    if gradC.size > 0:
        global i, listpatch, compprev, iplt
        i += 1
        print('\n--')
        print('Iter %3i' % i)
        # gradC[:] = optPB.compute_gradCompliance_FD(xC, eps=1.e-7)/c0
        # gradC[:] = optPB.compute_gradCompliance_cplgOnly_semiAN(xC,eps=1.e-8)/c0
        gradC[:] = optPB.compute_gradCompliance_cplgOnly_AN(xC)/c0
        # Ci-dessous : OK mais résultats à vérifier
        gradC[:] += optPB.compute_gradCompliance_AN(xC)/c0
        # gradC[:] = optPB.compute_gradCompliance_semiAN(xC)/c0

        if ci <= 1.01*compprev:
            optPB._coarseParametrization._ELT_TYPE_flat \
                = 'U1U0U0U00U00U00U00U4U4U0U00U00U00U00U4U4'
            pp.generatevtu(
                *optPB._coarseParametrization.get_inputs4postprocVTU(
                'opt2coarse%0.2d'%iplt,optPB._coarseParametrization._COORDS-optPB._initialCOORDS,
                nb_ref=np.array([5,5,1]),Flag=np.array([True, False, False])) )

            icps = optPB._coarseParametrization._indCPbyPatch[0]-1
            # np.savetxt('temp/cpsStiffRoof%0.2d.txt'%iplt,
            #           optPB._coarseParametrization._COORDS[:,icps].transpose(),delimiter=',')
            optPB._coarseParametrization.generate_vtk4controlMeshVisu('opt2coarse%0.2d'%iplt,0)

            SOL,u = rsol.reconstruction(
                **optPB._fineParametrization.get_inputs4solution(optPB._save_sol_fine))
            pp.generatevtu(*optPB._fineParametrization.get_inputs4postprocVTU(
                'opt2fine%0.2d'%iplt,  SOL.transpose(),
                nb_ref=np.array([3,3,1]), Flag=np.array([True, True, False])))

            optimHistory[i] = np.array([i,ci,optPB.compute_volume(xC, listpatch)/V0 - 1.])

            compprev = ci
            iplt += 1

    return ci


def vol(xV, gradV):
    """
    Volume and gradient of volume computation
    (relative to initial value)
    """
    print("xV in vol : ", xV)
    global listpatch
    if gradV.size > 0:
        gradV[:] = optPB.compute_gradVolume_AN(xV, listpatch)/V0
    return optPB.compute_volume(xV, listpatch)/V0 - 1.


# minimize = nlopt.opt(nlopt.LD_SLSQP, nb_varP)
minimize = nlopt.opt(nlopt.LD_SLSQP, nb_varS)

minimize.set_min_objective(comp)
minimize.add_inequality_constraint(vol, 1e-5)


minimize.set_ftol_rel(1.0e-06)
minimize.set_xtol_rel(1.0e-06)
minimize.set_maxeval(400)

minimize.set_lower_bounds(0.*np.ones(nb_varS))
minimize.set_upper_bounds(1.*np.ones(nb_varS))

x = minimize.optimize(x0)








# --- BROUILLON
if False:
    # generate optim history plot
    optimHistory[:,2] += 1
    imax = 120
    for i in range(1,imax):
        np.savetxt('temp/optimcvg%i.txt'%i,optimHistory[1:i+1],delimiter=',')



if False:
    # optimal result with out the stiffeners
    xplate = np.array(
        [ 0.09982427,  0.4259896 ,  0.60792969,  0.42598961,  0.09982426,
          0.09982427,  0.15235858,  0.37359429,  0.64076727,  0.37359431,
          0.1523586 ,  0.09982426,  0.4259896 ,  0.37359429,  0.51482169,
          0.65993333,  0.5148217 ,  0.3735943 ,  0.4259896 ,  0.60792969,
          0.64076727,  0.65993333,  0.69664074,  0.65993332,  0.64076727,
          0.60792969,  0.4259896 ,  0.37359429,  0.51482169,  0.65993333,
          0.51482169,  0.3735943 ,  0.4259896 ,  0.09982427,  0.15235857,
          0.37359429,  0.64076727,  0.3735943 ,  0.15235859,  0.09982426,
          0.09982427,  0.4259896 ,  0.60792969,  0.4259896 ,  0.09982427])
    cplate = comp(xplate,np.array([]))

if False:
    # display
    # Plate
    modeleIGA._ELT_TYPE_flat = 'U0U3U0U00U00U00U00U4U4U0U00U00U00U00U4U4'
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'plate',np.zeros_like(modeleIGA._COORDS),nb_ref=np.array([1,1,1]),
        Flag=np.array([False, False, False])) )
    np.savetxt('temp/plateCPs.txt',modeleIGA._COORDS[:,modeleIGA._indCPbyPatch[1]-1].transpose(),
               delimiter=',')
    # Mapping
    modeleIGA._ELT_TYPE_flat = 'U1U0U0U00U00U00U00U4U4U0U00U00U00U00U4U4'
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'mapping',np.zeros_like(modeleIGA._COORDS),nb_ref=np.array([1,1,1]),
        Flag=np.array([False, False, False])) )
    np.savetxt('temp/mappingCPs.txt',modeleIGA._COORDS[:,modeleIGA._indCPbyPatch[0]-1].transpose(),
               delimiter=',')
    modeleIGA.generate_vtk4controlMeshVisu('test',0)

    # Stiffeners
    modeleIGA._ELT_TYPE_flat = 'U0U0U30U00U00U00U00U4U4U30U00U00U00U00U4U4'
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'stiffeners',np.zeros_like(modeleIGA._COORDS),nb_ref=np.array([6,1,1]),
        Flag=np.array([False, False, False])) )
    modeleIGA._ELT_TYPE_flat = 'U0U0U3U00U00U00U00U4U4U3U00U00U00U00U4U4'
    pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
        'stiffenersPara',np.zeros_like(modeleIGA._COORDS),nb_ref=np.array([6,1,1]),
        Flag=np.array([False, False, False])) )
    np.savetxt('temp/stiffCPs.txt',
               modeleIGA._COORDS[:,np.concatenate((modeleIGA._indCPbyPatch[2],
                                                   modeleIGA._indCPbyPatch[9]))-1].transpose(),
               delimiter=',')

