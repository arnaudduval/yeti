# Python module
import numpy as np
import scipy.sparse as sp
import sys
import time

#IGA module
from yeti_iga.solver import PCPGortho
from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.preprocessing.igaparametrization import IGAsubdomain,FETI
from yeti_iga.preprocessing.igaparametrization import IGAmanip    as manip
from yeti_iga.preprocessing.igaparametrization import DDinterface as ddmanip
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp

completeIGA = IGAparametrization(filename='twoplatesDD')

nb_deg = np.zeros((3,completeIGA._nb_patch),dtype=np.intp)
nb_ref = np.zeros((3,completeIGA._nb_patch),dtype=np.intp)
additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}

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

# completeIGA.refine(nb_ref,nb_deg,additional_knots)
completeIGA._NBPINT[ np.where(completeIGA._ELT_TYPE == 'U00') ] = 6

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
    print('A')
    domain.set_stiffnessMATRIX()
    print('B')
    domain.set_couplingMATRIX()
    print('C')
    domain.set_factorizationMATRIX(tol=1.0e-6)
    print('D')
    domain.set_admissibleconstMATRIX()
    print('E')
print(' Local Assembly and Factorization\n (duration : %.2f s).' % (time.time() - ti))
exit()

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


# print(feti._matrixQ)
print(feti._matrixG)
exit()
feti.solve_coarsePB()
exit()
feti.set_projectorP()
print(' Build and resolution of feti coarse problem\n (duration : %.2f s).' % (time.time() - ti))
