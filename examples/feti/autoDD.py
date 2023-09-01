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

completeIGA.refine(nb_ref, nb_deg, additional_knots)
completeIGA._NBPINT[ np.where(completeIGA._ELT_TYPE == 'U00') ] = 6

multipleIGA = ddmanip.automaticDomainDecomposition(completeIGA)


dd = []
count = 0
for modeleIGA in multipleIGA:
    dd.append( IGAsubdomain(modeleIGA,count) )
    count += 1

for i, domain in enumerate(dd):
    print("Domaine ", i)
    domain.set_stiffnessMATRIX()
    print("Taille K : ", domain._K2solve.shape)
    print("Taille f : ", domain._f2solve.shape)
    domain.set_couplingMATRIX()
    print("Taille C : ", domain._C2solve.shape)
    domain.set_factorizationMATRIX(tol=1.0e-6)
    print("Factorisation")
    print("taille LU : ", domain._LU._LU.shape)
    print("taille P : ", domain._LU._P.shape)
    print("taille R : ", domain._LU.R.shape)
    domain.set_admissibleconstMATRIX()

# Coarse problem
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
