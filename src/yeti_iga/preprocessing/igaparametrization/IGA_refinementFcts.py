# Copyright 2018-2021 Thibaut Hirschler

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
***************************************************************************************************

CONTENU   Fonctions pour le raffinement par elevation de degree et par insertion de noeuds

Date      29/12/2016
Auteur    Thibaut HIRSCHLER

Fonctionnement :

          [ ... ]

--------------------------------------------------------------------------------------------------
'''



''' IMPORT MODULES '''

# Python modules
# --------------
from time import time
import numpy as np
from scipy.sparse import csc_matrix, eye, kron

# IGA Modules
# -----------
from .bsplineDegreeElevation import bsplineDegreeElevation
from .IGA_manipulation       import get_boundCPindice,find_boundNumber,get_edgeCPindice, \
                                   get_vertexCPindice,get_patchConnectionInfos, \
                                   get_boundCPindice_wEdges




def averageKnot_insertions(uKnots, nb_cp, p):
  '''
  Insertion des noeuds moyens sur les intervalles non-nulles d'un
  vecteur noeud unidimensionnel. Cette fonction est la base du
  raffinement des patch 2D et 3D.

  Entree :
   * vecteur noeuds,
   * nombre de points de controle,
   * degree des fcts de nurbs.
  Sortie :
   * vecteur noeuds raffine,
   * matrice pour obtenir les nouveaux pts de Controle.
  '''

  # Initialisation
  uStep = uKnots
  Mold  = np.identity(nb_cp)

  nb_knots = np.size(uKnots)

  # Insertion noeuds moyens
  nb_insertion = 0
  interval0 = p

  while interval0 < nb_knots-p-1:

    if uKnots[interval0] < uKnots[interval0+1]:

      nb_insertion += 1
      nb_cp_step = nb_cp + nb_insertion
      intervalStep = interval0 + nb_insertion - 1

      Mstep = np.zeros((nb_cp_step, nb_cp_step-1))
      u_moy = 0.5*(uKnots[interval0+1] + uKnots[interval0])

      i = 0
      while i<nb_cp_step:
        if i <= intervalStep-p:
          Mstep[i,i] = 1
        elif intervalStep-p+1 <= i and i <= intervalStep:
          alpha = (u_moy - uStep[i])/(uStep[i+p] - uStep[i])
          Mstep[i,i] = alpha
          Mstep[i,i-1] = 1-alpha
        else:
          Mstep[i,i-1] = 1

        i += 1

      # Actualisation Matrice raffinement et vecteur noeuds
      Mnew = np.dot(Mstep, Mold)
      Mold = Mnew

      uSave = uStep
      uStep = np.zeros(nb_cp_step+p+1);
      uStep[0:intervalStep+1] = uSave[0:intervalStep+1];
      uStep[intervalStep+1] = u_moy;
      uStep[intervalStep+2:]= uSave[intervalStep+1:];


    interval0 += 1

  uFinal = uStep
  Mfinal = csc_matrix(Mnew)

  return [uFinal, Mfinal]





def additionalKnot_insertions(uKnots, nb_cp, p, additional_knots_dir):
  '''
  Insertion de noeuds particuliers contenus dans le vecteur #additional_knots. Cela permet
  d'effectuer un raffinement en des zones precises.

  Entree :
   * vecteur noeuds initial,
   * nombre de points de controle,
   * degree des fcts de nurbs,
   * liste des noeuds a inserer
  Sortie :
   * vecteur noeuds raffine,
   * matrice pour obtenir les nouveaux pts de Controle.
  '''

  # Initialisation
  uStep = uKnots
  Mold  = np.identity(nb_cp)

  nb_knots = np.size(uKnots)

  # Insertion noeuds
  nb_insertion = np.size(additional_knots_dir)
  for num_insertion in np.arange(0,nb_insertion):
    nb_cp_step = nb_cp + num_insertion + 1

    u_ins = additional_knots_dir[num_insertion]
    intervalStep = np.searchsorted(uStep, u_ins, side='right') - 1

    Mstep = np.zeros((nb_cp_step, nb_cp_step-1))

    i = 0
    while i<nb_cp_step:
      if i <= intervalStep-p:
        Mstep[i,i] = 1
      elif intervalStep-p+1 <= i and i <= intervalStep:
        alpha = (u_ins - uStep[i])/(uStep[i+p] - uStep[i])
        Mstep[i,i] = alpha
        Mstep[i,i-1] = 1-alpha
      else:
        Mstep[i,i-1] = 1

      i += 1

    # Actualisation Matrice raffinement et vecteur noeuds
    Mnew = np.dot(Mstep, Mold)
    Mold = Mnew

    uSave = uStep
    uStep = np.zeros(nb_cp_step+p+1);
    uStep[0:intervalStep+1] = uSave[0:intervalStep+1];
    uStep[intervalStep+1] = u_ins;
    uStep[intervalStep+2:]= uSave[intervalStep+1:];

  uFinal = uStep
  Mfinal = csc_matrix(Mnew)

  return [uFinal, Mfinal]





def patchRefinement(nb_degreeElevationByDirection, nb_refinementByDirection, additional_knots,
                    Ukv_init, Nkv_init, Jpqr_init, mcrd_init, dim_patch, num_patch):
  '''
  Raffinement d'un patch Nurbs. Cette methode construit pour un patch,les
  nouveaux vecteurs noeuds ainsi que les matrices de passage de la parametrisation
  creuse vers la parametrisation fine, et ce selon le degree de raffinement par
  direction souhaite.

  Entree :
   * augmentation du degree par direction
   * niveau de raffinement par direction,
   * informations sur les bases Nurbs initiales (vecteurs nodaux, degree etc.),
   * dimension du problem (2D ou 3D),
   * numero du patch en question.
  Sortie :
   * vecteurs nodaux pour la parametrisation fine,
   * matrices de passage parametrisation creuse -> fine par direction.
  '''

  nb_knots_dir = Nkv_init[:, num_patch]
  nb_CP_dir = np.maximum(nb_knots_dir[:] - (Jpqr_init[:,num_patch]+1), np.ones(3))

  knotsFineByDirection = []            # refined knot vectors
  matrixRefinementByDirection = []     # transformation matrix, coarse -> fine parameterization

  # Degree elevation
  knotsDegElevByDirection = np.zeros(dim_patch, dtype='object')
  matrixDegElevByDirection= np.zeros(dim_patch, dtype='object')
  for direction in range(0,dim_patch):
    knotsCoarseThisDirection = Ukv_init[num_patch][direction]

    knotsDegElevByDirection[direction], matrixDegElevByDirection[direction] \
      = bsplineDegreeElevation(knotsCoarseThisDirection, Jpqr_init[direction,num_patch],
                               nb_degreeElevationByDirection[direction])
    matrixDegElevByDirection[direction] = csc_matrix(matrixDegElevByDirection[direction])

  Jpqr_ref = Jpqr_init[:,num_patch] + nb_degreeElevationByDirection


  # Special Knot insertion
  knotsInsertionByDirection       = np.copy(knotsDegElevByDirection)
  matrixKnotsInsertionByDirection = np.copy(matrixDegElevByDirection)
  if num_patch in additional_knots["patches"]:
    for direction in range(0,dim_patch):
      str_dir = '%i' % (direction+1)
      if np.size(additional_knots[str_dir]) > 0:
        u_befIns = knotsInsertionByDirection[direction]
        nb_CP_befIns = np.size(u_befIns) - (Jpqr_ref[direction]+1)

        knotsInsertionByDirection[direction], S_additional_knots \
          = additionalKnot_insertions(u_befIns, nb_CP_befIns, Jpqr_ref[direction],
                                      additional_knots[str_dir])
        matrixKnotsInsertionByDirection[direction] \
          = S_additional_knots * matrixKnotsInsertionByDirection[direction]


  # Average Knot insertion
  for direction in range(0,dim_patch):
    u_level = knotsInsertionByDirection[direction]
    #knotsDegElevByDirection[direction]
    S_level = matrixKnotsInsertionByDirection[direction]
    #csc_matrix(matrixDegElevByDirection[direction])

    for refinementLevel in range(0,nb_refinementByDirection[direction]):
      nb_CP_level = np.size(u_level) - (Jpqr_ref[direction]+1)
      [u_level, S_nextLevel] = averageKnot_insertions(u_level, nb_CP_level, Jpqr_ref[direction])
      S_level  = S_nextLevel * S_level

    knotsFineByDirection.append(u_level)
    matrixRefinementByDirection.append(S_level)

  if dim_patch<2:
    matrixRefinementByDirection.append(eye(nb_CP_dir[1]))
  if dim_patch<3:
    matrixRefinementByDirection.append(eye(nb_CP_dir[2]))

  return  knotsFineByDirection, matrixRefinementByDirection








def parameterizationRefinement(
    nb_degreeElevationByDirection, nb_refinementByDirection, additional_knots, Ukv_init, Nkv_init,
    Jpqr_init, vectWeight_init, mcrd_init, nb_patch, COORDS_init, indCPbyPatch_init,dim):

  '''
  Raffinement de la parameterisation entiere. Cette methode construit les nouveux
  vecteurs nodaux ainsi que les matrices de passage pour la parameterisation fine.
  De plus, le calcul des coordonnes des nouveaux points de controle par patch est
  effectue.
  Les matrices de passage permettent l'identification des points de controle sur
  lesquels s'appliquent les conditions aux limites et les chargements.

  Entree :
   * degree de raffinement par direction,
   * informations sur les bases Nurbs initiales (vecteurs nodaux, degree etc.),
   * dimension du problem (2D ou 3D),
   * nombre de patch,
   * coordonnes des points de controle,
   * indices des points de controle par patch (! Numerotation conventionnee !).
  Sortie :
   * vecteurs nodaux pour la parametrisation fine par patch,
   * matrices de passage parametrisation creuse -> fine par patch.
  '''

  # Appels successifs methode precedente `patchRefinement()'
  knotsFineByPatch = []
  matrixRefinementByPatch = []
  for num_patch in range(0,nb_patch):
    [u_temp, M_temp] = patchRefinement(
      nb_degreeElevationByDirection[:,num_patch], nb_refinementByDirection[:,num_patch],
      additional_knots, Ukv_init, Nkv_init, Jpqr_init, mcrd_init, dim[num_patch], num_patch)
    knotsFineByPatch.append(u_temp)
    matrixRefinementByPatch.append(M_temp)

  Jpqr_ref = Jpqr_init + nb_degreeElevationByDirection

  # Calcul nouveaux coordonnes par patch
  # Contruction matrices de passage
  COORDSfineByPatch = []
  weightfineByPatch = []
  saveS_xi  = []
  saveS_eta = []
  saveS_zeta= []
  nb_cpDbyPatch_fine = []

  for num_patch in range(0,nb_patch):

    nb_cp_thisPatch = np.size(indCPbyPatch_init[num_patch],0)
    Pwi_thisPatch = np.zeros((nb_cp_thisPatch, 4))
    Pwi_thisPatch[:,0] = COORDS_init[0,indCPbyPatch_init[num_patch]-1] \
                         * vectWeight_init[indCPbyPatch_init[num_patch]-1]
    Pwi_thisPatch[:,1] = COORDS_init[1,indCPbyPatch_init[num_patch]-1] \
                         * vectWeight_init[indCPbyPatch_init[num_patch]-1]
    Pwi_thisPatch[:,2] = COORDS_init[2,indCPbyPatch_init[num_patch]-1] \
                         * vectWeight_init[indCPbyPatch_init[num_patch]-1]
    Pwi_thisPatch[:,3] = vectWeight_init[indCPbyPatch_init[num_patch]-1]


    S_temp = matrixRefinementByPatch[num_patch]

    # get number of Control Points pre direction (coarse case)
    nb_knots_coarse = Nkv_init[:, num_patch]
    nb_cpD_coarse = np.maximum(nb_knots_coarse[:] - (Jpqr_init[:,num_patch]+1), np.ones(3))

    # get number of Control Points pre direction (refined case)
    nb_cpD_fine = np.ones(3)
    for direction in range(0,3):
      nb_cpD_fine[direction] = S_temp[direction].shape[0]


    # refinement direction xi
    Id1 = eye(nb_cpD_coarse[1])
    Id2 = eye(nb_cpD_coarse[2])
    S_xi = kron(Id2, kron(Id1, S_temp[0]) )
    Qwi_xi_thisPatch = S_xi * Pwi_thisPatch

    # refinement direction eta
    Id0 = eye(nb_cpD_fine[0])
    Id2 = eye(nb_cpD_coarse[2])
    S_eta = kron(Id2, kron(S_temp[1], Id0) )
    Qwi_xi_eta_thisPatch = S_eta * Qwi_xi_thisPatch

    # refinement direction zeta
    Id0 = eye(nb_cpD_fine[0])
    Id1 = eye(nb_cpD_fine[1])
    S_zeta = kron(S_temp[2], kron(Id1, Id0) )
    Qwi_thisPatch = S_zeta * Qwi_xi_eta_thisPatch

    nb_cp_fine_thisPatch = np.size(Qwi_thisPatch,0)
    Qi_thisPatch = np.zeros((nb_cp_fine_thisPatch,3))
    Qi_thisPatch[:,0] = Qwi_thisPatch[:,0] / Qwi_thisPatch[:,3]
    Qi_thisPatch[:,1] = Qwi_thisPatch[:,1] / Qwi_thisPatch[:,3]
    Qi_thisPatch[:,2] = Qwi_thisPatch[:,2] / Qwi_thisPatch[:,3]


    COORDSfineByPatch.append(Qi_thisPatch)
    weightfineByPatch.append(Qwi_thisPatch[:,3])

    saveS_xi.append(S_xi)
    saveS_eta.append(S_eta)
    saveS_zeta.append(S_zeta)

    nb_cpDbyPatch_fine.append(nb_cpD_fine)

  return [knotsFineByPatch, COORDSfineByPatch, weightfineByPatch, saveS_xi, saveS_eta, saveS_zeta,
          nb_cpDbyPatch_fine]









def build_coordsFine(COORDSfineByPatch, Nkv_fine, Jpqr_fine,
                     dim_init, indCPbyPatch_init, nb_patch, strongcplginfos):
  '''
  Contruction des nouveaux points de controle. Cette methode raccorde les points
  de controle issues du raffinement sur chacun des patchs.

  Entree :
   * coordonnes des points de controle par patch (para. fine),
   * indices des points de controle par patch (para. creuse, ! Numerotation conventionnee !),
   * tableau de definition des elements (para. creuse),
   * nombre d'elements par patch (para. creuse),
   * nombre de patch.
  Sortie :
   * coordonnes des points de controle pour la parameterisation fine,
   * nombre de points de controle pour la parameterisation fine.
  '''

  COORDSfine = COORDSfineByPatch[0]
  for num_patch in range(1,nb_patch):
    idx2rm = np.array([], dtype=np.intp)
    try:
      cplg_thispatch = np.where(strongcplginfos[:,2]==num_patch)[0]
      for i in cplg_thispatch:
        idx2rm = np.append(idx2rm,
                           get_boundCPindice_wEdges(Nkv_fine,Jpqr_fine,dim_init,
                                                    strongcplginfos[i,3],num_patch))
    except:
      None
    nb_cp  = COORDSfineByPatch[num_patch].shape[0]
    idx    = np.setxor1d(np.arange(0,nb_cp),idx2rm)
    COORDSfine = np.append(COORDSfine, COORDSfineByPatch[num_patch][idx], axis=0)

  nb_cp_fine = np.size(COORDSfine, 0)

  return COORDSfine, nb_cp_fine






def updateGeometricModel(nb_degreeElevationByDirection_byPatch,knotsFineByPatch,
                         Jpqr_init,nnode_init,mcrd_init,nb_patch,dim):

  '''
  Mise a jour du modele geometrique. Contruction des donnees definissant la geometrie,
  en accord avec la mise en donnees abqnurbs.

  Entree :
   * vecteurs nodaux pour chacun des patchs (para. fine),
   * degree des Nurbs initiales, nbre de points de controle par element,
   * dimension du pb, nombre de patch
  Sortie :
   * model geometric de la nouvelle parameterisation :
     - Nkv, Ukv, Nijk, Jpqr, weight,
     - elementsByPatch, nb_elem, nnode.
  '''

  # Knot information
  Jpqr_fine = Jpqr_init + nb_degreeElevationByDirection_byPatch
  list_shell_patch = np.where(dim==2)
  Jpqr_fine[2,list_shell_patch] = 0

  nnode_fine = np.prod(Jpqr_fine+1, axis=0, dtype=np.intp)
  Ukv_fine  = knotsFineByPatch
  Nkv_fine  = np.zeros((3,nb_patch), dtype=np.intp)
  for num_patch in range(0,nb_patch):
    for direction in range(0,dim[num_patch]):
      Nkv_fine[direction,num_patch] = Ukv_fine[num_patch][direction].size

  # Element information
  elementsByPatch_fine = np.zeros(nb_patch, dtype=np.intp)
  elem2rmByPatch_fine  = np.zeros(nb_patch, dtype=object)
  nb_elem2rmD          = np.zeros((3,nb_patch), dtype=np.intp)

  num_elem = 0
  for num_patch in range(0,nb_patch):
    list2     = np.array([0], dtype=np.intp)
    elem2rm_2 = np.array([ ], dtype=np.intp)
    if dim[num_patch]>2:
      k2    = np.arange(Jpqr_fine[2,num_patch],Nkv_fine[2,num_patch]-(Jpqr_fine[2,num_patch]+1))
      test2 = Ukv_fine[num_patch][2][k2] != Ukv_fine[num_patch][2][k2+1]
      list2 = Jpqr_fine[2,num_patch]+1 + np.where(test2)[0]
      elem2rm_2                = np.where(np.logical_not(test2))[0]
    nb_elem2rmD[2,num_patch] = elem2rm_2.size

    list1 = np.array([0], dtype=np.intp)
    elem2rm_1 = np.array([ ], dtype=np.intp)
    if dim[num_patch]>1:
      k1    = np.arange(Jpqr_fine[1,num_patch],Nkv_fine[1,num_patch]-(Jpqr_fine[1,num_patch]+1))
      test1 = Ukv_fine[num_patch][1][k1] != Ukv_fine[num_patch][1][k1+1]
      list1 = Jpqr_fine[1,num_patch]+1 + np.where(test1)[0]
      elem2rm_1                = np.where(np.logical_not(test1))[0]
    nb_elem2rmD[1,num_patch] = elem2rm_1.size

    k0    = np.arange(Jpqr_fine[0,num_patch],Nkv_fine[0,num_patch]-(Jpqr_fine[0,num_patch]+1))
    test0 = Ukv_fine[num_patch][0][k0] != Ukv_fine[num_patch][0][k0+1]
    list0 = Jpqr_fine[0,num_patch]+1 + np.where(test0)[0]
    elem2rm_0                = np.where(np.logical_not(test0))[0]
    nb_elem2rmD[0,num_patch] = elem2rm_0.size

    nijk0 = np.tile(list0, list1.size*list2.size)
    nijk1 = np.tile(np.repeat(list1, list0.size),list2.size)
    nijk2 = np.repeat(list2, list0.size*list1.size)
    Nijk_fine_patch = np.array([nijk0,nijk1,nijk2])

    if num_patch==0:
      Nijk_fine = Nijk_fine_patch.copy()
    else:
      Nijk_fine = np.append(Nijk_fine,Nijk_fine_patch,axis=1)

    elementsByPatch_fine[num_patch] = np.size(Nijk_fine_patch,1)
    elem2rmByPatch_fine[num_patch]  = [elem2rm_0,elem2rm_1,elem2rm_2]


  nb_elem_fine = np.sum(elementsByPatch_fine)
  return [Nkv_fine, Ukv_fine, Nijk_fine, elementsByPatch_fine, nb_elem_fine, Jpqr_fine,nnode_fine,
          nb_elem2rmD,elem2rmByPatch_fine]









def updateIndCPbyPatch(COORDSfineByPatch,Nkv_fine,Jpqr_fine,dim_init,strongcplginfos,nb_patch):
  '''
  Creation listes des points de controle par patch pour la parameterisation fine.
  Ces listes permettent en suite de contruite la table de connectivite pour la
  nouvelle parameterisation (i.e. liste des points de controle par element)

  Entree :
   * coordonnes des points de controle par patch,
   * coordonnes des points de controle pour toute la parameterisation fine,
   * nombre de patch.
  Sortie :
   * indices points de controle par patch.
  '''

  count_cp = 0
  indCPbyPatch_fine = np.zeros(nb_patch, dtype=object)
  for num_patch in range(0,nb_patch):
    nb_cp_patch = np.size(COORDSfineByPatch[num_patch],0)
    indCPbyPatch_fine[num_patch] = np.zeros(nb_cp_patch, dtype=np.intp)

    try:
      cplg_thispatch = np.where(strongcplginfos[:,2]==num_patch)[0]
      for i in cplg_thispatch:
        idx2rm = get_boundCPindice_wEdges(Nkv_fine,Jpqr_fine,dim_init,strongcplginfos[i,3],
                                          num_patch,num_orientation=strongcplginfos[i,4])
        idx2add= get_boundCPindice_wEdges(Nkv_fine,Jpqr_fine,dim_init,strongcplginfos[i,1],
                                          strongcplginfos[i,0])
        indCPbyPatch_fine[num_patch][idx2rm] = indCPbyPatch_fine[strongcplginfos[i,0]][idx2add]
    except:
      None
    idx = np.where(indCPbyPatch_fine[num_patch]==0)[0]
    indCPbyPatch_fine[num_patch][idx] = np.arange(1,idx.size+1) + count_cp
    count_cp = np.maximum(count_cp,np.max(indCPbyPatch_fine[num_patch]))

  return indCPbyPatch_fine







def updateConnectivity(Nkv_fine,Jpqr_fine,indCPbyPatch_fine,nb_elem2rmD,elem2rmByPatch,nb_patch):
  '''
  Reconstruction de la table de connectivite pour la parametrisation fine.

  Entree :
   * informations sur le nombre d'element, de points de controle etc.
   * liste des points de controle par patch
  Sortie :
   * table de connectivite : listes des points de controle par elements.
  '''
  IEN_fine = []
  for num_patch in range(0,nb_patch):
    nnodeD   = np.maximum(Jpqr_fine[:,num_patch]+1,1)
    nb_cpD   = np.maximum(Nkv_fine[:,num_patch]-(Jpqr_fine[:,num_patch]+1),       1)
    nb_elemDM= np.maximum(Nkv_fine[:,num_patch]-(Jpqr_fine[:,num_patch]+1)*2 + 1, 1)
    nb_elemD = nb_elemDM - nb_elem2rmD[:,num_patch]

    ien_xi    = np.tile(np.arange(0,nnodeD[0]), (nb_elemD[0],1)) \
                + np.tile(np.vstack(
                  np.delete(np.arange(0,nb_elemDM[0]),elem2rmByPatch[num_patch][0],axis=0)
                ), (1,nnodeD[0]))

    ien_eta   = np.tile(np.arange(0,nnodeD[1]), (nb_elemD[1],1)) \
                + np.tile(np.vstack(
                  np.delete(np.arange(0,nb_elemDM[1]),elem2rmByPatch[num_patch][1],axis=0)
                ), (1,nnodeD[1]))
    ien_eta  *= nb_cpD[0]

    ien_xieta = np.tile(ien_xi,(nb_elemD[1],nnodeD[1])) \
                + np.repeat(np.repeat(ien_eta,nnodeD[0],axis=1), nb_elemD[0], axis=0)

    ien_zeta  = np.tile(np.arange(0,nnodeD[2]), (nb_elemD[2],1)) \
                + np.tile(np.vstack(
                  np.delete(np.arange(0,nb_elemDM[2]),elem2rmByPatch[num_patch][2],axis=0)
                ), (1,nnodeD[2]))
    ien_zeta *= nb_cpD[0]*nb_cpD[1]


    ien_patch = np.tile(ien_xieta,(nb_elemD[2],nnodeD[2])) \
                + np.repeat(
                    np.repeat(ien_zeta,np.prod(nnodeD[:2]),axis=1),np.prod(nb_elemD[:2]),axis=0)

    IEN_fine.append(indCPbyPatch_fine[num_patch][ien_patch[:,::-1]])
  return IEN_fine












def updateBoundaryConditions(
    nb_bc_init,bc_values_init,bc_target_init, indCPbyPatch_init,Nkv_init,Jpqr_init,
    indCPbyPatch_fine,Nkv_fine,Jpqr_fine, nb_patch,dim):

  '''
  Mise a jour des conditions aux limites.

  Entree :
   * toutes les donnees concernant les conditions aux limites (para. creuse),
   * matrices passage parameterisation creuse/fine par patch,
   * indices des points de controle pour les deux parameterisations.
  Sortie :
   * toutes les donnees concernant les conditions aux limites (para. fine).
  '''

  nb_bc_fine     = nb_bc_init
  bc_values_fine = bc_values_init
  bc_target_fine = []
  bc_target_nbelem_fine = []

  num_bc=0
  for bc in bc_target_init:
    find_bc_target_fine = False
    while not find_bc_target_fine:
      i=0
      while i < nb_patch:
        if np.all(np.isin(bc, indCPbyPatch_init[i])):
          break
        i += 1

      # test bc target is the entire patch
      if np.all(np.isin(indCPbyPatch_init[i], bc)):
        bc_target_fine.append(indCPbyPatch_fine[i])
        find_bc_target_fine = True
        break

      # test bc target is face (bound)
      bound = find_boundNumber(bc,i,Nkv_init,Jpqr_init,indCPbyPatch_init)
      if bound:
        idx = get_boundCPindice(Nkv_fine,Jpqr_fine, bound,i)
        bc_target_fine.append(indCPbyPatch_fine[i][idx])
        find_bc_target_fine = True
        break

      # test bc target is edge
      bc_target_fine_thisbc = []
      edges = get_edgeCPindice(Nkv_init, Jpqr_init, dim, i)
      edges_fine = get_edgeCPindice(Nkv_fine, Jpqr_fine, dim, i)
      num_edge = 0
      for edge in edges:
        if np.all(np.isin(indCPbyPatch_init[i][edge],bc)):
          find_bc_target_fine = True
          bc_target_fine_thisbc.append(indCPbyPatch_fine[i][edges_fine[num_edge]])
        num_edge += 1
      if len(bc_target_fine_thisbc)>0:
        bc_target_fine.append(np.ravel(bc_target_fine_thisbc))
        find_bc_target_fine = True
        break

      # test bc target is vertex
      vertex      = get_vertexCPindice(Nkv_init,Jpqr_init,dim,num_patch=i)
      vertex_fine = get_vertexCPindice(Nkv_fine,Jpqr_fine,dim,num_patch=i)
      if np.all(np.isin(bc, indCPbyPatch_init[i][vertex])):
        idx = np.where(np.isin(indCPbyPatch_init[i][vertex], bc))[0]
        bc_target_fine.append(indCPbyPatch_fine[i][vertex_fine][idx])
        find_bc_target_fine = True
        break
      else:
        print('warning: cannot find equivalent bc target of num bc %i --> bc deleted.'%num_bc)
        nb_bc_fine -= 1
        np.delete(bc_values_fine, num_bc, axis=1)
        find_bc_target_fine = True
        break
    num_bc += 1

  bc_target_nbelem_fine = np.zeros(nb_bc_fine, dtype=np.intp)
  for num_bc in range(0,nb_bc_fine):
    bc_target_nbelem_fine[num_bc] = bc_target_fine[num_bc].size

  return nb_bc_fine, bc_values_fine, bc_target_fine, bc_target_nbelem_fine










def updateLoads(nb_load_init, JDLType_init, load_target_nbelem_init, indDLoad_init, ADLMAG_init,
                indCPbyPatch_init, Nkv_init, Jpqr_init, Nijk_init, elementsByPatch_init, Ukv_init,
                indCPbyPatch_fine, Nkv_fine, Jpqr_fine, Nijk_fine, elementsByPatch_fine, Ukv_fine,
                dim, nb_patch):

  '''
  Mise a jour des chargements.

  Entree :
   * toutes les donnees concernant les chargements (para. creuse),
   * matrices passage parameterisation creuse/fine par patch,
   * indices des points de controle pour les deux parameterisations.
  Sortie :
   * toutes les donnees concernant les chargements (para. fine).
  '''

  nb_load_fine = nb_load_init
  JDLType_fine = JDLType_init
  ADLMAG_fine  = ADLMAG_init
  load_target_nbelem_fine = np.zeros(nb_load_fine, dtype=np.intp)
  indDLoad_fine           = np.zeros(nb_load_fine, dtype=object)

  num_load  = np.arange(0,nb_load_fine, dtype=np.intp)
  num_dload = np.nonzero(np.floor_divide(JDLType_fine,10))[0]
  num_cload = np.setdiff1d(num_load,num_dload)

  #*Concentrated load
  for n in num_cload:
    cps = indDLoad_init[n]
    cps_fine = []
    # find patches
    list_patch = np.zeros_like(cps, dtype=np.intp)
    num_patch  = 0
    count = 0
    while num_patch < nb_patch and count<list_patch.size:
      test = np.isin(cps,indCPbyPatch_init[num_patch])
      if np.any(test):
        list_patch[test] = num_patch
        count += 1
      num_patch += 1

    # get new points
    i = 0
    for cp in cps:
      vertex      = get_vertexCPindice(Nkv_init,Jpqr_init,dim,list_patch[i])
      vertex_fine = get_vertexCPindice(Nkv_fine,Jpqr_fine,dim,list_patch[i])
      test        = np.where(indCPbyPatch_init[list_patch[i]][vertex]==cp)
      cps_fine.append(indCPbyPatch_fine[list_patch[i]][vertex_fine][test])
      i+=1

    # update
    cps_fine = np.ravel(cps_fine)
    if cps_fine.size==0:
      print(' warning: cannot find equivalence of concentrated load number %i. Deleted.' % n)
      np.delete(JDLType_fine, n)
      np.delete(ADLMAG_fine,  n)
      np.delete(indDLoad_fine,n)
      np.delete(load_target_nbelem_fine, n)
      nb_load_fine -= 1
    else:
      indDLoad_fine[n] = cps_fine
      load_target_nbelem_fine[n] = cps_fine.size

  #*Cas distributed load
  for n in num_dload:
    elems   = np.array(indDLoad_init[n])
    list_patch = np.searchsorted(np.cumsum(elementsByPatch_init),elems)

    if JDLType_init[n]>99:
      NumFace = 0 # body load
    else:
      NumFace = np.floor_divide(JDLType_init[n],10)

    i = 0
    indDLoad_thisLoad = np.array([], dtype=np.intp)
    for elem in elems:
      ukv_bounds = np.zeros((3,2))
      for j in np.arange(0,dim[list_patch[i]]):
        ukv_bounds[j,:]=Ukv_init[list_patch[i]][j][Nijk_init[j,elem-1]-1:Nijk_init[j,elem-1]+1]

      patch_bounds = np.insert(np.cumsum(elementsByPatch_fine),0,0)

      # body force
      if NumFace==0:
        elems_fine = np.arange(0,elementsByPatch_fine[list_patch[i]])
        for xi in np.arange(0,dim[list_patch[i]]):
          test4nijk_m  = np.where(Ukv_fine[list_patch[i]][xi] == ukv_bounds[xi,0])[0][-1]
          elems_fine_m = np.where(
            Nijk_fine[xi,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
            >= test4nijk_m+1)[0]
          elems_fine = np.intersect1d(elems_fine,elems_fine_m)

          test4nijk_M  = np.where(Ukv_fine[list_patch[i]][xi] == ukv_bounds[xi,1])[0][0]
          elems_fine_M = np.where(
            Nijk_fine[xi,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
            < test4nijk_M+1)[0]
          elems_fine = np.intersect1d(elems_fine,elems_fine_M)

      # pressure on face 1
      if NumFace==1:
        test4nijk  = np.where(Ukv_fine[list_patch[i]][0] == ukv_bounds[0,0])[0][-1]
        elems_fine = np.where(
          Nijk_fine[0,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
          == test4nijk+1)[0]
      # pressure on face 2
      elif NumFace==2:
        test4nijk  = np.where(Ukv_fine[list_patch[i]][0] == ukv_bounds[0,1])[0][0]
        elems_fine = np.where(
          Nijk_fine[0,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
          == test4nijk)[0]
      # pressure on face 3
      elif NumFace==3:
        test4nijk  = np.where(Ukv_fine[list_patch[i]][1] == ukv_bounds[1,0])[0][-1]
        elems_fine = np.where(
          Nijk_fine[1,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
          == test4nijk+1)[0]
      # pressure on face 4
      elif NumFace==4:
        test4nijk  = np.where(Ukv_fine[list_patch[i]][1] == ukv_bounds[1,1])[0][0]

        bound0 = np.where(Ukv_fine[list_patch[i]][0] == ukv_bounds[0,0])[0][-1]
        bound1 = np.where(Ukv_fine[list_patch[i]][0] == ukv_bounds[0,1])[0][0]

        elems_fine = np.where(np.all(np.array(
          [Nijk_fine[1,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]] == test4nijk,
           Nijk_fine[0,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]] >= bound0+1 ,
           Nijk_fine[0,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]] <= bound1 ])
                                     ,axis=0) )[0]

      # pressure on face 5
      elif NumFace==5:
        if dim[list_patch[i]]>2:
          test4nijk  = np.where(Ukv_fine[list_patch[i]][2] == ukv_bounds[2,0])[0][-1]
          elems_fine = np.where(
            Nijk_fine[2,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
            == test4nijk+1)[0]
      # pressure on face 6
      elif NumFace==6:
        if dim[list_patch[i]]>2:
          test4nijk  = np.where(Ukv_fine[list_patch[i]][2] == ukv_bounds[2,1])[0][0]
          elems_fine = np.where(
            Nijk_fine[2,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
            == test4nijk)[0]
      # face 5 and 6 for shells
      if NumFace>4 and  dim[list_patch[i]]==2:
          elems_fine = np.arange(0,elementsByPatch_fine[list_patch[i]])
          for xi in np.arange(0,2):
            test4nijk_m  = np.where(Ukv_fine[list_patch[i]][xi] == ukv_bounds[xi,0])[0][-1]
            elems_fine_m = np.where(
              Nijk_fine[xi,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
              >= test4nijk_m+1)[0]
            elems_fine = np.intersect1d(elems_fine,elems_fine_m)

            test4nijk_M  = np.where(Ukv_fine[list_patch[i]][xi] == ukv_bounds[xi,1])[0][0]
            elems_fine_M = np.where(
              Nijk_fine[xi,patch_bounds[list_patch[i]]:patch_bounds[list_patch[i]+1]]
              < test4nijk_M+1)[0]
            elems_fine = np.intersect1d(elems_fine,elems_fine_M)

      indDLoad_thisLoad = np.append(indDLoad_thisLoad, elems_fine+patch_bounds[list_patch[i]]+1)
      i += 1

    indDLoad_fine[n]           = indDLoad_thisLoad
    load_target_nbelem_fine[n] = indDLoad_thisLoad.size



  return nb_load_fine, JDLType_fine, load_target_nbelem_fine, indDLoad_fine, ADLMAG_fine









'''
---------------------------------------------------------------------------------------------------
                                                                                                '''
def iga_refinement(
    nb_degreeElevationByDirection, nb_refinementByDirection, mechanicalSettings, geometricSettings,
    indCPbyPatch_coarse,
    additional_knots = {"patches":np.array([]),"1":np.array([]),"2":np.array([]),"3":np.array([])}
):

  '''
  Methode principale, utilisant toutes les precedentes afin de raffiner par insertion
  de noeuds, une parameterisation et de recreer toute les connectivites.
  '''


  # Lecture des donnees initiales :
  # -----------------------------

  # Mechanical properties:
  #   parameters
  NBPINT_coarse = mechanicalSettings[0][1]
  mcrd_coarse   = mechanicalSettings[0][3]

  #   boundary conditions
  bc_target_coarse = mechanicalSettings[1][0]
  bc_values_coarse = mechanicalSettings[1][1]
  bc_target_nbelem_coarse = mechanicalSettings[1][2]
  nb_bc_coarse     = mechanicalSettings[1][3]

  #   loads
  indDLoad_coarse = mechanicalSettings[2][0]
  JDLType_coarse  = mechanicalSettings[2][1]
  ADLMAG_coarse   = mechanicalSettings[2][2]
  load_target_nbelem_coarse = mechanicalSettings[2][3]
  nb_load_coarse  = mechanicalSettings[2][4]

  #   nodes info
  COORDS_coarse   = mechanicalSettings[3][0]
  nb_cp_coarse    = mechanicalSettings[3][1]

  #   elements info
  IEN_coarse      = mechanicalSettings[4]
  #nb_elem_coarse = mechanicalSettings[4][1]


  # Geometric properties:
  dim_coarse      = geometricSettings[0]
  Nkv_coarse      = geometricSettings[1]
  Ukv_coarse      = geometricSettings[2]
  Jpqr_coarse     = geometricSettings[3]
  Nijk_coarse     = geometricSettings[4]
  weight_coarse   = geometricSettings[5]
  elementsByPatch_coarse = geometricSettings[6]
  nb_patch_coarse = geometricSettings[7]
  nnode_coarse    = geometricSettings[8]
  nb_elem_coarse  = geometricSettings[9]



  # Initialisation :
  # --------------

  new_mechanicalSettings = mechanicalSettings
  new_geometricSettings  = geometricSettings

  indCPbyPatch_fine = indCPbyPatch_coarse
  new_transformationMatrices = []
  for num_patch in range(0,nb_patch_coarse):
    Identity_thisPatch = eye(np.size(indCPbyPatch_coarse[num_patch],0))
    new_transformationMatrices.append([Identity_thisPatch, Identity_thisPatch, Identity_thisPatch])


  # Procedure de raffinement :
  # ------------------------
  if (np.sum(nb_degreeElevationByDirection) + np.sum(nb_refinementByDirection)
      +additional_knots['1'].size + additional_knots['2'].size + additional_knots['3'].size) > 0:

    # -------
    vectWeight_coarse = np.zeros(nb_cp_coarse)
    Jelem = 0
    for num_patch in range(0,nb_patch_coarse):
      for num_elem in range(0,elementsByPatch_coarse[num_patch]):
        vectWeight_coarse[IEN_coarse[num_patch][num_elem,:]-1 ]=weight_coarse[Jelem][:]
        Jelem += 1
    # -------


    # 1. Raffinement par patch .................................................................. #
    [knotsFineByPatch, COORDSfineByPatch, weightfineByPatch, saveS_xi, saveS_eta, saveS_zeta,
     nb_cpDbyPatch_fine] = parameterizationRefinement(
       nb_degreeElevationByDirection, nb_refinementByDirection, additional_knots, Ukv_coarse,
       Nkv_coarse, Jpqr_coarse, vectWeight_coarse, mcrd_coarse, nb_patch_coarse, COORDS_coarse,
       indCPbyPatch_coarse,dim_coarse)

    # 2. Construction donnees geometriques ...................................................... #
    [Nkv_fine, Ukv_fine, Nijk_fine, elementsByPatch_fine, nb_elem_fine, Jpqr_fine, nnode_fine,
     nb_elem2rmD,elem2rmByPatch_fine] = updateGeometricModel(
       nb_degreeElevationByDirection, knotsFineByPatch, Jpqr_coarse, nnode_coarse, mcrd_coarse,
       nb_patch_coarse,dim_coarse)

    # 3. Construction Polygone raffine .......................................................... #
    tab = get_patchConnectionInfos(Nkv_coarse,Jpqr_coarse,dim_coarse,indCPbyPatch_coarse,
                                   nb_patch_coarse)
    COORDSfine, nb_cp_fine = build_coordsFine(
      COORDSfineByPatch, Nkv_fine, Jpqr_fine, dim_coarse,indCPbyPatch_coarse, nb_patch_coarse, tab)

    # 4. Construction connectivite .............................................................. #
    indCPbyPatch_fine = updateIndCPbyPatch(COORDSfineByPatch,Nkv_fine,Jpqr_fine,dim_coarse,tab,
                                           nb_patch_coarse)
    IEN_fine = updateConnectivity(Nkv_fine,Jpqr_fine,indCPbyPatch_fine,
                                  nb_elem2rmD,elem2rmByPatch_fine,nb_patch_coarse)

    # -------
    vectWeight_fine = np.zeros(nb_cp_fine)
    for num_patch in range(0,nb_patch_coarse):
      vectWeight_fine[ indCPbyPatch_fine[num_patch][:]-1 ] = weightfineByPatch[num_patch][:]

    weight_fine = []
    for num_patch in np.arange(0,nb_patch_coarse):
      for num_elem in np.arange(0,elementsByPatch_fine[num_patch]):
        weight_fine.append(vectWeight_fine[IEN_fine[num_patch][num_elem]-1])
    # -------


    # 5. Mise a jour conditions aux limites ..................................................... #
    [nb_bc_fine, bc_values_fine, bc_target_fine, bc_target_nbelem_fine] \
      = updateBoundaryConditions(
        nb_bc_coarse, bc_values_coarse, bc_target_coarse, indCPbyPatch_coarse,Nkv_coarse,
        Jpqr_coarse,indCPbyPatch_fine,Nkv_fine, Jpqr_fine, nb_patch_coarse, dim_coarse)

    # 6. Mise a jour des chargements ............................................................ #
    [nb_load_fine, JDLType_fine, load_target_nbelem_fine, indDLoad_fine,ADLMAG_fine] \
      = updateLoads(
        nb_load_coarse, JDLType_coarse, load_target_nbelem_coarse,indDLoad_coarse,ADLMAG_coarse,
        indCPbyPatch_coarse,Nkv_coarse,Jpqr_coarse,Nijk_coarse,elementsByPatch_coarse,Ukv_coarse,
        indCPbyPatch_fine, Nkv_fine, Jpqr_fine, Nijk_fine, elementsByPatch_fine, Ukv_fine,
        dim_coarse,nb_patch_coarse)


    # ........................................................................................... #
    # FIN - Recapulatifs des donnees raffinees .................................................. #
    # ........................................................................................... #

    # Update mechanicalSettings
    #NBPINT_fine = np.ceil( (np.max(Jpqr_fine) + 1)/2.0 ) ** dim_coarse
    NBPINT_fine = (np.max(Jpqr_fine,axis=0)+1)**dim_coarse


    new_mechanicalSettings[0][1] = NBPINT_fine
    new_mechanicalSettings[1] = [bc_target_fine, bc_values_fine, bc_target_nbelem_fine, nb_bc_fine]
    new_mechanicalSettings[2] = [indDLoad_fine, JDLType_fine, ADLMAG_fine, load_target_nbelem_fine,
                                 nb_load_fine]
    new_mechanicalSettings[3] = [COORDSfine.transpose(), nb_cp_fine]
    new_mechanicalSettings[4] = IEN_fine


    # Update geometricSettings
    new_geometricSettings[1] = Nkv_fine
    new_geometricSettings[2] = Ukv_fine
    new_geometricSettings[3] = Jpqr_fine
    new_geometricSettings[4] = Nijk_fine
    new_geometricSettings[5] = weight_fine
    new_geometricSettings[6] = elementsByPatch_fine
    new_geometricSettings[8] = nnode_fine
    new_geometricSettings[9] = nb_elem_fine


    # Update transformation matrices
    new_transformationMatrices = []
    for num_patch in range(0,nb_patch_coarse):
      new_transformationMatrices.append(
        [saveS_xi[num_patch], saveS_eta[num_patch], saveS_zeta[num_patch] ])

    print(" Refinement has been successfully done.")


    #del load_target_nbelem_fine

  return [new_mechanicalSettings, new_geometricSettings, indCPbyPatch_fine,
          new_transformationMatrices]







