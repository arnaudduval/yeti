# Copyright 2018-2019 Thibaut Hirschler

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
Augmentation du degree des fonctions B-Spline.

Le resultat du calcul donne le nouveau vecteur noeuds ainsi
que la matrice qui permet la transformation des points de 
controle.
Cette fonction est la base de l'elevation des degrees des
patch 2D et 3D.

Entree :
 * vecteur noeuds,
 * nombre de points de controle,
 * degree des fcts de nurbs,
 * degrees supplementaires du raffinement.
Sortie : 
 * vecteur noeuds raffine,
 * matrice pour obtenir les nouveaux pts de Controle.
'''



# python modules
import numpy as np



# insertion d'un noeud$
def singleKnot_insertion(uKnots, nb_cp, knot):
  
  # Initialisation
  nb_knots = np.size(uKnots)
  p = nb_knots - nb_cp - 1
  refMatrix = np.zeros((nb_cp+1, nb_cp))
  
  idx = np.searchsorted(uKnots, knot, side='right');
  refKnots = np.insert(uKnots,idx,knot).astype(float)
  
  # Insertion noeud 'knot'
  interval = idx-1
  i = 0
  while i<nb_cp+1:
    if i <= interval-p:
      refMatrix[i,i] = 1
    elif interval-p+1 <= i and i <= interval:
      alpha = (knot - uKnots[i])/(uKnots[i+p] - uKnots[i])
      refMatrix[i,i] = alpha
      refMatrix[i,i-1] = 1-alpha
    else:
      refMatrix[i,i-1] = 1    
    
    i += 1  
  
  return [refKnots, refMatrix]

  

  
  



# Renvoie les composantes uniques et leur nombre d'apparition
def unique_count(vecteur):
  unique, inverse = np.unique(vecteur, return_inverse=True)
  count = np.zeros(len(unique), np.int)
  np.add.at(count, inverse, 1)
  return unique, count


# Decomposition d'un vecteur nodal -> insertion de noeuds jusqu'a continuite C0
def decomposition(uKnots, nb_cp, p):
  
  uniqueKnots, multiplicity = unique_count(uKnots)
  
  stepKnots  = np.copy(uKnots)
  step_nb_cp = nb_cp
  
  Msave = np.eye(nb_cp)
  
  nb_uniqueKnots = np.size(uniqueKnots)
  for num_knots in np.arange(1,nb_uniqueKnots-1):
    nb_insertion = 0
    while nb_insertion < p-multiplicity[num_knots]:
      [stepKnots, stepMatrix] = singleKnot_insertion(stepKnots, step_nb_cp, uniqueKnots[num_knots])
      
      step_nb_cp = np.size(stepMatrix,0)
      Msave = np.dot(stepMatrix, Msave)
      
      nb_insertion += 1
  
  return Msave

  
  
  
  
# Calcul coefficient binomial
def binomial(n, k):
  ''' Binomial coefficient, nCr, aka the "choose" function 
      n! / (r! * (n - r)!)
  '''
  p = 1    
  for i in range(1, min(int(k), int(n - k)) + 1):
      p *= n
      p //= i
      n -= 1
  return float(p)

# Elevation du degree fct bezier (p degree init et r l'augmentation)
def segmentBezierElevation(p, r):
  T = np.zeros((p+r+1,p+1))
  for i in np.arange(0,p+1):
    for j in np.arange(0,r+1):
      T[i+j,i] = binomial(p,i)*binomial(r,j)/binomial(p+r,i+j)
  return T
  
# Construction matrice d'elevation de tous les segments issus de la decomposition
def matrixBezierElevation(uKnots,p,r):

  T = segmentBezierElevation(p,r)
  nb_segment = np.size(np.unique(uKnots))-1
  # np.amax(uKnots).astype('int')-1
  Me = np.zeros((nb_segment*(p+r)+1, nb_segment*p+1))

  for k in np.arange(0,nb_segment):
    ind_row = k*(p+r)
    ind_col = k*p
    Me[ind_row:r+p+ind_row+1, ind_col:ind_col+p+1] = T[:,:]

  return Me

  
  
  
  
  
  
  
# Construction matrice de composition
def composition(uKnots,p,r):
  
  uniqueKnots, multiplicity = unique_count(uKnots)
  nb_uniqueKnots = np.size(uniqueKnots)
  
  uFinal = uniqueKnots[0]*np.ones(p+r+1)
  for k in np.arange(1,nb_uniqueKnots-1):
    uFinal = np.concatenate((uFinal, uniqueKnots[k]*np.ones(r+multiplicity[k])))
  uFinal = np.concatenate((uFinal, uniqueKnots[-1]*np.ones(p+r+1) ))
  
  
  stepKnots  = np.copy(uFinal)
  nb_cpFinal = np.size(stepKnots) - (p+r+1)
  step_nb_cp = nb_cpFinal
  
  Msave = np.eye(nb_cpFinal)
  
  nb_nodes2remove = p-multiplicity
  for num_knots in np.arange(1,nb_uniqueKnots-1):
    nb_insertion = 0
    while nb_insertion < nb_nodes2remove[num_knots]:
      [stepKnots, stepMatrix] = singleKnot_insertion(stepKnots, step_nb_cp, uniqueKnots[num_knots])
      
      step_nb_cp = np.size(stepMatrix,0)
      Msave = np.dot(stepMatrix, Msave)
      
      nb_insertion += 1
  
  MD = Msave
  MDt= np.transpose(MD)
  Mc = np.dot(np.linalg.inv(np.dot(MDt,MD)), MDt)
  
  return Mc, uFinal
  
  
  
  
  
  
  
  
# Construction Matrice totale pour l'elevation du degree
def bsplineDegreeElevation(uKnots,p,r):
  
  nb_cp = np.maximum(np.size(uKnots) - (p+1), 1)
  if r>0:
    Md = decomposition(uKnots, nb_cp, p)
    Me = matrixBezierElevation(uKnots,p,r)
    Mc,uFinal = composition(uKnots,p,r)
    
    Mfinal = np.dot(Mc,np.dot(Me,Md))
  
  else:
    uFinal = uKnots
    Mfinal = np.identity(nb_cp)
  
  return [uFinal, Mfinal]
  
  
 
 
