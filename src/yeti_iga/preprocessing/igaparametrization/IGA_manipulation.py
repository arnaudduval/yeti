# Copyright 2018-2020 Thibaut Hirschler
# Copyright 2020 Arnaud Duval

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

"""
Created on Fri Feb 23 2018

@author: thirschler

CONTAINS:
Utility functions to apply boundary conditions and for multi-patch coupling.
"""

import sys

# from .IGA_parametrization import IGAparametrization
import numpy as np
import scipy.sparse as sp
from scipy.linalg import blas

from ...fitting.interpolate import buildgrevinterpolmat, getgrevabscphysicalcoords
from .bsplineDegreeElevation import decomposition,decomposition_sparse,localExtraction1Dunique


def get_boundCPindice(Nkv, Jpqr, num_bound, num_patch=0, offset=0):
    """Get indices from boundary face control points of a given patch.

    Returns local indices of control points controlling the boundary indexed at
    `num_bound` of the patch indexed at `num_patch`. The `offset` argument
    enables to return the indices of control points forming a following row.


    Parameters
    ----------
    Nkv : list of int
        Number of knots by direction.
    Jpqr : list of int
        Degree by direction.
    num_bound : int
        Index of the boundary face.
        Must be ``1``, ``2``, ``3``, ``4``, ``5`` or ``6``.
    num_patch : int, optional
        Index of the patch. The default is 0.
    offset : int, optional
        Offset to get the indices of a following row of control points.
        The default is 0.

    Returns
    -------
    ind_face : list of int
        List of local control points indices on the boundary.

    """
    nb_cpPatchD = np.intp(
        np.maximum(Nkv[:, num_patch] - (Jpqr[:, num_patch] + 1), 1))

    edgeXi = np.arange(0, nb_cpPatchD[0], 1)
    edgeEta = np.arange(0, nb_cpPatchD[1], 1) * np.prod(nb_cpPatchD[0:1])
    edgeZeta = np.arange(0, nb_cpPatchD[2], 1) * np.prod(nb_cpPatchD[0:2])

    # Face s.t. xi = 0
    if num_bound == 1:
        ind_face1 = \
            np.kron(np.ones(nb_cpPatchD[2]), edgeEta + edgeXi[offset]) + \
            np.kron(edgeZeta, np.ones(nb_cpPatchD[1]))
        return ind_face1.astype(int)
    # Face s.t. xi = 1
    if num_bound == 2:
        ind_face2 = \
            np.kron(np.ones(nb_cpPatchD[2]), edgeEta + edgeXi[-1 - offset]) + \
            np.kron(edgeZeta, np.ones(nb_cpPatchD[1]))
        return ind_face2.astype(int)
    # Face s.t. eta = 0
    if num_bound == 3:
        ind_face3 = \
            np.kron(np.ones(nb_cpPatchD[2]), edgeXi + edgeEta[offset]) + \
            np.kron(edgeZeta, np.ones(nb_cpPatchD[0]))
        return ind_face3.astype(int)
    # Face s.t. eta = 1
    if num_bound == 4:
        ind_face4 = \
            np.kron(np.ones(nb_cpPatchD[2]), edgeXi + edgeEta[-1 - offset]) + \
            np.kron(edgeZeta, np.ones(nb_cpPatchD[0]))
        return ind_face4.astype(int)
    # Face s.t. zeta = 0
    if num_bound == 5:
        ind_face5 = \
            np.kron(np.ones(nb_cpPatchD[1]), edgeXi + edgeZeta[offset]) + \
            np.kron(edgeEta, np.ones(nb_cpPatchD[0]))
        return ind_face5.astype(int)
    # Face s.t. zeta = 1
    if num_bound == 6:
        ind_face6 = \
            np.kron(np.ones(nb_cpPatchD[1]),
                    edgeXi + edgeZeta[-1 - offset]) + \
            np.kron(edgeEta, np.ones(nb_cpPatchD[0]))
        return ind_face6.astype(int)


def get_boundCPindice_wEdges(Nkv, Jpqr, dim, num_bound, num_patch=0,
                             offset=0, num_orientation=0) -> np.ndarray:
    """Get indices from boundary face control points of a given patch.

    Returns local indices of control points controlling the boundary indexed at
    `num_bound` of the patch indexed at `num_patch`. The `offset` argument
    enables to return the indices of control points forming a following row.

    *Special cases:*
        If boundary number is between ``100`` and ``1000``: returns edges

        If boundary number is greater than ``1000``: returns vertices


    Parameters
    ----------
    Nkv : list of int
        Number of knots by direction.
    Jpqr : list of int
        Degree by direction.
    num_bound : int
        Index of the boundary face.
        Must be ``1``, ``2``, ``3``, ``4``, ``5`` or ``6``.
    num_patch : int, optional
        Index of the patch (python numbering). The default is 0.
    offset : int, optional
        Offset to get the indices of a following row of control points.
        The default is 0.
    num_orientation : int, optional
        Local ordering of the output array.\n
        It should be between ``0`` and ``7``, and defines the local directions as:\n
          0: ( 1, 2), 1: (-1, 2), 2: ( 1,-2), 3: (-1,-2),\n
          4: ( 2, 1), 5: (-2, 1), 6: ( 2,-1), 7: (-2,-1).
        For edges, if ``orientation`` is an odd number then the list is flipped, else
        it is outputed in the natural order.
        This argument has no effects for vertices.

        This option is useful when coupling strongly patches that have
        interface faces/edges with different orientations.

    Returns
    -------
    ind_face : list of int
        List of local control points indices on the boundary.
    edge : list of int
        List of local control points indices on the edge.
    vertex : list of int
        List of local control points indices on the vertex.

    """
    if not num_orientation in np.arange(8):
        raise ValueError("Wrong value for argument ``Orientation``. "\
                         "See description of the function.")
    localorientations = [
        ( 1, 2), (-1, 2), ( 1,-2), (-1,-2),
        ( 2, 1), (-2, 1), ( 2,-1), (-2,-1)]
    orientation = localorientations[num_orientation]

    nb_cpPatchD = np.intp(
        np.maximum(Nkv[:, num_patch] - (Jpqr[:, num_patch] + 1), 1))

    edgeXi = np.arange(0, nb_cpPatchD[0], 1)
    edgeEta = np.arange(0, nb_cpPatchD[1], 1) * np.prod(nb_cpPatchD[0:1])
    edgeZeta = np.arange(0, nb_cpPatchD[2], 1) * np.prod(nb_cpPatchD[0:2])

    # Face s.t. xi = 0
    if num_bound == 1:
        l1 = edgeEta + edgeXi[offset] if abs(orientation[0])==1 else edgeZeta
        l2 = edgeZeta if abs(orientation[1])==2 else edgeEta + edgeXi[offset]
        if orientation[0]<0: l1 = np.flip(l1)
        if orientation[1]<0: l2 = np.flip(l2)
        ind_face1 = np.kron(np.ones(l2.size), l1) + np.kron(l2, np.ones(l1.size))
        #ind_face1 = \
        #    np.kron(np.ones(nb_cpPatchD[2]), edgeEta + edgeXi[offset]) + \
        #    np.kron(edgeZeta, np.ones(nb_cpPatchD[1]))
        return ind_face1.astype(int)
    # Face s.t. xi = 1
    if num_bound == 2:
        l1 = edgeEta + edgeXi[-1 - offset] if abs(orientation[0])==1 else edgeZeta
        l2 = edgeZeta if abs(orientation[1])==2 else edgeEta + edgeXi[-1 - offset]
        if orientation[0]<0: l1 = np.flip(l1)
        if orientation[1]<0: l2 = np.flip(l2)
        ind_face2 = np.kron(np.ones(l2.size), l1) + np.kron(l2, np.ones(l1.size))
        # ind_face2 = \
        #     np.kron(np.ones(nb_cpPatchD[2]), edgeEta + edgeXi[-1 - offset]) + \
        #     np.kron(edgeZeta, np.ones(nb_cpPatchD[1]))
        return ind_face2.astype(int)
    # Face s.t. eta = 0
    if num_bound == 3:
        l1 = edgeXi + edgeEta[offset] if abs(orientation[0])==1 else edgeZeta
        l2 = edgeZeta if abs(orientation[1])==2 else edgeXi + edgeEta[offset]
        if orientation[0]<0: l1 = np.flip(l1)
        if orientation[1]<0: l2 = np.flip(l2)
        ind_face3 = np.kron(np.ones(l2.size), l1) + np.kron(l2, np.ones(l1.size))
        # ind_face3 = \
        #     np.kron(np.ones(nb_cpPatchD[2]), edgeXi + edgeEta[offset]) + \
        #     np.kron(edgeZeta, np.ones(nb_cpPatchD[0]))
        return ind_face3.astype(int)
    # Face s.t. eta = 1
    if num_bound == 4:
        l1 = edgeXi + edgeEta[-1 - offset] if abs(orientation[0])==1 else edgeZeta
        l2 = edgeZeta if abs(orientation[1])==2 else edgeXi + edgeEta[-1 - offset]
        if orientation[0]<0: l1 = np.flip(l1)
        if orientation[1]<0: l2 = np.flip(l2)
        ind_face4 = np.kron(np.ones(l2.size), l1) + np.kron(l2, np.ones(l1.size))
        # ind_face4 = \
        #     np.kron(np.ones(nb_cpPatchD[2]), edgeXi + edgeEta[-1 - offset]) + \
        #     np.kron(edgeZeta, np.ones(nb_cpPatchD[0]))
        return ind_face4.astype(int)
    # Face s.t. zeta = 0
    if num_bound == 5:
        l1 = edgeXi + edgeZeta[offset] if abs(orientation[0])==1 else edgeEta
        l2 = edgeEta if abs(orientation[1])==2 else edgeXi + edgeZeta[offset]
        if orientation[0]<0: l1 = np.flip(l1)
        if orientation[1]<0: l2 = np.flip(l2)
        ind_face5 = np.kron(np.ones(l2.size), l1) + np.kron(l2, np.ones(l1.size))
        # ind_face5 = \
        #     np.kron(np.ones(nb_cpPatchD[1]), edgeXi + edgeZeta[offset]) + \
        #     np.kron(edgeEta, np.ones(nb_cpPatchD[0]))
        return ind_face5.astype(int)
    # Face s.t. zeta = 1
    if num_bound == 6:
        l1 = edgeXi + edgeZeta[-1 - offset] if abs(orientation[0])==1 else edgeEta
        l2 = edgeEta if abs(orientation[1])==2 else edgeXi + edgeZeta[-1 - offset]
        if orientation[0]<0: l1 = np.flip(l1)
        if orientation[1]<0: l2 = np.flip(l2)
        ind_face6 = np.kron(np.ones(l2.size), l1) + np.kron(l2, np.ones(l1.size))
        # ind_face6 = \
        #     np.kron(np.ones(nb_cpPatchD[1]),
        #             edgeXi + edgeZeta[-1 - offset]) + \
        #     np.kron(edgeEta, np.ones(nb_cpPatchD[0]))
        return ind_face6.astype(int)
    # Edges
    if num_bound > 100 and num_bound < 1000:
        edges = get_edgeCPindice(Nkv, Jpqr, dim, num_patch=num_patch)
        e = np.mod(num_bound, 100)
        return edges[e - 1] if orientation[0]>0 else np.flip(edges[e - 1])

    # Vertices
    if num_bound > 1000:
        vertex = get_vertexCPindice(Nkv, Jpqr, dim, num_patch=num_patch)
        v = np.mod(num_bound, 1000)
        return vertex[v - 1]


def get_boundElement(igaPara, num_bound, num_patch=0):
    """Get boundary elements indices.

    Parameters
    ----------
    igaPara : IGA_parametrization
        NURBS parametrisation for IGA..
    num_bound : int
        Index of the boundary face.
        Must be ``1``, ``2``, ``3``, ``4``, ``5`` or ``6``.
    num_patch : int, optional
        Index of the patch. The default is 0.

    Returns
    -------
    elems : list
        List of elements indices.

    """
    patch_bounds = np.insert(np.cumsum(igaPara._elementsByPatch), 0, 0)
    elems = np.array([], dtype=np.intp)
    # Face s.t. xi = 0
    if num_bound == 1:
        elems = np.where(
            igaPara._Nijk[0, patch_bounds[num_patch]:patch_bounds[num_patch+1]]
            == igaPara._Jpqr[0, num_patch] + 1)[0]
    # Face s.t. xi = 1
    if num_bound == 2:
        elems = np.where(
            igaPara._Nijk[0, patch_bounds[num_patch]:patch_bounds[num_patch+1]]
            == igaPara._Nkv[0, num_patch] - (igaPara._Jpqr[0, num_patch]+1))[0]
    # Face s.t. eta = 0
    if num_bound == 3:
        elems = np.where(
            igaPara._Nijk[1, patch_bounds[num_patch]:patch_bounds[num_patch+1]]
            == igaPara._Jpqr[1, num_patch] + 1)[0]
    # Face s.t. eta = 1
    if num_bound == 4:
        elems = np.where(
            igaPara._Nijk[1, patch_bounds[num_patch]:patch_bounds[num_patch+1]]
            == igaPara._Nkv[1, num_patch] - (igaPara._Jpqr[1, num_patch]+1))[0]
    # Face s.t. zeta = 0
    if num_bound == 5:
        elems = np.where(
            igaPara._Nijk[2, patch_bounds[num_patch]:patch_bounds[num_patch+1]]
            == igaPara._Jpqr[2, num_patch]+1)[0]
    # Face s.t. zeta = 1
    if num_bound == 6:
        elems = np.where(
            igaPara._Nijk[2, patch_bounds[num_patch]:patch_bounds[num_patch+1]]
            == igaPara._Nkv[2, num_patch] - (igaPara._Jpqr[2, num_patch]+1))[0]
    return elems


def get_directionCP(igaPara, num_bound, num_patch=0, offset=0):
    """Get indices from boundary face control points of a given patch.

    Returns local indices of control points controlling the boundary indexed at
    `num_bound` of the patch indexed at `num_patch`. The `offset` argument
    enables to return the indices of control points forming a following row.


    Parameters
    ----------
    igaPara : IGA_parametrization
        NURBS parametrisation for IGA.
    num_bound : int
        Index of the boundary face.
        Must be ``1``, ``2``, ``3``, ``4``, ``5`` or ``6``.
    num_patch : int, optional
        Index of the patch. The default is 0.
    offset : int, optional
        Offset to get the indices of a following row of control points.
        The default is 0.

    Returns
    -------
    listCP : list of int
        List of local control points indices on the boundary.

    """
    indCP = get_boundCPindice(igaPara._Nkv, igaPara._Jpqr, num_bound,
                              num_patch, offset)
    listCP = igaPara._indCPbyPatch[num_patch][indCP]

    return listCP


def find_boundNumber(listCP,num_patch,Nkv,Jpqr,indCPbyPatch):
    '''
    Trouve la frontiere du patch #num_patch# contenu dans la liste des points de controle #listCP#

    Entree :
     * parametrisation Nurbs pour l'IGA
     * liste de points de controle
     * numero du patch
    Sortie :
     * retourne le numero de la frontiere concerne
    '''
    indices = np.where(np.isin(indCPbyPatch[num_patch], listCP))[0]
    bound = None
    if np.array_equal(indices, get_boundCPindice(Nkv, Jpqr, 3, num_patch)):
        bound = 3
    elif np.array_equal(indices, get_boundCPindice(Nkv, Jpqr, 4, num_patch)):
        bound = 4
    elif np.array_equal(indices, get_boundCPindice(Nkv, Jpqr, 1, num_patch)):
        bound = 1
    elif np.array_equal(indices, get_boundCPindice(Nkv, Jpqr, 2, num_patch)):
        bound = 2
    elif np.array_equal(indices, get_boundCPindice(Nkv, Jpqr, 5, num_patch)):
        bound = 5
    elif np.array_equal(indices, get_boundCPindice(Nkv, Jpqr, 6, num_patch)):
        bound = 6
    # else:
    #     print ' Warning: no bound has been found'
    return bound


def find_allBounds(listCP, Nkv, Jpqr, dim, indCPbyPatch, num_patch=0):
    '''
    Trouve les frontieres (faces, arretes, noeuds) du patch #num_patch# contenu dans la liste des
    points de controle #listCP#

    Entree :
     * parametrisation Nurbs pour l'IGA
     * liste de points de controle
     * numero du patch
    Sortie :
     * retourne les numeros des frontieres concernees
    '''
    indices = np.where(np.isin(indCPbyPatch[num_patch], listCP))[0]
    bound = []

    # face
    edgeIncluded = np.array([], dtype=np.intp)
    nodeIncluded = np.array([], dtype=np.intp)
    if dim[num_patch] == 3:
        for face in np.arange(0, 6) + 1:
            if np.all(np.isin(get_boundCPindice(Nkv, Jpqr, face,
                                                num_patch=num_patch), indices)
                      ):
                bound.append(face)
                if face == 1:
                    edgeIncluded = np.append(edgeIncluded,
                                             np.array([5, 7, 9, 11]))
                    nodeIncluded = np.append(nodeIncluded,
                                             np.array([1, 3, 5, 7]))
                elif face == 2:
                    edgeIncluded = np.append(edgeIncluded,
                                             np.array([6, 8, 10, 12]))
                    nodeIncluded = np.append(nodeIncluded,
                                             np.array([2, 4, 6, 8]))
                elif face == 3:
                    edgeIncluded = np.append(edgeIncluded,
                                             np.array([1, 3, 9, 10]))
                    nodeIncluded = np.append(nodeIncluded,
                                             np.array([1, 2, 5, 6]))
                elif face == 4:
                    edgeIncluded = np.append(edgeIncluded,
                                             np.array([2, 4, 11, 12]))
                    nodeIncluded = np.append(nodeIncluded,
                                             np.array([3, 4, 7, 8]))
                elif face == 5:
                    edgeIncluded = np.append(edgeIncluded,
                                             np.array([1, 2, 5, 6]))
                    nodeIncluded = np.append(nodeIncluded,
                                             np.array([1, 2, 3, 4]))
                elif face == 6:
                    edgeIncluded = np.append(edgeIncluded,
                                             np.array([3, 4, 7, 8]))
                    nodeIncluded = np.append(nodeIncluded,
                                             np.array([5, 6, 7, 8]))
    elif dim[num_patch] == 2:
        if np.all(np.isin(get_boundCPindice(Nkv, Jpqr, 6, num_patch=num_patch),
                          indices)):
            bound.append(6)
            edgeIncluded = np.append(edgeIncluded, np.array([1, 2, 3, 4]))
            nodeIncluded = np.append(nodeIncluded, np.array([1, 2, 3, 4]))

    # edge
    edges = get_edgeCPindice(Nkv, Jpqr, dim, num_patch=num_patch)
    if dim[num_patch] == 3:
        idx = np.setxor1d(np.arange(1, 12 + 1), edgeIncluded)
    elif dim[num_patch] == 2:
        idx = np.setxor1d(np.arange(1, 4 + 1), edgeIncluded)
    else:
        idx = np.array([1], dtype=np.intp)
    vertex2rm = np.array([], dtype=np.intp)
    for e in idx:
        if np.all(np.isin(edges[e - 1], indices)):
            bound.append(e + 100)
            vertex2rm = np.append(vertex2rm,
                                  np.array([edges[e - 1][0], edges[e - 1][-1]])
                                  )

    # vertex
    vertex = get_vertexCPindice(Nkv, Jpqr, dim, num_patch=num_patch)
    nodeIncluded = np.append(nodeIncluded,
                             np.where(np.isin(vertex, vertex2rm))[0] + 1)

    idx = np.setxor1d(np.arange(1, 2**dim[num_patch] + 1), nodeIncluded)
    for v in idx:
        if vertex[v - 1] in indices:
            bound.append(v + 1000)

    return bound


def get_edgeCPindice(Nkv, Jpqr, dim, num_patch=0):
    '''
    Retourne les indices des points de controle decrivant les arretes du patch #num_patch#
    '''
    nb_cpD = np.intp(np.maximum(Nkv[:, num_patch] - (Jpqr[:, num_patch] + 1),
                                np.ones(3)))

    edgeXi = np.arange(0, nb_cpD[0], 1, dtype=np.intp)
    edgeEta = np.arange(0, nb_cpD[1], 1) * np.prod(nb_cpD[0:1], dtype=np.intp)
    edgeZeta = np.arange(0, nb_cpD[2], 1) * np.prod(nb_cpD[0:2], dtype=np.intp)

    vertex = get_vertexCPindice(Nkv, Jpqr, dim, num_patch)

    edges = None
    if dim[num_patch] == 1:
        edges = [edgeXi]
    elif dim[num_patch] == 2:
        edges = [
            edgeXi,
            edgeXi + vertex[2],
            edgeEta,
            edgeEta + vertex[1]]
    elif dim[num_patch] == 3:
        edges = [
            edgeXi,
            edgeXi + vertex[2],
            edgeXi + vertex[4],
            edgeXi + vertex[6],
            edgeEta,
            edgeEta + vertex[1],
            edgeEta + vertex[4],
            edgeEta + vertex[5],
            edgeZeta,
            edgeZeta + vertex[1],
            edgeZeta + vertex[2],
            edgeZeta + vertex[3]]
    return edges


def get_vertexCPindice(Nkv, Jpqr, dim, num_patch=0):
    """
    Return indices of control points describing the corners of a given patch

    Parameters
    ----------
    Nkv : 2D array of float
        knot vectors size for all patches
    Jpqr : 2D array of int
        degrees for all patches
    dim : array of int
        dimension for all patches
    num_patch : int, optional
        Index of the patch. The default is 0

    Returns
    -------
    vertex : array
        corner points of the patch
    """

    nb_cpD = np.maximum(Nkv[:, num_patch] - (Jpqr[:, num_patch] + 1), 1)
    nb_cpD01 = nb_cpD[0] * nb_cpD[1]

    vertex_all = np.array([0,
                           nb_cpD[0] - 1,
                           (nb_cpD[1] - 1)*nb_cpD[0],
                           nb_cpD[1]*nb_cpD[0] - 1,
                           (nb_cpD[2] - 1)*nb_cpD01,
                           (nb_cpD[2] - 1)*nb_cpD01 + nb_cpD[0] - 1,
                           np.prod(nb_cpD) - nb_cpD[0],
                           np.prod(nb_cpD) - 1],
                          dtype=np.intp)

    if dim[num_patch] == 1:
        vertex = vertex_all[0:2]
    elif dim[num_patch] == 2:
        vertex = vertex_all[0:4]
    else:
        vertex = vertex_all
    return vertex


def add_displacementBC(igaPara, listCP, direction, value):
    """Add a displacement boundary condition.

    Adds displacement boundary condition for the control points listed in
    `listCP` in the direction given in `direction` with magnitude of `value`.


    Parameters
    ----------
    igaPara : IGA_parametrization
        NURBS parametrisation for IGA.
    listCP : list
        List of control points indices to apply boundary condition.
    direction : int
        Direction of the boundary condition. Must be ``1``, ``2`` or ``3``.
    value : float
        Boundary condition magnitude.

    Returns
    -------
    None. The IGA_parametrization object is updated with given information.

    """
    # num_bc = igaPara._nb_bc
    igaPara._nb_bc += 1

    igaPara._bc_values = \
        np.concatenate((igaPara._bc_values, np.vstack([direction, value])),
                       axis=1)

    nb_cp_thisBC = np.size(listCP)
    igaPara._bc_target.append(listCP)
    igaPara._bc_target_nbelem = \
        np.append(igaPara._bc_target_nbelem, nb_cp_thisBC)

    igaPara._update_dof_info()

    print(' Add displacement of {} in direction {} to {} CPs.'.format(
        value, direction, nb_cp_thisBC))

    return None

def find_orientation(test_cps:np.ndarray, Nkv:np.ndarray, Jpqr:np.ndarray, dim:np.ndarray,
                     indCPbyPatch:list, num_bound:int, num_patch=0, offset=0):
    """Find the orientation that orders the control points located at the bound ``num_bound``
    of patch ``num_patch`` has given by ``test_cps``.

    Parameters
    ----------
    Nkv : list of int
        Number of knots by direction.
    Jpqr : list of int
        Degree by direction.
    dim : array of int
        dimension for all patches
    indCPbyPatch : list of array of int
        Control point indices for all patches
    num_bound : int
        Index of the boundary entities (faces, edges, or vertices).
    num_patch : int, optional
        Index of the patch (python numbering). The default is 0.
    offset : int, optional
        Offset to get the indices of a following row of control points.
        The default is 0.

    Returns
    -------
    num_orientation : int
        Integer between 0 and 7.
    """
    icps = get_boundCPindice_wEdges(Nkv, Jpqr, dim, num_bound, num_patch, offset)
    cps = indCPbyPatch[num_patch][icps]
    if not np.all(np.isin(cps, test_cps)):
        raise ValueError("``test_cps`` are not defining the bound \
                         {num_bound} for the patch {num_patch}")
    if np.all(cps==test_cps):
        return 0
    argsort1 = np.argsort(test_cps)
    argsort2 = np.argsort(cps)
    icps_ordered = icps[argsort2][np.argsort(argsort1)]
    for num_orientation in np.arange(8)[1:]:
        test_icps = get_boundCPindice_wEdges(Nkv, Jpqr, dim, num_bound, num_patch,
                                             offset, num_orientation=num_orientation)
        if np.all(test_icps==icps_ordered):
            return num_orientation
    raise ValueError("No appropriate orientation has been found.")

def get_interface(num_patch:int, Nkv:np.ndarray,Jpqr:np.ndarray, dim:np.ndarray, indCPbyPatch:list):
    '''Identifies the strong coupling interfaces between ``num_patch`` and all other patches.

    Parameters
    ----------
    num_patch : int
        Index of the patch (python numbering).
    Nkv : list of int
        Number of knots by direction.
    Jpqr : list of int
        Degree by direction.
    dim : array of int
        dimension for all patches
    indCPbyPatch : list of array of int
        Control point indices for all patches

    Return
    ------
    tab : array of int
        A table with 5 colums and as many rows as found interfaces with other patches.
        Columns are of type:
        [index master patch, master bound, index slave patch, slave bound, orientation].
    '''
    tab_thispatch = []
    for i in range(0, num_patch):
        cps = np.intersect1d(indCPbyPatch[i], indCPbyPatch[num_patch])
        if np.size(cps):
            boundi = find_allBounds(cps, Nkv, Jpqr, dim, indCPbyPatch, i)
            bounds = find_allBounds(cps, Nkv, Jpqr, dim, indCPbyPatch,
                                    num_patch)

            n = np.argmax(np.array([np.size(boundi), np.size(bounds)]))
            if n == 0:
                for bnd in boundi:
                    icps_bound = get_boundCPindice_wEdges(Nkv, Jpqr, dim, bnd,
                                                          i)
                    cps = indCPbyPatch[i][icps_bound]
                    bounds = find_allBounds(cps, Nkv, Jpqr, dim, indCPbyPatch,
                                            num_patch)

                    orientation = find_orientation(
                        cps, Nkv, Jpqr, dim, indCPbyPatch, bounds[0], num_patch)

                    tab_thispatch.append([i, bnd, num_patch, bounds[0], orientation])
            else:
                for bnd in bounds:
                    icps_bound = get_boundCPindice_wEdges(Nkv, Jpqr, dim, bnd,
                                                          num_patch)
                    cps = indCPbyPatch[num_patch][icps_bound]
                    boundi = find_allBounds(cps, Nkv, Jpqr, dim,
                                            indCPbyPatch, i)

                    orientation = find_orientation(
                        cps, Nkv, Jpqr, dim, indCPbyPatch, boundi[0], i)

                    tab_thispatch.append([i, boundi[0], num_patch,bnd,orientation])
            # j = 0
            # for bnd in bounds[:n]:
            #     tab_thispatch.append([i,boundi[j],num_patch,bnd])
            #     j += 1

    return np.array(tab_thispatch, dtype=np.intp)

def get_patchConnectionInfos(Nkv, Jpqr, dim, indCPbyPatch, nb_patch):
    '''
    Boucle sur les patch pour trouver les connections entre chacun d'eux. Ces informations sont
    retournees sous forme d'un tableau. Il est utilise dans les procedures de raffinement pour
    raccorder les patch une fois le raffinement effectue.
    '''
    tab = np.array([], dtype=np.intp)
    for i in range(1, nb_patch):
        tab_thispatch = get_interface(i, Nkv, Jpqr, dim, indCPbyPatch)
        if np.size(tab_thispatch) > 0:
            try:
                tab = np.concatenate((tab, tab_thispatch), axis=0)
            except:
                tab = tab_thispatch
    return tab.astype(np.intp)


def build_rotRGmodesBYinterpolation(igapara,activePatch=None):
    pts = getgrevabscphysicalcoords( *igapara.get_inputs4grevphyscoords(activePatch=activePatch) )
    Idata,Iindices,Iindptr = buildgrevinterpolmat(
        *igapara.get_inputs4interpolation(activePatch=activePatch) )
    IMATRX = sp.csr_matrix((Idata,Iindices,Iindptr))

    mcrd  = igapara._mcrd
    nbpts = pts.shape[1]
    if mcrd == 2:
        RGmodes = np.zeros((nbpts*mcrd,1),dtype=np.float64)
    else:
        RGmodes = np.zeros((nbpts*mcrd,mcrd),dtype=np.float64)

    nbcpD = np.maximum(igapara._Nkv - (igapara._Jpqr+1),1)
    nbcp  = np.prod(nbcpD,axis=0)
    counts= 0; countp=0
    Urot  = []
    for patch in np.where(activePatch==1)[0]:

        pts2interpol = pts[:,countp:countp+nbcp[patch]]
        countp += nbcp[patch]

        Rrotz = np.cross(np.array([0,0,1]),pts2interpol.T)[:,:2]
        if mcrd==3:
            Rrotx  = np.cross(np.array([1,0,0]),pts2interpol.T)[:,(1,2)]
            Rroty  = np.cross(np.array([0,1,0]),pts2interpol.T)[:,(0,2)]
            Rrot   = np.concatenate((Rrotx,Rroty,Rrotz),axis=1)
        else:
            Rrot = np.vstack(Rrotz.flatten())


        IMAT0 = IMATRX[counts:counts+nbcpD[0,patch],:nbcpD[0,patch]].tocsc()
        LU0   = sp.linalg.splu(IMAT0)
        counts += nbcpD[0,patch]
        if igapara._dim[patch]>1:
            IMAT1 = IMATRX[counts:counts+nbcpD[1,patch],:nbcpD[1,patch]].tocsc()
            LU1   = sp.linalg.splu(IMAT1)
            counts += nbcpD[1,patch]
        if igapara._dim[patch]>2:
            IMAT2 = IMATRX[counts:counts+nbcpD[2,patch],:nbcpD[2,patch]].tocsc()
            LU2   = sp.linalg.splu(IMAT2)
            counts += nbcpD[2,patch]

        if igapara._dim[patch]==1:
            # curve
            Upatch = LU0.solve([Rrot])
        if igapara._dim[patch]==2:
            # surface
            indx = np.arange(0,nbcpD[0,patch])
            indy = np.arange(0,nbcpD[1,patch])
            Ustep1 = np.zeros_like(Rrot)
            Upatch = np.zeros_like(Rrot)
            for j in indy:
                ind = indx + j*indx.size
                Ustep1[ind,:] = LU0.solve(Rrot[ind,:])
            for i in indx:
                ind = indy*indx.size + i
                Upatch[ind,:] = LU1.solve(Ustep1[ind,:])

        if mcrd == 3:
            Upatch = np.insert(Upatch,[0,3,6],0.,axis=1)
            Urx = np.ravel(Upatch[:,0:3])
            Ury = np.ravel(Upatch[:,3:6])
            Urz = np.ravel(Upatch[:,6:9])
            Upatch = np.array([Urx,Ury,Urz]).T

        Urot.append([Upatch])

    return np.block(Urot)


def get_slaveCP4rotationBC(igaPara, bound, num_patch):
    listCP = igaPara._indCPbyPatch[num_patch][
        get_boundCPindice(igaPara, bound, num_patch, 1)]

    return listCP


def add_symmetryBC(igaPara, listCP, plan):
    '''
    Ajout d'une condition aux limites de type symmetry par rapport au plan #plan# pour les points
    de controle listes dans #listCP#.

    Entree :
     * parametrisation Nurbs pour l'IGA
     * liste de points de controle
     * plan de symmetry (1:Oxy, 2:Oxz, 3:Oyz)
    Sortie :
     * parametrisation Nurbs avec une condition de symmetry
    '''

    print('\nAdd symmetry condition.')

    listSlaveCP = np.array([])
    for num_patch in np.arange(0, igaPara._nb_patch):
        bound = find_boundNumber(igaPara, listCP, num_patch)
        if bound:
            listSlaveCP = np.union1d(
                listSlaveCP, get_slaveCP4rotationBC(igaPara, bound, num_patch))

    num_bc = igaPara._nb_bc
    igaPara._nb_bc += 3

    nb_cp_thisBC = np.size(listCP)
    igaPara._bc_target[num_bc,0:nb_cp_thisBC] = listCP[:]
    igaPara._bc_target_nbelem[num_bc] = nb_cp_thisBC
    if plan == 1:   # sym Oxy
        igaPara._bc_values[num_bc,:]   = np.array([3., 0.])
        igaPara._bc_values[num_bc+1,:] = np.array([7., 0.])
        igaPara._bc_values[num_bc+2,:] = np.array([8., 0.])
    elif plan == 2: # sym Oxz
        igaPara._bc_values[num_bc,:]   = np.array([2., 0.])
        igaPara._bc_values[num_bc+1,:] = np.array([7., 0.])
        igaPara._bc_values[num_bc+2,:] = np.array([9., 0.])
    elif plan == 3: # sym Oyz
        igaPara._bc_values[num_bc,:]   = np.array([1., 0.])
        igaPara._bc_values[num_bc+1,:] = np.array([8., 0.])
        igaPara._bc_values[num_bc+2,:] = np.array([9., 0.])

    igaPara._bc_target[num_bc+1, 0:nb_cp_thisBC] = listSlaveCP[:]
    igaPara._bc_target_nbelem[num_bc+1] = nb_cp_thisBC
    igaPara._bc_target[num_bc+2, 0:nb_cp_thisBC] = listSlaveCP[:]
    igaPara._bc_target_nbelem[num_bc+2] = nb_cp_thisBC

    igaPara._bc_target[num_bc+1 + 200, 0:nb_cp_thisBC] = listCP[:]
    igaPara._bc_target_nbelem[num_bc+1 + 200] = nb_cp_thisBC
    igaPara._bc_target[num_bc+2 + 200, 0:nb_cp_thisBC] = listCP[:]
    igaPara._bc_target_nbelem[num_bc+2 + 200] = nb_cp_thisBC

    return None


def clean_MasterSlave_BC(igaPara):
    '''
    Nettoyage des conditions aux limites type condition de symmetrie. Test sur les ddl esclaves et
    les ddl maitres pour trouver ceux qui sont a la fois maitre et esclave. On enleve leur role de
    ddl maitre par transitivite. Si ddl A est maitre de B et esclave de C, alors C devient maitre
    de B (et de A).

    Entree :
     * parametrisation Nurbs pour l'IGA
    Sortie :
     * parametrisation nettoyee : suppression des BCs incompatibles
    '''

    print('\nClean BCs')

    numBC_MtrSlv = np.where( igaPara._bc_values[:,0] > 6 )[0]

    M_cp_slave  = igaPara._bc_target[numBC_MtrSlv,:]
    M_cp_master = igaPara._bc_target[numBC_MtrSlv+200,:]

    M_ddl_slave  = (M_cp_slave-1)*igaPara._mcrd \
                   + np.transpose(np.kron(np.ones((1000,1)), igaPara._bc_values[numBC_MtrSlv,0]-6))
    M_ddl_slave  = np.maximum(M_ddl_slave, 0)
    M_ddl_master = (M_cp_master-1)*igaPara._mcrd \
                   + np.transpose(np.kron(np.ones((1000,1)), igaPara._bc_values[numBC_MtrSlv,0]-6))
    M_ddl_master = np.maximum(M_ddl_master, 0)

    ddl_slave = np.unique(
        np.reshape(M_ddl_slave,  np.size(M_ddl_slave ))[ np.flatnonzero(M_ddl_slave ) ])
    ddl_master= np.unique(
        np.reshape(M_ddl_master, np.size(M_ddl_master))[ np.flatnonzero(M_ddl_master) ])


    ddl_masterANDslave = np.intersect1d(ddl_slave, ddl_master)
    for ddl in ddl_masterANDslave:
        indSlave_needNewMaster = np.where( M_ddl_master == ddl)
        indNewMaster = np.where(M_ddl_slave == ddl)

        print('   Replace', \
            igaPara._bc_target[numBC_MtrSlv[indSlave_needNewMaster[0]]+200,
                               indSlave_needNewMaster[1] ], \
            'with', igaPara._bc_target[numBC_MtrSlv[indNewMaster[0]]+200, indNewMaster[1] ])

        igaPara._bc_target[numBC_MtrSlv[indSlave_needNewMaster[0]]+200,
                           indSlave_needNewMaster[1] ] \
            = igaPara._bc_target[numBC_MtrSlv[indNewMaster[0]]+200, indNewMaster[1] ]

    return None


def add_bendingStrip(igaPara, patch1, patch2):
    '''
    Ajout d'un bending strip pour le couplage de patch de type shell Kirchhoff-Love.
    '''

    # points de controle sur la frontiere entre les patchs
    cp_bound = np.intersect1d(igaPara._indCPbyPatch[patch1], igaPara._indCPbyPatch[patch2])

    if not cp_bound.size:
        print('Patch %i and %i are not neighbours. No bending strip added.' % (patch1,patch2))
    else:
        # numero des frontieres
        num_bound1 = find_boundNumber(igaPara, cp_bound, patch1)
        num_bound2 = find_boundNumber(igaPara, cp_bound, patch2)

        # numero des cp definissant l'element de flexion
        ind_cp_patch0 = get_boundCPindice(igaPara,num_bound1,patch1,0)
        ind_cp_patch1 = get_boundCPindice(igaPara,num_bound1,patch1,1)
        ind_cp_patch2 = get_boundCPindice(igaPara,num_bound2,patch2,1)

        # points de controle definissant l'element de flexion
        cp_patch0 = igaPara._indCPbyPatch[patch1][ind_cp_patch0]
        cp_patch1 = igaPara._indCPbyPatch[patch1][ind_cp_patch1]
        cp_patch2 = igaPara._indCPbyPatch[patch2][ind_cp_patch2]

        # ajout definition de l'element de flexion
        numbc = igaPara._nb_bc
        igaPara._bc_values[numbc,0] = 11
        n = cp_bound.size
        igaPara._bc_target_nbelem[numbc] = n
        igaPara._bc_target[numbc,0:n] = cp_patch0[:]
        igaPara._bc_target[numbc+500,0:n] = cp_patch1[:]
        igaPara._bc_target[numbc+600,0:n] = cp_patch2[:]
        igaPara._nb_bc += 1

        print('\nBending strip added between patch %i and %i.' % (patch1,patch2))

    return None


def add_strongCoupling(igaPara, patch1, patch2):
    '''
    Ajout d'un couplage fort entre le patch #patch1# et le #patch2#. Couplage possible pour des
    maillages et geometries compatibles, et des raccords tangents ou droits uniquement !

    Entree :
     * parametrisation Nurbs pour l'IGA
     * patch 1
     * patch 2
    Sortie :
     * parametrisation Nurbs avec un couplage fort
    '''

    # points de controle sur la frontiere entre les patchs
    cp_bound = np.intersect1d(igaPara._indCPbyPatch[patch1], igaPara._indCPbyPatch[patch2])

    if not cp_bound.size:
        print('Patch %i and %i are not neighbours. No coupling added.' % (patch1,patch2))
    else:
        # numero des frontieres
        num_bound1 = find_boundNumber(igaPara, cp_bound, patch1)
        num_bound2 = find_boundNumber(igaPara, cp_bound, patch2)

        # numero des cp definissant l'element de flexion
        ind_cp_patch1 = get_boundCPindice(igaPara,num_bound1,patch1,1)
        ind_cp_patch2 = get_boundCPindice(igaPara,num_bound2,patch2,1)

        # points de controle definissant l'element de flexion
        cp_patch1 = igaPara._indCPbyPatch[patch1][ind_cp_patch1]
        cp_patch2 = igaPara._indCPbyPatch[patch2][ind_cp_patch2]

        # ajout definition de l'element de flexion
        numbc = igaPara._nb_bc
        igaPara._bc_values[numbc,0] = 10
        n = cp_bound.size
        igaPara._bc_target_nbelem[numbc] = n
        igaPara._bc_target[numbc,0:n] = cp_patch2[:] # ----- slave
        igaPara._bc_target[numbc+300,0:n] = cp_bound[:] # -- master 1
        igaPara._bc_target[numbc+400,0:n] = cp_patch1[:] # - master 2
        igaPara._nb_bc += 1

    print('\nStrong coupling added between patch %i and %i.' % (patch1,patch2))

    return None


def add_slaveCPs(igaPara, slaveCPs, masterCP):
    '''

    '''
    # Modification table des elements
    for slaveCP in slaveCPs:
        ind4IEN = np.where(igaPara._IEN == slaveCP)
        igaPara._IEN[ind4IEN] = masterCP

    # Propagation des CL des cps esclaves au cp maitre
    for slaveCP in slaveCPs:
        ind4BCs = np.where(igaPara._bc_target == slaveCP)
        igaPara._bc_target[ind4BCs] = masterCP

    # Ajout des CL pour supprimer les d.o.f. des points de controle escalves
    for direction in list([1,2,3]):
        add_displacementBC(igaPara, slaveCPs, direction, 0.)


    return None


def get_tabCLG(igaPara):
    num_lgrge= np.where(igaPara._ELT_TYPE == 'U4' )[0] + 1
    nb_lgrge = np.size(num_lgrge)

    tabCPLG      = np.zeros((nb_lgrge,5),dtype=np.intp)
    tabCPLG[:,0] = num_lgrge[:]

    num_crve = np.where(igaPara._ELT_TYPE == 'U00')[0]
    for crv in num_crve:
        patch,lgrge,ismaster = np.intp(igaPara._PROPS[crv][1:])
        i = np.where(num_lgrge == lgrge)
        tabCPLG[i,2-ismaster] = patch
        tabCPLG[i,4-ismaster] = crv+1
    return tabCPLG


def cplginterface_infos(igaPara,GPs):
    num_patch  = -1
    face_found = False
    while num_patch < igaPara._nb_patch-1 and not face_found:
        num_patch += 1
        vertex = get_vertexCPindice(igaPara._Nkv,igaPara._Jpqr,igaPara._dim,num_patch)
        face = 0
        while face < 6 and not face_found:
            face += 1

            zeta = 0.
            if face == 1:
                A = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[0]]-1]
                B = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[2]]-1]
                C = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[4]]-1]
                direction = list([1,2,0])
                #if not igaPara._ELT_TYPE[num_patch] == 'U3':
                zeta = np.min(igaPara._Ukv[num_patch][0])
            if face == 2:
                A = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[1]]-1]
                B = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[3]]-1]
                C = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[5]]-1]
                direction = list([1,2,0])
                #if not igaPara._ELT_TYPE[num_patch] == 'U3':
                zeta = np.max(igaPara._Ukv[num_patch][0])
            if face == 3:
                A = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[0]]-1]
                B = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[1]]-1]
                C = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[4]]-1]
                direction = list([0,2,1])
                #if not igaPara._ELT_TYPE[num_patch] == 'U3':
                zeta = np.min(igaPara._Ukv[num_patch][1])
            if face == 4:
                A = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[2]]-1]
                B = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[4]]-1]
                C = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[6]]-1]
                direction = list([0,2,1])
                #if not igaPara._ELT_TYPE[num_patch] == 'U3':
                zeta = np.max(igaPara._Ukv[num_patch][1])
            if face == 5:
                A = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[0]]-1]
                B = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[1]]-1]
                C = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[2]]-1]
                direction = list([0,1,2])
                if not igaPara._ELT_TYPE[num_patch] == 'U3':
                    zeta = np.min(igaPara._Ukv[num_patch][2])
            if face == 6:
                A = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[4]]-1]
                B = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[5]]-1]
                C = igaPara._COORDS[:,igaPara._indCPbyPatch[num_patch][vertex[6]]-1]
                direction = list([0,1,2])
                if not igaPara._ELT_TYPE[num_patch] == 'U3':
                    zeta = np.max(igaPara._Ukv[num_patch][2])

            zeta *= np.ones(np.size(GPs,1))

            if face < 5 and igaPara._ELT_TYPE[num_patch] == 'U3' or igaPara._mcrd == 2:
                AG = GPs - np.repeat(np.vstack(A),np.size(GPs,1),axis=1)
                AB = np.repeat(np.vstack(B-A),np.size(GPs,1),axis=1)
                if np.all(np.isclose(
                        np.linalg.norm(np.cross(AB.transpose(),AG.transpose()),axis=1),0.)):
                    face_found = True

                    xi_m  = np.min(igaPara._Ukv[num_patch][direction[0]])
                    xi_M  = np.max(igaPara._Ukv[num_patch][direction[0]])
                    Lx    = np.linalg.norm(B-A)
                    xi    = (xi_M-xi_m)/Lx * np.linalg.norm(AG,axis=0) + xi_m
                    eta   = np.zeros_like(xi)

            else:
                n = np.cross(B-A,C-A)
                n/= np.linalg.norm(n)

                AG = GPs - np.repeat(np.vstack(A),np.size(GPs,1),axis=1)
                if np.all(np.isclose(np.dot(n,AG),0.)):
                    face_found = True

                    xi_m  = np.min(igaPara._Ukv[num_patch][direction[0]])
                    xi_M  = np.max(igaPara._Ukv[num_patch][direction[0]])
                    eta_m = np.min(igaPara._Ukv[num_patch][direction[1]])
                    eta_M = np.max(igaPara._Ukv[num_patch][direction[1]])
                    Lx    = np.linalg.norm(B-A)
                    Ly    = np.linalg.norm(C-A)
                    Vi    = (B-A)/Lx
                    Vjb   = np.cross(n,Vi)
                    Vj    = (C-A)/Ly
                    Vib   = np.cross(n,Vj)

                    #xi    = ( xi_M -  xi_m)/Lx * np.dot((B-A)/Lx,AG) +  xi_m
                    #eta   = (eta_M - eta_m)/Ly * np.dot((C-A)/Ly,AG) + eta_m

                    xi    = ( xi_M -  xi_m)/Lx * np.dot(Vib,AG)/np.dot(Vi,Vib) +  xi_m
                    eta   = (eta_M - eta_m)/Ly * np.dot(Vjb,AG)/np.dot(Vj,Vjb) + eta_m

    GaussCoords = np.zeros_like(GPs)
    GaussCoords[direction,:] = np.array([xi,eta,zeta])

    u0 = igaPara._Ukv[num_patch][0][
        igaPara._Jpqr[0,num_patch]+1:igaPara._Nkv[0,num_patch]-igaPara._Jpqr[0,num_patch]-1]
    u1 = igaPara._Ukv[num_patch][1][
        igaPara._Jpqr[1,num_patch]+1:igaPara._Nkv[1,num_patch]-igaPara._Jpqr[1,num_patch]-1]
    u2 = np.array([])
    if igaPara._dim[num_patch]>2:
        u2 = igaPara._Ukv[num_patch][2][
            igaPara._Jpqr[2,num_patch]+1:igaPara._Nkv[2,num_patch]-igaPara._Jpqr[2,num_patch]-1]
    ukv = np.array([u0,u1,u2])

    test4nijk0 = np.searchsorted(ukv[0], GaussCoords[0,:]) + igaPara._Jpqr[0,num_patch] + 1
    test4nijk1 = np.searchsorted(ukv[1], GaussCoords[1,:]) + igaPara._Jpqr[1,num_patch] + 1
    test4nijk2 = np.zeros_like(test4nijk0)
    if igaPara._dim[num_patch]>2:
        test4nijk2 = np.searchsorted(ukv[2], GaussCoords[2,:]) + igaPara._Jpqr[2,num_patch] + 1
    test4nijk = np.array([test4nijk0,test4nijk1,test4nijk2], dtype=np.intp)

    offset = np.insert(np.cumsum(igaPara._elementsByPatch),0,0)[num_patch] + 1
    nb_gps = np.size(GPs,1)
    IGJelem = np.zeros(nb_gps,dtype=np.intp)
    i = 0
    for test in test4nijk.transpose():
        IGJelem[i] = np.where(np.all(test == igaPara._Nijk.transpose(), axis=1))[0][0] + offset
        i += 1

    return GaussCoords,IGJelem


def bezier_decomposition_patch(igapara,numpatch=0,return_ien=False,
                               return_mat=False,return_invmat=False):
    '''
    Beizer decomposition of a NURBS patch.
    Inputs:
    * igapara --- an iga parametrization
    * numpatch -- patch number (optional, default=0)
    Outputs:
    * cps per element (coords + weights)
    * degree per element
    '''

    M = []
    ukv = igapara._Ukv[numpatch]
    jpqr= igapara._Jpqr[:,numpatch]
    ncpF= np.ones(3,dtype=int)
    nelF= np.ones(3,dtype=int)
    for d in range(len(ukv)):
        p = igapara._Jpqr[d,numpatch]
        u = ukv[d]
        n = u.size - (p+1)
        M.append(sp.csc_matrix(decomposition(u,n,p)))

        nelF[d] = np.unique(u).size-1
        ncpF[d] = nelF[d]*p + 1

    cpc = igapara._COORDS[:,igapara._indCPbyPatch[numpatch]-1].T
    cpw = np.vstack(igapara._vectWeight[igapara._indCPbyPatch[numpatch]-1])
    ncpD= np.maximum(igapara._Nkv[:,numpatch]-(igapara._Jpqr[:,numpatch]+1),1)

    M0 = sp.kron(sp.eye(ncpD[2]), sp.kron(sp.eye(ncpD[1]), M[0]) )
    cpc0 = M0 * np.block([cpc*cpw,cpw])
    M1 = sp.kron(sp.eye(ncpD[2]), sp.kron(M[1], sp.eye(ncpF[0])) )
    cpc1 = M1 * cpc0
    M2 = sp.kron(M[2], sp.kron(sp.eye(ncpF[1]), sp.eye(ncpF[0])) )
    cpc2 = M2 * cpc1
    cpc2[:,:-1] /= np.vstack(cpc2[:,-1])

    nnode=jpqr+1
    icpD = []
    for d in range(len(ukv)):
        icp = np.tile(np.arange(nnode[d]),nelF[d]) \
            + np.repeat(np.arange(nelF[d])*(nnode[d]-1),nnode[d])
        icpD.append(np.reshape(icp,(nelF[d],nnode[d])))

    nijk = np.block([[np.tile(np.arange(nelF[0]),nelF[1]*nelF[2])],
                     [np.tile(np.repeat(np.arange(nelF[1]),nelF[0]),nelF[2])],
                     [np.repeat(np.arange(nelF[2]),nelF[0]*nelF[1])]]).T

    ien = np.tile(np.arange(nnode[0]),(np.prod(nelF),np.prod(nnode[1:]))) \
        + np.tile(np.repeat(np.arange(nnode[1]),nnode[0]),(np.prod(nelF),nnode[2]))*ncpF[0]\
        + np.repeat(np.arange(nnode[2]),np.prod(nnode[:2]))*np.prod(ncpF[:2])

    offset = nijk[:,0]*(nnode[0]-1) \
        + nijk[:,1]*(ncpF[0]*(nnode[1]-1)) \
        + nijk[:,2]*ncpF[0]*ncpF[1]*(nnode[2]-1)
    ien += np.vstack(offset)

    cp_decomp  = cpc2[ien]
    jpqr_decomp= np.tile(jpqr,(offset.size,1))
    toreturn = [cp_decomp,jpqr_decomp]
    if return_ien:
        toreturn.append(ien)
    Mtot = M2*M1*M0
    if return_mat:
        toreturn.append(Mtot)
    if return_invmat:
        invM = []
        for Mi in M:
            MiTMi = sp.csc_matrix(Mi.T*Mi)
            invM.append(sp.linalg.inv(MiTMi)*Mi.T)
        invM0 = sp.kron(sp.eye(ncpD[2]), sp.kron(sp.eye(ncpD[1]), invM[0]) )
        invM1 = sp.kron(sp.eye(ncpD[2]), sp.kron(invM[1], sp.eye(ncpF[0])) )
        invM2 = sp.kron(invM[2], sp.kron(sp.eye(ncpF[1]), sp.eye(ncpF[0])) )
        invMtot = invM0*invM1*invM2
        toreturn.append(invMtot)
    return toreturn


def bezier_decomposition_patch_discontinuous(igapara,numpatch=0,return_ien=False,
                                             return_mat=False,return_invmat=False):
    '''
    Bezier decomposition of a NURBS patch.
    Inputs:
    * igapara --- an iga parametrization
    * numpatch -- patch number (optional, default=0)
    Outputs:
    * cps per element (coords + weights)
    * degree per element
    '''

    M = []
    dimPatch = igapara._dim[numpatch]
    ukv = igapara._Ukv[numpatch]
    jpqr= igapara._Jpqr[:,numpatch]
    ncpF= np.ones(3, dtype=int)
    nelF= np.ones(3, dtype=int)
    for d in range(dimPatch):
        p = igapara._Jpqr[d,numpatch]
        u = ukv[d]
        n = u.size - (p+1)
        M.append(decomposition_sparse(u,n,p+1))

        nelF[d] = np.unique(u).size-1
        ncpF[d] = nelF[d]*(p+1)

    if dimPatch == 2:
        M.append(sp.eye(1).tocsc())

    cpc = igapara._COORDS[:,igapara._indCPbyPatch[numpatch]-1].T
    cpw = np.vstack(igapara._vectWeight[igapara._indCPbyPatch[numpatch]-1])
    ncpD= np.maximum(igapara._Nkv[:,numpatch]-(igapara._Jpqr[:,numpatch]+1),1)

    M0 = sp.kron(sp.eye(ncpD[2]), sp.kron(sp.eye(ncpD[1]), M[0]) )
    cpc0 = M0 * np.block([cpc*cpw,cpw])
    M1 = sp.kron(sp.eye(ncpD[2]), sp.kron(M[1], sp.eye(ncpF[0])) )
    cpc1 = M1 * cpc0

    M2 = sp.kron(M[2], sp.kron(sp.eye(ncpF[1]), sp.eye(ncpF[0])) )
    cpc2 = M2 * cpc1
    cpc2[:,:-1] /= np.vstack(cpc2[:,-1])


    nnode=jpqr+1
    icpD = []
    for d in range(len(ukv)):
        icp = np.tile(np.arange(nnode[d]),nelF[d]) \
            + np.repeat(np.arange(nelF[d])*(nnode[d]-1),nnode[d])
        icpD.append(np.reshape(icp,(nelF[d],nnode[d])))

    nijk = np.block([[np.tile(np.arange(nelF[0]),nelF[1]*nelF[2])],
                     [np.tile(np.repeat(np.arange(nelF[1]),nelF[0]),nelF[2])],
                     [np.repeat(np.arange(nelF[2]),nelF[0]*nelF[1])]]).T

    ien = np.tile(np.arange(nnode[0]),(np.prod(nelF),np.prod(nnode[1:]))) \
        + np.tile(np.repeat(np.arange(nnode[1]),nnode[0]),(np.prod(nelF),nnode[2]))*ncpF[0]\
        + np.repeat(np.arange(nnode[2]),np.prod(nnode[:2]))*np.prod(ncpF[:2])

    offset = nijk[:,0]*nnode[0] \
        + nijk[:,1]*(ncpF[0]*nnode[1]) \
        + nijk[:,2]*ncpF[0]*ncpF[1]*nnode[2]
    ien += np.vstack(offset)

    cp_decomp  = cpc2[ien]
    jpqr_decomp= np.tile(jpqr,(offset.size,1))
    toreturn = [cp_decomp,jpqr_decomp]
    if return_ien==True:
        toreturn.append(ien)
    Mtot = M2 @ M1 @ M0
    if return_mat==True:
        toreturn.append(Mtot.tocsr())
    if return_invmat==True:
        invM = []
        for Mi in M:
            MiTMi = sp.csc_matrix(Mi.T*Mi)
            invM.append(sp.linalg.inv(MiTMi)*Mi.T)
        invM0 = sp.kron(sp.eye(ncpD[2]), sp.kron(sp.eye(ncpD[1]), invM[0]) )
        invM1 = sp.kron(sp.eye(ncpD[2]), sp.kron(invM[1], sp.eye(ncpF[0])) )
        invM2 = sp.kron(invM[2], sp.kron(sp.eye(ncpF[1]), sp.eye(ncpF[0])) )
        invMtot = invM0*invM1*invM2
        toreturn.append(invMtot.tocsr())
    return toreturn


def bezier_decomposition_elem(igapara,jelem):
    '''
    Beizer decomposition of a NURBS patch.
    Inputs:
    * igapara -- an iga parametrization
    * jelem ---- element number
    Outputs:
    * cps (coords + weights)
    * degree
    '''
    if (jelem<0 or jelem>=igapara._nb_elem):
        raise ValueError("jelem should be an integer value in [0,{igapara._nb_elem}[")

    testiel  = jelem-np.cumsum(np.block([0,igapara._elementsByPatch]))
    numpatch = np.where(testiel>=0)[0][0]
    numelem  = testiel[numpatch]
    cp_decomp,jpqr_decomp = bezier_decomposition_patch(igapara,numpatch)
    return cp_decomp[numelem],jpqr_decomp[numelem]


class localExtractionPatch:
    '''This class is made for the function `bezier_decomposition`. It provides
    basis algorithms for the bezier decomposition of spline patches.
    '''
    def __init__(self, Ukv:list[np.ndarray], Jpqr:np.ndarray, dim:int) -> None:
        self._dim = dim
        self._Jpqr = Jpqr.copy()
        self._Ukv = Ukv.copy()

        uniqueC1D = []
        self._alltouniqueC1D = []
        for idim in range(dim):
            uniqueCi, alltouniqueCi = localExtraction1Dunique(Ukv[idim], Jpqr[idim])
            uniqueC1D.append(uniqueCi)
            self._alltouniqueC1D.append(alltouniqueCi)

        self._uniqueC = {}
        if dim == 1:
            for i0 in range(len(uniqueC1D[0])):
                C0 = uniqueC1D[0][i0]
                key = '%i'%i0
                self._uniqueC[key] = np.asfortranarray(C0)
        elif dim == 2:
            for i1 in range(len(uniqueC1D[1])):
                C1 = uniqueC1D[1][i1]
                for i0 in range(len(uniqueC1D[0])):
                    C0 = uniqueC1D[0][i0]
                    key = '%i,%i'%(i0,i1)
                    self._uniqueC[key] = np.kron(C1,C0).T
        elif dim == 3:
            for i2 in range(len(uniqueC1D[2])):
                C2 = uniqueC1D[2][i2]
                for i1 in range(len(uniqueC1D[1])):
                    C1 = uniqueC1D[1][i1]
                    C2C1 = np.kron(C2,C1)
                    for i0 in range(len(uniqueC1D[0])):
                        C0 = uniqueC1D[0][i0]
                        key = '%i,%i,%i'%(i0,i1,i2)
                        self._uniqueC[key] = np.kron(C2C1,C0).T

    def getC(self, Nijk:np.ndarray[int]):
        key = '%i' % self._alltouniqueC1D[0][Nijk[0]]
        for i in range(1,self._dim):
            key += ',%i' % self._alltouniqueC1D[i][Nijk[i]]
        return self._uniqueC[key]

    def setSparseMats(self, ncp:int):
        self._uniqueCsparse = {}
        for key,value in self._uniqueC.items():
            data = np.ravel(value,order='C')
            ien = np.arange(value.shape[1])
            indices = np.tile(ien,ien.size)
            indptr  = np.arange(ien.size+1)*ien.size
            self._uniqueCsparse[key] = sp.csr_matrix((data,indices,indptr),shape=(ien.size,ncp))

    def getCsparse(self, Nijk:np.ndarray[int]):
        if not hasattr(self, '_uniqueCsparse'):
            raise AttributeError("First one needs to set the unique sparse operators by calling `setSparseMats`.")
        key = '%i' % self._alltouniqueC1D[0][Nijk[0]]
        for i in range(1,self._dim):
            key += ',%i' % self._alltouniqueC1D[i][Nijk[i]]
        return self._uniqueCsparse[key]

# WARNING Operator | needs Python 3.10+
# def bezier_decomposition(igapara, listelem:None|np.ndarray=None, return_mat=False):
from typing import Optional
def bezier_decomposition(igapara, listelem: Optional[np.ndarray] = None, return_mat=False):
    '''Bezier decomposition of elements of a IGA parameterization.

    Parameters
    ----------
    igapara : IGAparametrization
        A IGA model, can be multipatch.
    listelem : None or Array of floats, optional
        An array containing the indices of the elements for which a Bezier expression
        is required. If None (default value), all the elements are decomposed.
    return_mat : bool, optional.
        If True, all the sparse extraction operators are outputed.

    Returns
    -------
    bezierDecomp : dict
        Dictionary that provides the Bezier representation of the elements.
        Keys are the element indices. Values are a list of two arrays:
        - First array gives the Bezier control point coordinates and weights,
        - Second array gives the element degrees per direction.
    bezierDecomp_mat : dict
        Returned if `return_mat is True`.
        Keys are the element indices. Values are a sparse matrices.
    '''
    if listelem is None:
        listelemsort = np.arange(igapara._nb_elem)
    else:
        listelemsort = np.sort(listelem)

    CPspline = np.block([[igapara._COORDS],[igapara._vectWeight]])
    CPspline[:-1] *= CPspline[-1]
    CPspline = np.asfortranarray(CPspline.T)

    localExtractionByPatches = {}
    bezierDecomp_mat = {}
    bezierDecomp = {}
    elemoffset = np.cumsum(np.block([0,igapara._elementsByPatch]))
    ipatch = 0
    for iel in listelemsort:
        while iel >= elemoffset[ipatch+1]:
            ipatch += 1
        if ipatch not in localExtractionByPatches.keys():
            localExtractionByPatches[ipatch] = localExtractionPatch(
                igapara._Ukv[ipatch], igapara._Jpqr[:,ipatch], igapara._dim[ipatch])
            if return_mat:
                localExtractionByPatches[ipatch].setSparseMats(igapara._nb_cp)

        nijk = igapara._Nijk[:,iel]
        ien = igapara._IEN[ipatch][iel-elemoffset[ipatch]][::-1]-1


        Cel = localExtractionByPatches[ipatch].getC(nijk)
        CPsplineEl = CPspline[ien,:]
        wSplineEl = CPsplineEl[:,-1]
        CPbezier = Cel @ CPsplineEl
        wBezierInv = 1/CPbezier[:,-1]
        CPbezier[:,:-1] *= np.vstack(wBezierInv)
        bezierDecomp[iel] = [CPbezier,igapara._Jpqr[:,ipatch]]

        if return_mat:
            bezierDecomp_mat[iel] = localExtractionByPatches[ipatch].getCsparse(nijk).copy()
            bezierDecomp_mat[iel].indices[:] = np.tile(ien,ien.size)
            data = np.ravel(Cel * blas.dger(1.0, wBezierInv, wSplineEl), order='C')
            bezierDecomp_mat[iel].data[:] = data[:]

    if return_mat:
        return bezierDecomp, bezierDecomp_mat
    else:
        return bezierDecomp
