"""
Linear analysis of a clamped square in 2D
Computation is made for different degrees and refinement levels
"""


import time

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp

# Number of samples in each direction for solution
N_SAMPLE = 500
# Degree for referenc solution
DEG_REFERENCE = 7
# Number of refinements for reference solution
REF_REFERENCE = 7


def compute(degree, refinement):
    """
    Compute clamped square case for a given degree and refinement level

    Parameters
    ----------
    degree: int
        degree in both directions
    refinement: int
        refinement level in both directions
    """

    # Read data and create IGAparametrization object
    modeleIGA = IGAparametrization(filename='square')

    # Refine modele
    nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)

    nb_ref[:, 0] = np.array([refinement, refinement, 0])
    nb_deg[:, 0] = np.array([degree-1, degree-1, 0])

    # Initial refinement (none)
    modeleIGA.refine(nb_ref, nb_deg)

    # Matrix assembly
    ndof = modeleIGA._nb_dof_free
    idof = modeleIGA._ind_dof_free[:ndof]-1

    t1 = time.time()

    data, row, col, Fb = build_stiffmatrix(
                            *modeleIGA.get_inputs4system_elemStorage())

    Kside = sp.coo_matrix((data, (row, col)),
                        shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                        dtype='float64').tocsc()
    Ktot = Kside + Kside.transpose()
    K2solve = Ktot[idof, :][:, idof]
    del Kside, data, row, col
    print(('\n Time to build stiffness matrix : %.2f s' % (time.time() - t1)))

    # RESOLUTION : direct solver
    t2 = time.time()
    x  = sp.linalg.spsolve(K2solve, Fb[idof])
    # LU = sp.linalg.splu(K2solve)
    # x = LU.solve(Fb[idof])
    print(('\n Time for direct solving :    %.2f s\n\n' % (time.time() - t2)))

    # Solution reconstruction
    SOL, _ = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

    resmul = pp.evaldispmulti(*modeleIGA.get_inputs4evaldispmulti(
        SOL.transpose(), N_SAMPLE, N_SAMPLE, 1))

    return resmul[:, :, 0, :2]


reference = compute(DEG_REFERENCE, REF_REFERENCE)

with open("convergence.txt", "w", encoding='UTF-8') as file:
    for degree in range(1, 7):
        file.write(str(degree)+'\t')
        for ref in range(8):
            result = compute(degree, ref)
            delta = result - reference
            error = np.linalg.norm(delta[:, :, None, :]@delta[:, :, :, None])
            file.write(str(error)+'\t')
        file.write('\n')
