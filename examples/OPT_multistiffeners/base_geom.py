"""
TODO
"""

import numpy as np
import scipy.sparse as sp

from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip


iga_model = IGAparametrization(filename='base_geom')

# Refine model
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

nb_deg[:, 0] = np.array([2, 1, 0])
additional_knots = {'patches': np.array([0]),
                    '1': np.array([0.14]),
                    '2': np.array([0.1, 0.3]),
                    '3': np.array([])}

iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

# Add loading on a single element
iga_model._indDLoad = np.array([np.array([3]),np.array([3])])
iga_model._JDLType = np.array([61,63])
iga_model._ADLMAG = np.array([1.190, 7.143])
iga_model._load_target_nbelem = np.array([1, 1])
iga_model._nb_load = 2
iga_model._indDLoad_flat = np.array([], dtype=np.intp)
for load in iga_model._indDLoad:
    iga_model._indDLoad_flat = np.hstack((iga_model._indDLoad_flat, load))

print(iga_model.get_mechanicalSettings()[2])

# Refinement to get elements with approx the same size
additional_knots = {'patches': np.array([0]),
                    '1': np.array([0.28, 0.42, 0.7, 0.84]),
                    '2': np.array([0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
                    '3': np.array([])}

nb_deg[:, 0] = np.array([0, 0, 0])
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

# OTHER REFINEMENTS

# END OF OTHER REFINEMENTS

# Matrix assembly
ndof = iga_model._nb_dof_free
idof = iga_model._ind_dof_free[:ndof]-1



data, row, col, Fb = build_stiffmatrix(
                            *iga_model.get_inputs4system_elemStorage())
Kside = sp.coo_matrix((data, (row, col)),
                      shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside + Kside.transpose()
K2solve = Ktot[idof, :][:, idof]
Fbsolve = Fb[idof]
del Kside, data, row, col




# Resolution
x = sp.linalg.spsolve(K2solve, Fbsolve)

# Solution reconstruction
SOL, u = rsol.reconstruction(**iga_model.get_inputs4solution(x))
pp.generatevtu(*iga_model.get_inputs4postprocVTU(
    'Ariane',
    SOL.transpose(),
    nb_ref=np.array([4, 4, 0]),
    Flag=np.array([True, False, False])))