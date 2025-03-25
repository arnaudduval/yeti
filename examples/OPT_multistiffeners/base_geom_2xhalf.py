"""
TODO
"""

import numpy as np
import scipy.sparse as sp

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp
from yeti_iga.preprocessing.igaparametrization import IGAmanip as manip
from yeti_iga.coupling.cplgmatrix import cplg_matrix


iga_model = IGAparametrization(filename='base_geom_2xhalf')
iga_model._NBPINT[
    np.where(iga_model._ELT_TYPE == 'U00')] = 6

# Refine model
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

nb_deg[:, 0] = np.array([1, 2, 0])
nb_deg[:, 1] = np.array([1, 2, 0])
additional_knots = {'patches': np.array([0]),
                    '1': np.array([0.43, 0.57]),
                    '2': np.array([0.1, 0.3]),
                    '3': np.array([])}
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
additional_knots = {'patches': np.array([1]),
                    '1': np.array([0.43, 0.57]),
                    '2': np.array([0.1, 0.3]),
                    '3': np.array([])}
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

# Add loading on a single element
iga_model._indDLoad = np.array([ np.array([6,7]), np.array([6,7]), np.array([18,19]), np.array([18,19]) ])
iga_model._JDLType = np.array([62,63,62,63])
iga_model._ADLMAG = np.array([1.190, 7.143, -1.190, 7.143])
iga_model._load_target_nbelem = np.array([2, 2, 2, 2])
iga_model._nb_load = 4
iga_model._indDLoad_flat = np.array([], dtype=np.intp)
for load in iga_model._indDLoad:
    iga_model._indDLoad_flat = np.hstack((iga_model._indDLoad_flat, load))


# Refinement to get elements with approx the same size
additional_knots = {'patches': np.array([0]),
                    '1': np.array([0.1075, 0.215, 0.3325, 0.6775, 0.785, 0.8925]),
                    '2': np.array([0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
                    '3': np.array([])}
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)
additional_knots = {'patches': np.array([1]),
                    '1': np.array([0.1075, 0.215, 0.3325, 0.6775, 0.785, 0.8925]),
                    '2': np.array([0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
                    '3': np.array([])}
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

# OTHER REFINEMENTS
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
# Curves
nb_ref[0, 2:10] = 3

# Lagrange patches
nb_deg[0, 10] = 1
nb_deg[0, 11] = 0
nb_deg[0, 12] = 1
nb_deg[0, 13] = 0

nb_ref[0, 10] = 3
nb_ref[0, 11] = 3
nb_ref[0, 12] = 3
nb_ref[0, 13] = 3
iga_model.refine(nb_ref, nb_deg)

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

Cdata,Crow,Ccol = cplg_matrix( *iga_model.get_inputs4cplgmatrix() )
Cside = sp.coo_matrix((Cdata,(Crow,Ccol)), shape=(iga_model._nb_dof_tot, iga_model._nb_dof_tot),
                      dtype='float64').tocsc()
Ctot  = Cside + Cside.transpose()
del Cdata,Crow,Ccol,Cside
C2solve = Ctot[idof,:][:,idof] * K2solve.max()


# Resolution
x = sp.linalg.spsolve(K2solve + C2solve, Fbsolve)

# Solution reconstruction
SOL, u = rsol.reconstruction(**iga_model.get_inputs4solution(x))
pp.generatevtu(*iga_model.get_inputs4postprocVTU(
    'Ariane',
    SOL.transpose(),
    nb_ref=np.array([4, 4, 0]),
    Flag=np.array([True, False, False])))