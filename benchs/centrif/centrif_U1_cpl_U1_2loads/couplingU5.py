import sys
import time
import pickle

import numpy as np
from numpy.lib.arraysetops import intersect1d
import scipy.sparse as sp

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
from yeti_iga.coupling.cplgmatrix import cplg_matrixu5 as cplg_matrixU5
from yeti_iga.coupling.cplgmatrix import cplg_matrixu5_wdatabase as cplg_matrixU5_wdatabase
from yeti_iga.utils import gausspts
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp
from yeti_iga.preprocessing.igaparametrization import IGAmanip as manip


t0 = time.time()

# =============================================================================
# Import file & refine
# =============================================================================
FILENAME = 'centrif_U1_cpl_U1_2loads'

modeleIGA = IGAparametrization(filename=FILENAME)

# Initialise refinement arrays
nb_deg = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, modeleIGA._nb_patch), dtype=np.intp)
additional_knots = {"patches": np.array([]),
                    "1": np.array([]), "2": np.array([]), "3": np.array([])}

# Set refinement arrays

# ---- U5 w\multiple patches - blade & platform 123456 & fillet 1234 + pres

# U1 patch
nb_ref[:, 0] = np.array([2, 2, 2])
nb_deg[:, 0] = np.array([1, 1, 1])
# U1 patch
nb_ref[:, 1] = np.array([2, 2, 2])
nb_deg[:, 1] = np.array([1, 1, 1])
# coupling
nb_ref[:, 2] = np.array([2, 2, 0])
nb_deg[:, 2] = np.array([1, 1, 0])

# Refine
modeleIGA.refine(nb_ref, nb_deg, additional_knots)

# =============================================================================
# Stiffness matrix
# =============================================================================
ndof = modeleIGA._nb_dof_free
idof = modeleIGA._ind_dof_free[:ndof] - 1

# Time for assembling K - start
t_start_k = time.time()

# Build stiffness matrix
data, row, col, Fb = build_stiffmatrix(
    *modeleIGA.get_inputs4system_elemStorage())
print(np.amax(row), np.amin(row), np.amax(col), np.amin(col))
# Print size of K
print('\nNb. DoF K: {}'.format(len(np.unique(row))))

Kside = sp.coo_matrix((data, (row, col)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside+Kside.transpose()
del Kside, data, row, col

# Time for assembling K - end
t_end_k = time.time()

# Print elapsed time for assembly
frmt_time = time.gmtime(t_end_k - t_start_k)
nb_hours = frmt_time.tm_hour
nb_mins = frmt_time.tm_min
nb_secs = frmt_time.tm_sec
str_time = str()
if nb_hours:
    plural_mark = 's' if nb_hours > 1 else ''
    str_time += '{} hour{} '.format(nb_hours, plural_mark)
if nb_mins:
    plural_mark = 's' if nb_mins > 1 else ''
    str_time += '{} minute{} '.format(nb_mins, plural_mark)
if nb_secs:
    plural_mark = 's' if nb_secs > 1 else ''
    str_time += '{} second{} '.format(nb_secs, plural_mark)

print("\nStiffness matrix assembly done in {}.\n".format(str_time[:-1]))

# =============================================================================
# Coupling matrix
# =============================================================================
# Time for assembling C - start
t_start_c = time.time()

# Build coupling matrix
Cdata, Crow, Ccol = cplg_matrixU5(*modeleIGA.get_inputs4cplgmatrixU5())
print(np.amax(Crow), np.amin(Crow), np.amax(Ccol), np.amin(Ccol))
# Cdata, Crow, Ccol = cplg_matrixU5_wdatabase(*modeleIGA.get_inputs4cplgmatrixU5())
Cside = sp.coo_matrix((Cdata, (Crow, Ccol)),
                      shape=(modeleIGA._nb_dof_tot, modeleIGA._nb_dof_tot),
                      dtype='float64').tocsc()

# Print size of C
print('\nNb. DoF C: {}'.format(len(np.unique(Ccol))))

Ctot = Cside + Cside.transpose()
del Cdata, Crow, Ccol, Cside

# Time for assembling C - end
t_end_c = time.time()

# Print elapsed time for assembly
frmt_time = time.gmtime(t_end_c - t_start_c)
nb_hours = frmt_time.tm_hour
nb_mins = frmt_time.tm_min
nb_secs = frmt_time.tm_sec
str_time = str()
if nb_hours:
    plural_mark = 's' if nb_hours > 1 else ''
    str_time += '{} hour{} '.format(nb_hours, plural_mark)
if nb_mins:
    plural_mark = 's' if nb_mins > 1 else ''
    str_time += '{} minute{} '.format(nb_mins, plural_mark)
if nb_secs:
    plural_mark = 's' if nb_secs > 1 else ''
    str_time += '{} second{} '.format(nb_secs, plural_mark)

print("\nCoupling matrix assembly done in {}.\n".format(str_time[:-1]))

# exit()

# =============================================================================
# Solve - monolithic
# =============================================================================
print('\nStart resolution...')
t2 = time.time()
K2solve = Ktot[idof, :][:, idof]
C2solve = Ctot[idof, :][:, idof] * K2solve.max()

# LU = sp.linalg.splu(K2solve + C2solve)
# x = LU.solve(Fb[idof])
x = sp.linalg.spsolve(K2solve + C2solve, Fb[idof])
t3 = time.time()

# =============================================================================
# Print elapsed time for resolution
# =============================================================================
frmt_time = time.gmtime(t3-t2)
nb_hours = frmt_time.tm_hour
nb_mins = frmt_time.tm_min
nb_secs = frmt_time.tm_sec
str_time = str()
if nb_hours:
    plural_mark = 's' if nb_hours > 1 else ''
    str_time += '{} hour{} '.format(nb_hours, plural_mark)
if nb_mins:
    plural_mark = 's' if nb_mins > 1 else ''
    str_time += '{} minute{} '.format(nb_mins, plural_mark)
if nb_secs:
    plural_mark = 's' if nb_secs > 1 else ''
    str_time += '{} second{} '.format(nb_secs, plural_mark)

print("Resolution done in {}.\n".format(str_time[:-1]))

# =============================================================================
# Solution reconstruction
# =============================================================================
SOL, u = rsol.reconstruction(**modeleIGA.get_inputs4solution(x))

# =============================================================================
# Save objects for postproc.
# =============================================================================
# IGAparametrization
file_name = 'modeleIGA'
with open(file_name, 'wb') as file:
    my_pickler = pickle.Pickler(file)
    my_pickler.dump(modeleIGA)
print('\n---> Dumped {}'.format(file_name))
# Solution
file_name = 'solution'
with open(file_name, 'wb') as file:
    my_pickler = pickle.Pickler(file)
    my_pickler.dump(SOL)
print('\n---> Dumped {}'.format(file_name))

# =============================================================================
# Postprocessing
# =============================================================================

# Classical postprocessing
# ------------------------
pp.generatevtu(*modeleIGA.get_inputs4postprocVTU(
    FILENAME,
    SOL.transpose(),
    nb_ref=np.array([2, 2, 2]),
    Flag=np.array([True, False, False])))

# Get solution at interpolating control points and compare with Abaqus reference result
# Tolerance = max 2% on quasi radial component
pt_coords = np.array([12.0, 1.5, 0.0])
ref_aba = np.array([1.60785E-6, -1.94110E-8, 3.15129E-7])

found = False

for idx in range(np.shape(modeleIGA._COORDS)[1]):
    if np.all(modeleIGA._COORDS[:,idx] == pt_coords):
        found = True
        break

if not found:
    sys.exit(-1)

print("yeti solution : ", SOL[idx,:])
print("reference solution : ", ref_aba)
print("error / component : ", np.abs((SOL[idx,:]-ref_aba)/ref_aba))

if np.abs((SOL[idx,:]-ref_aba)/ref_aba)[0] < 0.02:
    sys.exit(0)
else:
    sys.exit(-1)

