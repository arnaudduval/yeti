"""
TODO
"""

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import csv
from datetime import datetime

from yeti_iga.preprocessing.igaparametrization import IGAparametrization
from yeti_iga.stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import yeti_iga.reconstructionSOL as rsol
import yeti_iga.postprocessing.postproc as pp
from yeti_iga.preprocessing.igaparametrization import IGAmanip as manip
from yeti_iga.coupling.cplgmatrix import cplg_matrix
from yeti_iga.preprocessing.igaparametrization import OPTmodelling

from stiffened_cylinder import TEMPLATE_INP_FILE, NB_STIFFENERS_BY_HALF, \
    set_X_stiffeners, set_H_stiffeners, write_inp_file, write_nb_file
from tools import len_stiffeners, grad_len_stiffeners, smoothmin




ALLOWED_VOLUME_VAR = 0.05
ALLOWED_VOLUME = 1.7e6
ALLOWED_COMPLIANCE = 0.8e6
MIN_STIFF_LEN = 100.
MIN_P_ORDER = 20

OPTIM = True

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

write_inp_file(f'stiff_geom_2xhalf_{timestamp}.inp')
write_nb_file(f'stiff_geom_2xhalf_{timestamp}.NB')


# write NB file
with open(f'stiff_geom_2xhalf_{timestamp}.NB', 'w') as nbfile:
    nbfile.write('*Dimension\n')
    nbfile.write('2,2,1,1,1,1,1,1,1,1,1,1,1,1,3,3')
    for istiff in range(NB_STIFFENERS_BY_HALF):
        nbfile.write(',2,2,1,1,1,1,1,1,1,1,1,1,1,1')
    nbfile.write('\n')

    nbfile.write('*Number of CP by element\n')
    nbfile.write('6,6,2,2,2,2,2,2,2,2,1,1,1,1,12,12')
    for istiff in range(NB_STIFFENERS_BY_HALF):
        nbfile.write(',4,4,2,2,2,2,2,2,2,2,1,1,1,1')
    nbfile.write('\n')

    nbfile.write('*Number of patch\n')
    nbfile.write(f'{16+NB_STIFFENERS_BY_HALF*14}\n')

    nbfile.write('*Total number of element\n')
    nbfile.write(f'{20+NB_STIFFENERS_BY_HALF*14}\n')

    nbfile.write('*Number of element by patch\n')
    nbfile.write('2,2,1,1,1,1,1,1,1,1,1,1,1,1,2,2')
    for istiff in range(NB_STIFFENERS_BY_HALF):
        nbfile.write(',1,1,1,1,1,1,1,1,1,1,1,1,1,1')
    nbfile.write('\n')

    nbfile.write('*Patch(1)\n7\n0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(2)\n7\n0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(3)\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(4)\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(5)\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(6)\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(7)\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(8)\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(9)\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(10)\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(11)\n2\n0.0, 1.0\n')
    nbfile.write('*Patch(12)\n2\n0.0, 1.0\n')
    nbfile.write('*Patch(13)\n2\n0.0, 1.0\n')
    nbfile.write('*Patch(14)\n2\n0.0, 1.0\n')
    nbfile.write('*Patch(15)\n7\n0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
    nbfile.write('*Patch(16)\n7\n0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
    for istiff in range(NB_STIFFENERS_BY_HALF):
        nbfile.write(f'*Patch({17+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write(f'*Patch({18+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')

        nbfile.write(f'*Patch({19+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write(f'*Patch({20+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write(f'*Patch({21+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write(f'*Patch({22+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write(f'*Patch({23+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write(f'*Patch({24+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write(f'*Patch({25+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write(f'*Patch({26+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')

        nbfile.write(f'*Patch({27+istiff*14})\n2\n0.0, 1.0\n')
        nbfile.write(f'*Patch({28+istiff*14})\n2\n0.0, 1.0\n')
        nbfile.write(f'*Patch({29+istiff*14})\n2\n0.0, 1.0\n')
        nbfile.write(f'*Patch({30+istiff*14})\n2\n0.0, 1.0\n')

    nbfile.write('*Jpqr\n')
    nbfile.write('2,1\n2,1\n1\n1\n1\n1\n1\n1\n1\n1\n0\n0\n0\n0\n2,1,1\n2,1,1\n')
    for istiff in range(NB_STIFFENERS_BY_HALF):
        nbfile.write('1,1\n1,1\n1\n1\n1\n1\n1\n1\n1\n1\n0\n0\n0\n0\n')
    nbfile.write('*Nijk\n')
    nbfile.write('1, 3, 2\n2, 4, 2\n3, 3, 2\n4, 4, 2\n5, 2\n6, 2\n7, 2\n8, 2\n9, 2\n10, 2\n11, 2\n12, 2\n13, 1\n14, 1\n15, 1\n16, 1\n17, 3, 2, 2\n18, 4, 2, 2\n19, 3, 2, 2\n20, 4, 2, 2\n')
    for istiff in range(NB_STIFFENERS_BY_HALF):
        nbfile.write(f'{21+istiff*14}, 2, 2\n')
        nbfile.write(f'{22+istiff*14}, 2, 2\n')

        nbfile.write(f'{23+istiff*14}, 2\n')
        nbfile.write(f'{24+istiff*14}, 2\n')
        nbfile.write(f'{25+istiff*14}, 2\n')
        nbfile.write(f'{26+istiff*14}, 2\n')
        nbfile.write(f'{27+istiff*14}, 2\n')
        nbfile.write(f'{28+istiff*14}, 2\n')
        nbfile.write(f'{29+istiff*14}, 2\n')
        nbfile.write(f'{30+istiff*14}, 2\n')

        nbfile.write(f'{31+istiff*14}, 1\n')
        nbfile.write(f'{32+istiff*14}, 1\n')
        nbfile.write(f'{33+istiff*14}, 1\n')
        nbfile.write(f'{34+istiff*14}, 1\n')

    nbfile.write('*Weight\n')
    nbfile.write('1, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0\n')
    nbfile.write('2, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5\n')
    nbfile.write('3, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0\n')
    nbfile.write('4, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5\n')
    nbfile.write('5, 1.0, 1.0\n')
    nbfile.write('6, 1.0, 1.0\n')
    nbfile.write('7, 1.0, 1.0\n')
    nbfile.write('8, 1.0, 1.0\n')
    nbfile.write('9, 1.0, 1.0\n')
    nbfile.write('10, 1.0, 1.0\n')
    nbfile.write('11, 1.0, 1.0\n')
    nbfile.write('12, 1.0, 1.0\n')
    nbfile.write('13, 1.0\n')
    nbfile.write('14, 1.0\n')
    nbfile.write('15, 1.0\n')
    nbfile.write('16, 1.0\n')
    nbfile.write('17, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0\n')
    nbfile.write('18, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5\n')
    nbfile.write('19, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0\n')
    nbfile.write('20, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5\n')
    for istiff in range(NB_STIFFENERS_BY_HALF):
        nbfile.write(f'{21+istiff*14}, 1.0, 1.0, 1.0, 1.0\n')
        nbfile.write(f'{22+istiff*14}, 1.0, 1.0, 1.0, 1.0\n')
        nbfile.write(f'{23+istiff*14}, 1.0, 1.0\n')
        nbfile.write(f'{24+istiff*14}, 1.0, 1.0\n')
        nbfile.write(f'{25+istiff*14}, 1.0, 1.0\n')
        nbfile.write(f'{26+istiff*14}, 1.0, 1.0\n')
        nbfile.write(f'{27+istiff*14}, 1.0, 1.0\n')
        nbfile.write(f'{28+istiff*14}, 1.0, 1.0\n')
        nbfile.write(f'{29+istiff*14}, 1.0, 1.0\n')
        nbfile.write(f'{30+istiff*14}, 1.0, 1.0\n')
        nbfile.write(f'{31+istiff*14}, 1.0\n')
        nbfile.write(f'{32+istiff*14}, 1.0\n')
        nbfile.write(f'{33+istiff*14}, 1.0\n')
        nbfile.write(f'{34+istiff*14}, 1.0\n')



iga_model = IGAparametrization(filename=f'stiff_geom_2xhalf_{timestamp}')
iga_model._NBPINT[
    np.where(iga_model._ELT_TYPE == 'U00')] = 6

# Refine model
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

print('Raffinement initial + découpage patch 1')
nb_deg[:, 0] = np.array([0, 1, 0])
nb_deg[:, 1] = np.array([0, 1, 0])
additional_knots = {'patches': np.array([0]),
                    '1': np.array([0.43, 0.57]),
                    '2': np.array([0.1, 0.3]),
                    '3': np.array([])}
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

print('découpage patch 2')
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

additional_knots = {"patches": np.array([]),
                    "1": np.array([]), "2": np.array([]), "3": np.array([])}

# Refinement to get elements with approx the same size
print("Uniformisation maillage patch 1")
additional_knots = {'patches': np.array([0]),
                    '1': np.array([0.1075, 0.215, 0.3325, 0.6775, 0.785, 0.8925]),
                    '2': np.array([0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
                    '3': np.array([])}
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)
print("Uniformisation maillage patch 2")
additional_knots = {'patches': np.array([1]),
                    '1': np.array([0.1075, 0.215, 0.3325, 0.6775, 0.785, 0.8925]),
                    '2': np.array([0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
                    '3': np.array([])}
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

def movestiffeners(coords0, igapara, var):
    igapara._COORDS[:, :] = coords0[:, :]

    assert var.shape[0] == NB_STIFFENERS_BY_HALF*4

    for istiff in range(NB_STIFFENERS_BY_HALF):
        igapara._COORDS[0, igapara._indCPbyPatch[16+istiff*14][:]-1] = var[[0 + istiff*4, 1 + istiff*4, 0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[17+istiff*14][:]-1] = var[[0 + istiff*4, 1 + istiff*4, 0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[18+istiff*14][[0,1]]-1] = var[[0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[20+istiff*14][[0,1]]-1] = var[[0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[22+istiff*14][[0,1]]-1] = var[[0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[24+istiff*14][[0,1]]-1] = var[[0 + istiff*4, 1 + istiff*4]]

        igapara._COORDS[1, igapara._indCPbyPatch[16+istiff*14][:]-1] = var[[2 + istiff*4, 3 + istiff*4, 2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[17+istiff*14][:]-1] = var[[2 + istiff*4, 3 + istiff*4, 2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[18+istiff*14][[0,1]]-1] = var[[2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[20+istiff*14][[0,1]]-1] = var[[2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[22+istiff*14][[0,1]]-1] = var[[2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[24+istiff*14][[0,1]]-1] = var[[2 + istiff*4, 3 + istiff*4]]

    return None

# x0 = set_X_stiffeners(NB_STIFFENERS_BY_HALF)
x0 = set_H_stiffeners(NB_STIFFENERS_BY_HALF)

movestiffeners(iga_model._COORDS.copy(), iga_model, x0)

# OTHER REFINEMENTS For optim problem
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

nb_ref[:, 0] = [1, 1, 0]
nb_ref[:, 1] = [1, 1, 0]

# Hull
nb_deg[:, 14] = [0, 1, 1]
nb_deg[:, 15] = [0, 1, 1]

# Curves 2 halfs
nb_ref[0, 2:10] = 3
# Lagrange patches 2 halfs
nb_deg[0, 10] = 1
nb_deg[0, 11] = 0
nb_deg[0, 12] = 1
nb_deg[0, 13] = 0

nb_ref[0, 10] = 3
nb_ref[0, 11] = 3
nb_ref[0, 12] = 3
nb_ref[0, 13] = 3

# Stiffeners list
for istiff in range(NB_STIFFENERS_BY_HALF):
    nb_deg[:, 16+istiff*14] = [1, 1, 0]
    nb_deg[:, 17+istiff*14] = [1, 1, 0]
    nb_ref[:, 16+istiff*14] = [4, 2, 0]
    nb_ref[:, 17+istiff*14] = [4, 2, 0]

    nb_ref[0, 18+istiff*14] = 4
    nb_ref[0, 19+istiff*14] = 4
    nb_ref[0, 20+istiff*14] = 4
    nb_ref[0, 21+istiff*14] = 4
    nb_ref[0, 22+istiff*14] = 4
    nb_ref[0, 23+istiff*14] = 4
    nb_ref[0, 24+istiff*14] = 4
    nb_ref[0, 25+istiff*14] = 4

    nb_deg[0, 26+istiff*14] = 1
    nb_deg[0, 27+istiff*14] = 0
    nb_deg[0, 28+istiff*14] = 1
    nb_deg[0, 29+istiff*14] = 0

    nb_ref[0, 26+istiff*14] = 3
    nb_ref[0, 27+istiff*14] = 3
    nb_ref[0, 28+istiff*14] = 3
    nb_ref[0, 29+istiff*14] = 3

# END OF OTHER REFINEMENTS

# Define optim problem
opt_pb = OPTmodelling(iga_model, NB_STIFFENERS_BY_HALF*4, movestiffeners,
                    nb_degreeElevationByDirection=nb_deg,
                    nb_refinementByDirection=nb_ref)

# List of patches for volume computation
listpatch = np.zeros(opt_pb._coarseParametrization._nb_patch, dtype=np.intp)
listpatch[0] = 1
listpatch[1] = 1
for istiff in range(NB_STIFFENERS_BY_HALF):
    listpatch[16+istiff*14] = 1
    listpatch[17+istiff*14] = 1

if OPTIM:
    # Optimization

    # Referenc values for volume ad compliance
    V0 = opt_pb.compute_volume(x0, listpatch)
    C0 = opt_pb.compute_compliance_discrete(x0)

    # Optimization fonction (relative compliance and volume variation)
    def var_vol(xk):
        return ALLOWED_VOLUME_VAR - abs(opt_pb.compute_volume(xk, listpatch) - V0) / V0

    def grad_var_vol(xk):
        vk = opt_pb.compute_volume(xk, listpatch)
        return np.sign(vk - V0)*opt_pb.compute_gradVolume_AN(xk, listpatch)/V0

    def imposed_volume(xk):
        vk = opt_pb.compute_volume(xk, listpatch)
        return (vk - ALLOWED_VOLUME)/ALLOWED_VOLUME

    def grad_imposed_volume(xk):
        return opt_pb.compute_gradVolume_AN(xk, listpatch)/ALLOWED_VOLUME

    def volume(xk):
        vk = opt_pb.compute_volume(xk, listpatch)
        return vk/V0

    def grad_volume(xk):
        return opt_pb.compute_gradVolume_AN(xk, listpatch)/V0

    def compliance(xk):
        return opt_pb.compute_compliance_discrete(xk)/C0

    def grad_compliance(xk):
        return (opt_pb.compute_gradCompliance_AN(xk)+ opt_pb.compute_gradCompliance_cplgOnly_AN(xk))/C0

    def imposed_compliance(xk):
        return (opt_pb.compute_compliance_discrete(xk) - ALLOWED_COMPLIANCE)/ALLOWED_COMPLIANCE

    def grad_imposed_compliance(xk):
        return (opt_pb.compute_gradCompliance_AN(xk)+ opt_pb.compute_gradCompliance_cplgOnly_AN(xk))/ALLOWED_COMPLIANCE

    def min_stiff_len(xk):
        lengths = len_stiffeners(xk)
        min_len, _ = smoothmin(lengths, p=MIN_P_ORDER)
        # print(f'{min_len = }')
        return min_len - MIN_STIFF_LEN

    def grad_min_stiff_len(xk):
        lengths = len_stiffeners(xk)
        min_len, grad_min_len = smoothmin(lengths, p=MIN_P_ORDER)
        return grad_min_len @ grad_len_stiffeners(xk)

    iopt = 0

    compliances = []
    volumes = []

    def save_xk(xk):
        global iopt
        print((f'\nIteration {iopt:03}'))
        pp.generatevtu(*opt_pb.coarse_parametrization.get_inputs4postprocVTU(
            f'opt-coarse{iopt:02}', np.zeros_like(opt_pb.coarse_parametrization.coords),
            nb_ref=4*np.array([1, 1, 0]),Flag=np.array([False]*3)))
        opt_pb.coarse_parametrization.generate_vtk4controlMeshVisu(f'opt-coarse{iopt:02}',0)

        sol, _ = rsol.reconstruction(
                **opt_pb.fine_parametrization.get_inputs4solution(opt_pb.save_sol_fine))
        pp.generatevtu(*opt_pb.fine_parametrization.get_inputs4postprocVTU(
                        f'opt-fine{iopt:02}',  sol.transpose(),
                        nb_ref=4*np.array([1, 1, 0]),
                        Flag=np.array([True, False, False])))
        compliances.append(opt_pb.compute_compliance_discrete(xk))
        volumes.append(opt_pb.compute_volume(xk, listpatch))
        iopt += 1

    save_xk(x0)

    # Bounds for design variables
    bounds = ((0., 1.),) * NB_STIFFENERS_BY_HALF*4
    constraints = ({'type': 'eq', 'fun': imposed_volume, 'jac': grad_imposed_volume},
                   {'type': 'ineq', 'fun': min_stiff_len, 'jac': grad_min_stiff_len})
    # constraints = ({'type': 'eq', 'fun': imposed_compliance, 'jac': grad_imposed_compliance},
    #                {'type': 'ineq', 'fun': min_stiff_len, 'jac': grad_min_stiff_len})
    res = minimize(compliance, x0, method='SLSQP',
                jac=grad_compliance, bounds=bounds,
                constraints=constraints, callback=save_xk,
                tol=5.E-3)


    print('Optimization succes : ', res['success'])
    print(res['message'])
    # print('Result : ', res['x'])
    save_xk(res['x'])
    print('Objective function value : ', res['fun'])
    print('# of evaluations of objective function : ', res['nfev'])
    print('# of evaluations of jacobian : ', res['njev'])
    v_f = opt_pb.compute_volume(res['x'], listpatch)
    c_f = opt_pb.compute_compliance_discrete(res['x'])
    print('Volume: ', V0, '->', v_f, (100.*(v_f - V0)/V0), ' %')
    print('Compliance: ', C0, '->',  c_f, (100.*(c_f - C0)/C0), ' %')
    # with open(f'optim_imposed_comp_{ALLOWED_COMPLIANCE}.csv', mode='w', newline='') as f:
    with open(f'optim_imposed_volume_{ALLOWED_VOLUME}.csv', mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Volume', 'Compliance'])
        for v, c in zip(volumes, compliances):
            writer.writerow([v, c])
    print(volumes)
    print(compliances)


else:
    output_file = f'results_{timestamp}.csv'
    nb_design_variables = NB_STIFFENERS_BY_HALF*4
    def compute_volume_and_compliance(x):
        # NOTE same as save_kk ??
        volume = opt_pb.compute_volume(x, listpatch)
        compliance = opt_pb.compute_compliance_discrete(x)
        print(volume, compliance)
        return volume, compliance

    def save_results(x, volume, compliance, writer):
        """Write results in CSV file"""
        row = [volume, compliance] + list(x)
        writer.writerow(row)

    volumes = []
    compliances = []
    header = ["volume", "compliance"] + [f"x{i}" for i in range(nb_design_variables)]
    with open(output_file, mode='w', newline='')as f:
        writer = csv.writer(f)
        writer.writerow(header)

        x_ref = set_X_stiffeners(NB_STIFFENERS_BY_HALF)
        try:
            v, c = compute_volume_and_compliance(x_ref)
            volumes.append(v)
            compliances.append(c)
            save_results(x_ref, v, c, writer)
        except Exception as e:
            print(f"Erreur ignorée : {e}")

        for i in range(10):
            print(i)
            try:
                min_len = -1.
                while min_len < MIN_STIFF_LEN:
                    x_rand = np.random.rand(nb_design_variables)
                    lengths = len_stiffeners(x_rand)
                    print(f'Lg : {lengths}')
                    min_len, _ = smoothmin(lengths, p=MIN_P_ORDER)
                    print(f"Lg mini : {min_len}, {'pas OK' if min_len < MIN_STIFF_LEN else 'OK'}")
                v, c = compute_volume_and_compliance(x_rand)
                volumes.append(v)
                compliances.append(c)
                save_results(x_rand, v, c, writer)
            except Exception as e:
                print(f"Erreur ignorée : {e}")
                continue

