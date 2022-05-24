# SIMULATIONS
# ==============================================
# Author : Joaquin CORNEJO

# Python libraries
import os
import matplotlib.pyplot as plt 
import numpy as np

# My libraries
from pysrc.lib.create_model import read_text_file, plot_iterative_solver

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder):
    os.mkdir(folder)

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                    '#f781bf', '#a65628', '#984ea3',
                    '#999999', '#e41a1c', '#dede00']

marker_cycle = ['o', 'v', 'X', 's', '+']

def compute_acceleration(inputs, method_list, ax, marker, color, label=None, PosMethod=0): 
    "We assume that the comparison is between the first method with the others"

    # Define inputs
    residue = inputs["Res"]
    time_noiter = inputs["TimeNoIter"]
    time_iter = inputs["TimeIter"]

    # Define outputs
    acceleration = []

    # Define reference
    residue_t = np.asarray(residue[PosMethod])
    res0 = residue_t[len(residue_t[residue_t>0.0])-1]
    time0 = time_iter[PosMethod]

    # Define others
    for method in range(PosMethod+1, len(method_list)): 
        # Raw time and residue (with all iterations needed)
        time_noiter_t = time_noiter[method]
        time_iter_t = time_iter[method]
        residue_t = np.asarray(residue[method])

        # Find how many iterations were necessary to find res0
        nbIter_t1 = len(residue_t[residue_t>0.0])-1
        nbIter_t2 = 0
        for _ in residue_t:
            if _ >= res0: nbIter_t2 += 1
            else:  break

        # Find time only iterations
        time_onlyIter = time_iter_t - time_noiter_t

        # Scale time to find the same error
        time_t = time_noiter_t + time_onlyIter/nbIter_t1*nbIter_t2

        # Define acceleration
        acceleration.append(time0/time_t)

    # Plots results
    positions = tuple([i for i in range(len(method_list)-PosMethod-1)])
    labels = tuple(method_list[PosMethod+1:])

    ax.plot(np.arange(len(positions)), acceleration, 
            color=color, marker=marker, linestyle="None", label=label)
    ax.plot(np.arange(len(positions)), np.ones(len(positions)), "k--")
    ax.set_xlabel('Methods', fontsize=16)
    ax.set_ylabel('Acceleration', fontsize=16)
    ax.set_xticks(positions, labels)
    ax.tick_params(axis='x', labelsize=13)
    ax.tick_params(axis='y', labelsize=13)
    ax.legend(bbox_to_anchor= (1.05, 1.0), loc= 'upper left')

    return 

# # ====================
# # POST TREATEMENT
# # ====================
# for CUTS in range(4, 6):
#     for IS_IGA_GALERKIN in [False]:
#         for GEOMETRY_CASE in ['CB', 'VB', 'TR', 'RQA']:

#             if IS_IGA_GALERKIN: is_cg_list = [True]
#             else: is_cg_list = [True, False]
        
#             for IS_CG in is_cg_list:
#                 for DEGREE in range(3, 7):

#                     # Recreate file name
#                     txtname = GEOMETRY_CASE + '_p' + str(DEGREE) + '_nbel' + str(2**CUTS)
#                     if IS_IGA_GALERKIN: txtname += '_IGAG'
#                     else: txtname += '_IGAWQ'
#                     if IS_CG: txtname += '_CG'
#                     else: txtname += '_BiCG'

#                     # Extract results
#                     txtname = folder + txtname
#                     try: 
#                         inputs = read_text_file(txtname + '.txt')
#                         method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]
#                         plot_iterative_solver(txtname, inputs, method_list)
#                     except: pass

for GEOMETRY_CASE in ['CB', 'VB', 'TR', 'RQA']:
    fig, ax = plt.subplots()
    for i, CUTS in enumerate(range(4, 6)):
        color = CB_color_cycle[i]
        for j, DEGREE in enumerate(range(3, 7)):
            marker = marker_cycle[j]
            for IS_IGA_GALERKIN in [False]:
            
                if IS_IGA_GALERKIN: is_cg_list = [True]
                else: is_cg_list = [True]
            
                for IS_CG in is_cg_list:
                    # Recreate file name
                    txtname = GEOMETRY_CASE + '_p' + str(DEGREE) + '_nbel' + str(2**CUTS)
                    if IS_IGA_GALERKIN: txtname += '_IGAG'
                    else: txtname += '_IGAWQ'
                    if IS_CG: txtname += '_CG'
                    else: txtname += '_BiCG'

                    # Extract results
                    txtname = folder + txtname
                    inputs = read_text_file(txtname + '.txt')
                    method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"] # !!!!!!! This must be equal to simulations
                    label = 'nbel = ' + str(2**CUTS) + ", p = " + str(DEGREE)
                    compute_acceleration(inputs, method_list, 
                                        ax, marker, color, label=label, 
                                        PosMethod= 0)
                    
    plt.tight_layout()
    plt.savefig(folder + GEOMETRY_CASE + '_Acc' +'.png')

