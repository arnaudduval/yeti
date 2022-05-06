# SIMULATIONS
# ==============================================
# Author : Joaquin CORNEJO

# Python libraries
import matplotlib.pyplot as plt 
import numpy as np
import os

# My libraries
from pysrc.lib.create_model import read_text_file

def plot_iterative_solver(filename, inputs, method_list, extension ='.png'):
    # Define color
    CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                    '#f781bf', '#a65628', '#984ea3',
                    '#999999', '#e41a1c', '#dede00']

    # # Define inputs
    # time_assembly = inputs["TimeAssembly"]
    # time_direct = inputs["TimeDirect"]
    # memory_direct = inputs["MemDirect"]
    # time_iter = inputs["TimeIter"]
    residue = inputs["Res"]
    error = inputs["Error"]
    # memory_iter = inputs["MemIter"]

    # Select important values
    tol = 1.e-16
    
    # Set figure parameters
    fig, axs = plt.subplots(1, 2)
    fig.set_size_inches(12.5, 6.5)

    for i, pcgmethod in enumerate(method_list):
        error_method = np.asarray(error[i])
        residue_method = np.asarray(residue[i])

        error_method = error_method[residue_method>tol]
        residue_method = residue_method[residue_method>tol]
        
        axs[0].plot(np.arange(len(error_method)), error_method*100, label=pcgmethod, color=CB_color_cycle[i])
        axs[1].plot(np.arange(len(residue_method)), residue_method*100, label=pcgmethod, color=CB_color_cycle[i])

    # Set properties
    axs[0].set_yscale("log")
    axs[0].set_xlabel('Number of iterations', fontsize=16)
    axs[0].set_ylabel('Relative error (%)', fontsize=16)
    axs[0].tick_params(axis='x', labelsize=16)
    axs[0].tick_params(axis='y', labelsize=16)
    axs[0].legend()
    axs[0].grid()

    axs[1].set_yscale("log")
    axs[1].set_xlabel('Number of iterations', fontsize=16)
    axs[1].set_ylabel('Relative residue (%)', fontsize=16)
    axs[1].tick_params(axis='x', labelsize=16)
    axs[1].tick_params(axis='y', labelsize=16)
    axs[1].legend()
    axs[1].grid()

    # Save figure
    fig.tight_layout()
    fig.savefig(filename + extension)
    plt.close(fig)

    return

# ====================
# POST TREATEMENT
# ====================
# Constants
GEOMETRY_CASE = 1
IS_IGA_GALERKIN = False

full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'

for CUTS in range(2, 6):
    for IS_IGA_GALERKIN in [False, True]:
        for GEOMETRY_CASE in range(4):
            if IS_IGA_GALERKIN: is_cg_list = [True]
            else: is_cg_list = [True, False]
        
            for IS_CG in is_cg_list:
                for DEGREE in range(3, 7):

                    # Choose name
                    if GEOMETRY_CASE == 0:   txtname = 'CB' 
                    elif GEOMETRY_CASE == 1: txtname = 'VB' 
                    elif GEOMETRY_CASE == 2: txtname = 'TR' 
                    elif GEOMETRY_CASE == 3: txtname = 'RQA'
                    else: raise Warning('Geometry does not exist')

                    # Recreate file name
                    txtname += '_p' + str(DEGREE) + '_nbel' + str(2**CUTS)

                    if IS_IGA_GALERKIN: txtname += '_IGAG'
                    else: txtname += '_IGAWQ'
                
                    if IS_CG: txtname += '_CG'
                    else: txtname += '_BiCG'

                    # Extract results
                    txtname = folder + txtname
                    try: 
                        inputs = read_text_file(txtname + '.txt')

                        # Plot results
                        method_list = ["WP", "FD-C", "FD-TDS", "FD-JM", "FD-TD", "FD-JMS"]
                        # method_list = ["WP", "FD-C", "FD-TDS", "FD-JM"]
                        plot_iterative_solver(txtname, inputs, method_list)
                    except: 
                        # print(txtname)
                        pass



