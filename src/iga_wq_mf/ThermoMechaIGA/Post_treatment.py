# SIMULATIONS
# ==============================================
# Author : Joaquin CORNEJO

# Python libraries
import matplotlib.pyplot as plt 
import numpy as np
import os

# My libraries
from pysrc.lib.create_model import read_text_file, plot_iterative_solver

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
        for GEOMETRY_CASE in ['CB', 'VB', 'TR', 'RQA']:
            
            if IS_IGA_GALERKIN: is_cg_list = [True]
            else: is_cg_list = [True, False]
        
            for IS_CG in is_cg_list:
                for DEGREE in range(3, 7):

                    # Recreate file name
                    txtname = GEOMETRY_CASE + '_p' + str(DEGREE) + '_nbel' + str(2**CUTS)
                    if IS_IGA_GALERKIN: txtname += '_IGAG'
                    else: txtname += '_IGAWQ'
                    if IS_CG: txtname += '_CG'
                    else: txtname += '_BiCG'

                    # Extract results
                    txtname = folder + txtname
                    try: 
                        inputs = read_text_file(txtname + '.txt')
                        method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]
                        plot_iterative_solver(txtname, inputs, method_list)
                    except: pass



