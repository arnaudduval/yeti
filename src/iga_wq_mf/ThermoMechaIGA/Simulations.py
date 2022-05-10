# SIMULATIONS
# ==============================================
# Author : Joaquin CORNEJO

# Python libraries
from datetime import datetime
import scipy
from scipy import sparse as sp
import numpy as np
import os, sys, tracemalloc
import time

# My libraries
from pysrc.lib.create_geomdl import create_geometry
from pysrc.lib.fortran_mf_iga import fortran_mf_iga
from pysrc.lib.fortran_mf_wq import fortran_mf_wq
from pysrc.lib.physics import (powden_cube, 
                            powden_prism,
                            powden_thickring, 
                            powden_rotring, 
                            temperature
)
from pysrc.lib.create_model import write_text_file

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder):
    os.mkdir(folder)

# Enable and disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

def run_simulation(degree, cuts, geometry_case, funpowden, funtemp, isiga, 
                    method_list, iscg, thermalblockedboundaries= None):
    
    # Direct solver
    # -------------
    tracemalloc.start()

    # Define geometry 
    modelGeo = create_geometry(degree, cuts, geometry_case)

    # Create thermal model object
    if isiga: Model1 = fortran_mf_iga(modelGeo, thermalblockedboundaries= thermalblockedboundaries)
    else: Model1 = fortran_mf_wq(modelGeo, thermalblockedboundaries= thermalblockedboundaries)
    del modelGeo

    # Block boundaries
    dof = Model1._thermal_dof
    
    # Assemble conductivity matrix K
    start = time.time()
    K2solve = Model1.eval_conductivity_matrix(indi=dof, indj=dof)
    stop = time.time()
    time_assembly = stop - start

    # Assemble source vector F
    try: 
        dod = Model1._thermal_dod
        T_cp, Td = Model1.MSE_ControlPoints(funtemp)
        F2solve = Model1.eval_source_vector(funpowden, indi=dof, indj=dod, Td=Td)
        del dod, Td, T_cp
    except: 
        F2solve = Model1.eval_source_vector(funpowden, indi=dof)
    
    # Solve system
    start = time.time()
    # LU = sp.linalg.splu(K2solve)
    # sol_direct = LU.solve(F2solve)
    sol_direct = scipy.linalg.solve(K2solve.todense(), F2solve)
    stop = time.time()
    time_direct = stop - start
    time.sleep(1)
    _, memory_direct = tracemalloc.get_traced_memory()
    memory_direct /= 1024*1024
    # del K2solve, F2solve, LU, Model1, dof
    del K2solve, F2solve, Model1, dof

    # Recursive solver 
    # ----------------
    tracemalloc.clear_traces()

    # Define geometry 
    modelGeo = create_geometry(degree, cuts, geometry_case)

    # Create thermal model object
    if isiga: Model1 = fortran_mf_iga(modelGeo, thermalblockedboundaries= thermalblockedboundaries)
    else: Model1 = fortran_mf_wq(modelGeo, thermalblockedboundaries= thermalblockedboundaries)
    del modelGeo

    # Block boundaries
    dof = Model1._thermal_dof

    # Assemble source vector F
    try: 
        dod = Model1._thermal_dod
        T_cp, Td = Model1.MSE_ControlPoints(funtemp)
        F2solve = Model1.eval_source_vector(funpowden, indi=dof, indj=dod, Td=Td)
        del dod, Td, T_cp
    except: 
        F2solve = Model1.eval_source_vector(funpowden, indi=dof)

    time.sleep(1)
    _, memory_iter_base = tracemalloc.get_traced_memory()
    memory_iter_base /= 1024*1024

    # Only compute time to prepare method before iterations
    time_noiter, memory_noiter = [], []
    epsilon  = 1e-14 
    iterations = 0
    tracemalloc.clear_traces()
    for name in method_list:
        start = time.time()
        Model1.mf_conj_grad(F2solve, dof, iterations, epsilon, name, sol_direct, iscg)
        stop = time.time()
        time_noiter_t = stop - start 
        time.sleep(1)
        _, memory_noiter_t = tracemalloc.get_traced_memory()
        memory_noiter_t /= 1024*1024
        tracemalloc.clear_traces()

        # Save data
        time_noiter.append(time_noiter_t)
        memory_noiter.append(memory_noiter_t)

    # With and without preconditioner
    # Initialize
    time_iter, residue, error, memory_iter = [], [], [], []
    epsilon  = 1e-14 
    iterations = 100
    tracemalloc.clear_traces()
    for name in method_list:
        start = time.time()
        _, residue_t, error_t = Model1.mf_conj_grad(F2solve, dof, iterations, epsilon, name, sol_direct, iscg)
        stop = time.time()
        time_iter_t = stop - start 
        time.sleep(1)
        _, memory_iter_sup = tracemalloc.get_traced_memory()
        memory_iter_sup /= 1024*1024
        tracemalloc.clear_traces()
    
        # Save data
        time_iter.append(time_iter_t)
        residue.append(residue_t)
        error.append(error_t)
        memory_iter.append(memory_iter_base+memory_iter_sup)

    del Model1, F2solve, dof
    tracemalloc.stop()
    
    output = {"TimeAssembly": time_assembly, "TimeDirect": time_direct, "MemDirect": memory_direct, 
                "TimeNoIter":time_noiter, "TimeIter": time_iter, "Res": residue, 
                "Error": error, "MemNoIter": memory_noiter, "MemIter": memory_iter}

    return output

# Constants
BlockedBoundaries = [[1, 1], [1, 1], [1, 1]]
for CUTS in range(3, 5):
    for IS_IGA_GALERKIN in [False]:
        for GEOMETRY_CASE in range(3):

            if IS_IGA_GALERKIN: is_cg_list = [True]
            else: is_cg_list = [False]
        
            for IS_CG in is_cg_list:
                for DEGREE in range(3, 5):
                    # Print current time
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    print("Current Time =", current_time)
                    
                    # Get file name
                    if GEOMETRY_CASE == 0: 
                        txtname = 'CB' 
                        funpow = powden_cube 
                        funtemp = None
                    elif GEOMETRY_CASE == 1: 
                        txtname = 'VB' 
                        funpow = powden_prism 
                        funtemp = None
                    elif GEOMETRY_CASE == 2: 
                        txtname = 'TR' 
                        funpow = powden_thickring 
                        funtemp = None
                    elif GEOMETRY_CASE == 3: 
                        txtname = 'RQA' 
                        funpow = powden_rotring 
                        funtemp = temperature
                    else: raise Warning('Geometry does not exist')
                    txtname += '_p' + str(DEGREE) + '_nbel' + str(2**CUTS)

                    if IS_IGA_GALERKIN: txtname += '_IGAG'
                    else: txtname += '_IGAWQ'

                    if IS_CG: txtname += '_CG'
                    else: txtname += '_BiCG'

                    print([DEGREE, 2**CUTS, txtname])

                    # Run simulation
                    blockPrint()
                    method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]
                    inputs_export = run_simulation(DEGREE, CUTS, GEOMETRY_CASE, funpow, funtemp, IS_IGA_GALERKIN, 
                                    method_list, IS_CG, thermalblockedboundaries= BlockedBoundaries)
                    enablePrint()

                    # Export results
                    txtname = folder + txtname + '.txt'
                    write_text_file(txtname, method_list, inputs_export)

