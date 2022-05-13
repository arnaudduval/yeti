"""
.. Test of simulation
.. We test a single example of what is meant to do simulation script
.. Joaquin Cornejo 
"""

# Python libraries
import os, sys, tracemalloc
import scipy, numpy as np
import time

# My libraries
from lib.create_geomdl import create_geometry
from lib.fortran_mf_iga import fortran_mf_iga
from lib.fortran_mf_wq import fortran_mf_wq
from lib.physics import (powden_cube, powden_prism,
                        powden_thickring, powden_rotring, 
                        temperature
)
from lib.create_model import write_text_file, read_text_file, plot_iterative_solver

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
                    method_list, iscg, thermalblockedboundaries= None, isDirect=None):
    
    if isDirect is None: 
        doDirect = True
        doIterative = True
    elif isDirect is True: 
        doDirect = True
        doIterative = False
    elif isDirect is False: 
        doDirect = False
        doIterative = True

    # Define solution 
    sol_direct = None
    
    # Direct solver
    # -------------
    if doDirect :
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
        if funtemp is not None: 
            dod = Model1._thermal_dod
            T_cp, Td = Model1.MSE_ControlPoints(funtemp)
            F2solve = Model1.eval_source_vector(funpowden, indi=dof, indj=dod, Td=Td)
            del dod, Td, T_cp
        else: 
            F2solve = Model1.eval_source_vector(funpowden, indi=dof)
        
        # Solve system
        start = time.time()
        sol_direct = scipy.linalg.solve(K2solve.todense(), F2solve)
        stop = time.time()
        time_direct = stop - start
        time.sleep(1)
        _, memory_direct = tracemalloc.get_traced_memory()
        memory_direct /= 1024*1024
        del K2solve, F2solve, Model1, dof

    # Recursive solver 
    # ----------------
    if doIterative:
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
        if funtemp is not None:  
            dod = Model1._thermal_dod
            T_cp, Td = Model1.MSE_ControlPoints(funtemp)
            F2solve = Model1.eval_source_vector(funpowden, indi=dof, indj=dod, Td=Td)
            del dod, Td, T_cp
        else:
            F2solve = Model1.eval_source_vector(funpowden, indi=dof)

        time.sleep(1)
        _, memory_iter_base = tracemalloc.get_traced_memory()
        memory_iter_base /= 1024*1024

        enablePrint()
        if sol_direct is None: 
            sol_direct = np.ones(len(F2solve))
            print("Direct solution not known. Defalult: ones chosen. Be aware of residue results")
        blockPrint()

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

    # Results 
    # ----------------
    if isDirect is None: 
        output = {"TimeAssembly": time_assembly, "TimeDirect": time_direct, "MemDirect": memory_direct, 
                    "TimeNoIter":time_noiter, "TimeIter": time_iter, "Res": residue, 
                    "Error": error, "MemNoIter": memory_noiter, "MemIter": memory_iter}
    elif isDirect is True: 
        output = {"TimeAssembly": time_assembly, "TimeDirect": time_direct, "MemDirect": memory_direct}
    elif isDirect is False: 
        time_assembly = 1e3
        time_direct = 1e3
        memory_direct = 1e3
        output = {"TimeAssembly": time_assembly, "TimeDirect": time_direct, "MemDirect": memory_direct, 
                    "TimeNoIter":time_noiter, "TimeIter": time_iter, "Res": residue, 
                    "Error": error, "MemNoIter": memory_noiter, "MemIter": memory_iter}
    
    return output

# Some constants
FileExist = False
BlockedBoundaries = [[1, 1], [1, 1], [1, 1]]
DEGREE = 3
CUTS = 3
IS_IGA_GALERKIN = False
GEOMETRY_CASE = 2
if IS_IGA_GALERKIN: is_cg_list = [True]
else: is_cg_list = [False]
isDirect = False

for IS_CG in is_cg_list:   
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
    
    # Get text file name
    txtname += '_p' + str(DEGREE) + '_nbel' + str(2**CUTS)
    if IS_IGA_GALERKIN: txtname += '_IGAG'
    else: txtname += '_IGAWQ'
    if IS_CG: txtname += '_CG'
    else: txtname += '_BiCG'
    txtname = folder + txtname 

    if not FileExist:
        # Run simulation
        blockPrint()
        method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]
        inputs_export = run_simulation(DEGREE, CUTS, GEOMETRY_CASE, funpow, funtemp, IS_IGA_GALERKIN, 
                        method_list, IS_CG, thermalblockedboundaries= BlockedBoundaries, isDirect=isDirect)
        enablePrint()

        # Export results
        write_text_file(txtname+'.txt', method_list, inputs_export)

    else :
        try: 
            # Plot results
            inputs = read_text_file(txtname +'.txt')
            method_list = ["WP", "FD-C", "FD-TDS", "FD-JM", "FD-TD", "FD-JMS"]
            plot_iterative_solver(txtname, inputs, method_list)
        except: 
            print("File does not exist")