# SIMULATIONS
# ==============================================
# Author : Joaquin CORNEJO

# Python libraries
import os, tracemalloc
import scipy, numpy as np
import time
from datetime import datetime

# My libraries
from pysrc.lib import enablePrint, blockPrint
from pysrc.lib.create_geomdl import create_geometry
from pysrc.lib.fortran_mf_iga import fortran_mf_iga
from pysrc.lib.fortran_mf_wq import fortran_mf_wq
from pysrc.lib.physics import (powden_cube, 
                            powden_prism,
                            powden_thickring, 
                            powden_rotring, 
                            temperature_rotring
)
from pysrc.lib.create_model import write_text_file

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder):
    os.mkdir(folder)

def run_simulation(degree, cuts, geometry_case, funpowden, funtemp, isiga, 
                    method_list, iscg, isOnlyDirect=False, isOnlyIter=False):
    
    # Define actions 
    doDirect, doIterative = True, True
    if isOnlyDirect is True: doIterative = False
    if isOnlyIter is True: 
        doDirect = False
        time_assembly = time_direct = memory_direct = -1e5

    # Define solution 
    sol_direct = None

    # Direct solver
    # -------------
    if doDirect :
        tracemalloc.start()

        # Define geometry 
        modelGeo = create_geometry(degree, cuts, geometry_case)

        # Create thermal model object
        if isiga: Model1 = fortran_mf_iga(modelGeo)
        else: Model1 = fortran_mf_wq(modelGeo)
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
        if isiga: Model1 = fortran_mf_iga(modelGeo)
        else: Model1 = fortran_mf_wq(modelGeo)
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
            print("Direct solution unknown. Default: ones chosen. Be aware of residue results")
        blockPrint()

        # Only compute time to prepare method before iterations
        time_noiter, memory_noiter = [], []
        epsilon  = 1e-14 
        iterations = 0
        tracemalloc.clear_traces()
        for name in method_list:
            start = time.time()
            Model1.mf_conj_grad(F2solve, iterations, epsilon, name, sol_direct, iscg)
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
            _, residue_t, error_t = Model1.mf_conj_grad(F2solve, iterations, epsilon, name, sol_direct, iscg)
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
for CUTS in range(3, 6):
    for IS_IGA_GALERKIN in [False]:
        for GEOMETRY_CASE in ['CB', 'VB', 'TR', 'RQA']:

            if IS_IGA_GALERKIN: is_cg_list = [True]
            else: is_cg_list = [True]
        
            for IS_CG in is_cg_list:
                for DEGREE in range(3, 7):
                    # Print current time
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    print("Current Time =", current_time)
                    
                    # Get file name
                    if GEOMETRY_CASE == 'CB': funpow, funtemp = powden_cube, None 
                    elif GEOMETRY_CASE == 'VB': funpow, funtemp = powden_prism, None 
                    elif GEOMETRY_CASE == 'TR': funpow, funtemp = powden_thickring, None 
                    elif GEOMETRY_CASE == 'RQA': funpow, funtemp = powden_rotring, temperature_rotring 
                    else: raise Warning('Geometry does not exist')
                                    
                    # Get text file name
                    txtname = GEOMETRY_CASE + '_p' + str(DEGREE) + '_nbel' + str(2**CUTS)
                    if IS_IGA_GALERKIN: txtname += '_IGAG'
                    else: txtname += '_IGAWQ'
                    if IS_CG: txtname += '_CG'
                    else: txtname += '_BiCG'

                    print([DEGREE, 2**CUTS, txtname])

                    # Run simulation
                    method_list = ["WP", "C", "TDS", "JM", "TD", "JMS"]
                    blockPrint()
                    inputs_export = run_simulation(DEGREE, CUTS, GEOMETRY_CASE, funpow, funtemp, IS_IGA_GALERKIN, 
                                    method_list, IS_CG, isOnlyIter= True)
                    enablePrint()

                    # Export results
                    txtname = folder + txtname + '.txt'
                    write_text_file(txtname, method_list, inputs_export)

