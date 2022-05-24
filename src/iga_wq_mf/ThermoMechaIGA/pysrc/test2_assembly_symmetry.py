"""
.. Test of assembly and symmetry 
.. We test if the assembly by fortran and python are the same
.. We also test how asymetric are K and C matrices
.. Joaquin Cornejo 
"""

# Python libraries
import sys, os
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# My libraries
from lib import enablePrint, blockPrint
from lib.create_geomdl import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq
from lib.methods_wq import WQ
from lib.fortran_mf_iga import fortran_mf_iga
from lib.methods_iga import IGA
from lib.physics import power_density

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

isIGA = False
CONSTRUCTION = True
SYMMETRY = False

if isIGA: 
    classfortran = fortran_mf_iga
    classpython = IGA
else: 
    classfortran = fortran_mf_wq
    classpython = WQ

# ====================================================================
# TEST ERROR CONSTRUCTION
# ====================================================================
if CONSTRUCTION: 
    # Set degree and number of divisions
    for varName in ['K', 'C', 'F']:
        for GEOMETRY_CASE in ['CB', 'VB', 'TR', 'RQA']:
            for DEGREE in range(3, 6):
                norm = []; ddl =[]
                for CUTS in range(1, 4): 
                    print(DEGREE, CUTS)

                    blockPrint()
                    # Get file name
                    if GEOMETRY_CASE == 'CB': funpow = power_density 
                    elif GEOMETRY_CASE == 'VB': funpow = power_density 
                    elif GEOMETRY_CASE == 'TR': funpow = power_density 
                    elif GEOMETRY_CASE == 'RQA': funpow = power_density

                    # Define geometry 
                    modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)
                
                    # Creation of thermal model object
                    Model1 = classfortran(modelGeo)
                    Model2 = classpython(modelGeo)

                    if varName == "K": 
                        var1 = Model1.eval_conductivity_matrix()
                        var2 = Model2.eval_conductivity_matrix()

                    elif varName == "C": 
                        var1 = Model1.eval_capacity_matrix()
                        var2 = Model2.eval_capacity_matrix()
                    
                    elif varName == 'F':
                        var1 = Model1.eval_source_vector(funpow)
                        var2 = Model2.eval_source_vector(funpow)
                
                    enablePrint()

                    # Compare results 
                    error = var1 - var2
                    try: norm_temp = sparse.linalg.norm(error, np.inf)/sparse.linalg.norm(var1, np.inf)
                    except: norm_temp = np.linalg.norm(error, np.inf)/np.linalg.norm(var1, np.inf)
                    if norm_temp > 1e-5:
                        raise Warning("Something happend. Fortran and Python give different results")
                    norm.append(norm_temp)

                    # Set number of elements
                    nbel = 2 ** CUTS
                    ddl.append(nbel)

                # Change type 
                norm = np.asarray(norm)
                ddl = np.asarray(ddl)

                # Figure 
                plt.figure(1)
                plt.plot(ddl, norm*100, label='degree p = ' + str(DEGREE))

            # Properties
            plt.grid()
            plt.xscale("log")
            plt.xlabel("Number of elements $nb_{el}$", fontsize= 16)
            plt.yscale("log")
            plt.ylabel("Relative error (%)", fontsize= 16)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.xlim(1, 100)
            plt.legend()
            plt.tight_layout()
            plt.savefig(folder + 'Error_constructionI_' + GEOMETRY_CASE + '_' + varName + '.png')
            plt.figure(1).clear()

# ====================================================================
# TEST ERROR SYMMETRY
# ====================================================================
if SYMMETRY:
    # Set degree and number of divisions
    for varName in ['K', 'C']:
        for GEOMETRY_CASE in ['CB', 'VB', 'TR', 'RQA']:
            for DEGREE in range(3, 6):
                norm = []; ddl =[]
                for CUTS in range(1, 5): 
                    print(DEGREE, CUTS)
                    
                    blockPrint()

                    # Define geometry 
                    modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)

                    # Creation of thermal model object
                    Model1 = classfortran(modelGeo)
                    del modelGeo

                    if varName == "K": var1 = Model1.eval_conductivity_matrix()
                    elif varName == "C": var1 = Model1.eval_capacity_matrix()
                    del Model1
                    enablePrint()

                    # Compare results 
                    error = var1.transpose() - var1
                    norm_temp = sparse.linalg.norm(error, np.inf)/sparse.linalg.norm(var1, np.inf)
                    norm.append(norm_temp)

                    # Set number of elements
                    nbel = 2 ** CUTS
                    ddl.append(nbel/(DEGREE+1))

                # Change type 
                norm = np.asarray(norm)
                ddl = np.asarray(ddl)

                # Figure 
                plt.figure(1)
                plt.plot(ddl, norm*100, label='p = ' + str(DEGREE))

            # Properties
            plt.grid()
            plt.xscale("log")
            plt.xlabel("(Parametric support width)" + r"$^{-1}$", fontsize= 16)
            plt.yscale("log")
            plt.ylabel("Relative error (%)", fontsize= 16)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.xlim(0.1, 100)
            plt.legend()
            plt.tight_layout()
            plt.savefig(folder + 'Error_symmetry_' + GEOMETRY_CASE + '_' + varName + '.png')
            plt.figure(1).clear()
