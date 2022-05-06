"""
.. Test of assembly and symmetry 
.. We test how asymetric are K and C matrices
.. Joaquin Cornejo 
"""

# Python libraries
import sys, os
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# My libraries
from lib.geomdl_geometry import create_geometry
from lib.fortran_mf_wq import fortran_mf_wq
from lib.methods_wq import WQ
from lib.fortran_mf_iga import fortran_mf_iga
from lib.methods_iga import IGA
from lib.physics import (powden_cube, 
                        powden_prism,
                        powden_thickring
)

# Enable and disable print
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'

CONSTRUCTION = True
SYMMETRY = False

# ====================================================================
# TEST ERROR CONSTRUCTION
# ====================================================================
if CONSTRUCTION: 
    # Set degree and number of divisions
    for varName in ['K', 'C', 'F']:
        for GEOMETRY_CASE in range(1, 3):
            for DEGREE in range(3, 5):
                norm = []; ddl =[]
                for CUTS in range(1, 4): 
                    print(DEGREE, CUTS)

                    blockPrint()
                    # Get file name
                    if GEOMETRY_CASE == 0: 
                        txtname = 'CB' 
                        funpow = powden_cube 
                    elif GEOMETRY_CASE == 1: 
                        txtname = 'VB' 
                        funpow = powden_prism 
                    elif GEOMETRY_CASE == 2: 
                        txtname = 'TR' 
                        funpow = powden_thickring 

                    # Define geometry 
                    modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)
                
                    # Creation of thermal model object
                    # Model1 = fortran_mf_wq(modelGeo)
                    # Model2 = WQ(modelGeo)

                    Model1 = fortran_mf_iga(modelGeo)
                    Model2 = IGA(modelGeo)

                    if varName == "K": 
                        var1 = Model1.eval_conductivity_matrix()
                        var2 = Model2.eval_conductivity_matrix()

                    elif varName == "C": 
                        var1 = Model1.eval_capacity_matrix()
                        var2 = Model2.eval_capacity_matrix()
                    
                    elif varName == 'F':
                        var1 = Model1.eval_source_vector(funpow)
                        var2 = Model2.eval_source_vector(funpow)

                    elif varName == 'J':
                        var1 = Model1._detJ
                        var2 = Model2._detJ

                    elif varName == 'QP': 
                        var1 = Model1._qp_PS[0, :, :]
                        var2 = Model2._qp_PS[0, :, :]
                
                    enablePrint()

                    # Compare results 
                    error = var1 - var2
                    try: norm_temp = sparse.linalg.norm(error, np.inf)/sparse.linalg.norm(var1, np.inf)
                    except: norm_temp = np.linalg.norm(error, np.inf)/np.linalg.norm(var1, np.inf)
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
            plt.savefig(folder + 'Error_constructionI_' + txtname + '_' + varName + '.png')
            plt.figure(1).clear()

# ====================================================================
# TEST ERROR SYMMETRY
# ====================================================================
if SYMMETRY:
    # Set degree and number of divisions
    for varName in ['K', 'C']:
        for GEOMETRY_CASE in range(3):
            for DEGREE in range(3, 6):
                norm = []; ddl =[]
                for CUTS in range(1, 5): 
                    print(DEGREE, CUTS)
                    
                    blockPrint()
                    # Get file name
                    if GEOMETRY_CASE == 0: txtname = 'CB' 
                    elif GEOMETRY_CASE == 1: txtname = 'VB' 
                    elif GEOMETRY_CASE == 2: txtname = 'TR' 

                    # Define geometry 
                    modelGeo = create_geometry(DEGREE, CUTS, GEOMETRY_CASE)

                    # Creation of thermal model object
                    Model1 = fortran_mf_wq(modelGeo)
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
            plt.savefig(folder + 'Error_symmetry_' + txtname + '_' + varName + '.png')
            plt.figure(1).clear()
