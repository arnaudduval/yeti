"""
.. Test of assembly and symmetry 
.. We test if the assembly by fortran and python are the same
.. We also test how asymetric are K and C matrices
.. Joaquin Cornejo 
"""

# Python libraries
import os
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# My libraries
from lib.__init__ import blockPrint, enablePrint
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq
from lib.python_wq import WQ
from lib.fortran_mf_iga import fortran_mf_iga
from lib.python_iga import IGA
from lib.physics import power_density

# Choose folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test2/'
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

# =========================
# TEST ERROR CONSTRUCTION
# =========================
if CONSTRUCTION: 
    # Set degree and number of divisions
    for geoName in ['VB', 'TR', 'RQA', 'CB']:
        for varName in ['C', 'K', 'F', 'J', 'QP']:
            for degree in range(3, 6):
                norm = []; ddl =[]
                for cuts in range(1, 4): 
                    print(degree, cuts, varName)

                    blockPrint()
                    # Get file name
                    funpow = power_density 

                    # Define geometry 
                    geometry = {'degree':[degree, degree, degree]}
                    modelGeo = geomdlModel(geoName, **geometry)
                    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                                np.array([cuts, cuts, cuts]))
                
                    # Creation of thermal model object
                    modelPhy1 = classfortran(modelIGA)
                    modelPhy2 = classpython(modelIGA)

                    # Set physical properties
                    material = {'capacity':1, 'conductivity':np.eye(3)}
                    modelPhy1._set_material(material)
                    modelPhy2._set_material(material)

                    if varName == "K": 
                        var1 = modelPhy1.eval_conductivity_matrix()
                        var2 = modelPhy2.eval_conductivity_matrix()

                    elif varName == "C": 
                        var1 = modelPhy1.eval_capacity_matrix()
                        var2 = modelPhy2.eval_capacity_matrix()
                    
                    elif varName == 'F':
                        var1 = modelPhy1.eval_source_vector(funpow)
                        var2 = modelPhy2.eval_source_vector(funpow)

                    elif varName == 'J':
                        var1 = modelPhy1._detJ
                        var2 = modelPhy2._detJ
                    
                    elif varName == 'QP':
                        var1 = modelPhy1._qp_PS[:,:]
                        var2 = modelPhy2._qp_PS[:,:]
                
                    enablePrint()

                    # Compare results 
                    error = var1 - var2
                    try: norm_temp = sparse.linalg.norm(error, np.inf)/sparse.linalg.norm(var1, np.inf)
                    except: norm_temp = np.linalg.norm(error, np.inf)/np.linalg.norm(var1, np.inf)
                    if norm_temp > 1e-5:
                        raise Warning("Something happend. Fortran and Python give different results")
                    norm.append(norm_temp)

                    # Set number of elements
                    nbel = 2 ** cuts
                    ddl.append(nbel)

                # Change type 
                norm = np.asarray(norm)
                ddl = np.asarray(ddl)

                # Figure 
                plt.figure(1)
                plt.plot(ddl, norm*100, label='degree p = ' + str(degree))

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
            plt.savefig(folder + 'Error_constructionI_' + geoName + '_' + varName + '.png')
            plt.figure(1).clear()
        print('----')

# =========================
# TEST ERROR SYMMETRY
# =========================
if SYMMETRY:
    # Set degree and number of divisions
    for geoName in ['CB', 'VB', 'TR', 'RQA']:
        for varName in ['K', 'C']:
            for degree in range(3, 6):
                norm = []; ddl =[]
                for cuts in range(1, 5): 
                    print(degree, cuts)
                    
                    blockPrint()

                    # Define geometry 
                    geometry = {'degree', [degree, degree, degree]}
                    modelGeo = geomdlModel(geoName, **geometry)
                    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                                np.array([cuts, cuts, cuts]))

                    # Creation of thermal model object
                    modelPhy1 = classfortran(modelIGA)
                    del modelGeo, modelIGA

                    if varName == "K": var1 = modelPhy1.eval_conductivity_matrix()
                    elif varName == "C": var1 = modelPhy1.eval_capacity_matrix()
                    del modelPhy1
                    enablePrint()

                    # Compare results 
                    error = var1.transpose() - var1
                    norm_temp = sparse.linalg.norm(error, np.inf)/sparse.linalg.norm(var1, np.inf)
                    norm.append(norm_temp)

                    # Set number of elements
                    nbel = 2 ** cuts
                    ddl.append(nbel/(degree+1))

                # Change type 
                norm = np.asarray(norm)
                ddl = np.asarray(ddl)

                # Figure 
                plt.figure(1)
                plt.plot(ddl, norm*100, label='p = ' + str(degree))

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
            plt.savefig(folder + 'Error_symmetry_' + geoName + '_' + varName + '.png')
            plt.figure(1).clear()
        print('----')
