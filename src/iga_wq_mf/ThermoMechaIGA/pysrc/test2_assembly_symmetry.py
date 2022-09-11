"""
.. Test of assembly and symmetry 
.. We test if the assembly by fortran and python are the same
.. We also test how asymetric are K and C matrices
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_iga import fortran_mf_iga
from lib.fortran_mf_wq import fortran_mf_wq
from lib.python_iga import IGA
from lib.python_wq import WQ
from lib.physics import power_density

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test2/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
isIGA = False
doConstruction = True
doSymetry = False

if isIGA: cfortran = fortran_mf_iga; cpython = IGA
else: cfortran = fortran_mf_wq; cpython = WQ

# ------------------------
# TEST ERROR CONSTRUCTION
# ------------------------
if doConstruction: 
    for geoName in ['VB', 'TR', 'CB']:
        for varName in ['K', 'C', 'F', 'J', 'QP']:

            # Create plot
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,4))

            for degree in range(3, 6):

                # Initialize
                norm = []; nbel_list =[]

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
                    modelPhy1 = cfortran(modelIGA)
                    modelPhy2 = cpython(modelIGA)

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
                    try: norm_temp = sp.linalg.norm(error, np.inf)/sp.linalg.norm(var1, np.inf)*100
                    except: norm_temp = np.linalg.norm(error, np.inf)/np.linalg.norm(var1, np.inf)*100
                    if norm_temp > 1e-5: raise Warning("Fortran and Python give different results")
                    norm.append(norm_temp)

                    # Set number of elements
                    nbel_list.append(2**cuts)
                
                # Plot 
                strlabel = 'Degree p = ' + str(degree)
                ax.loglog(nbel_list, norm, label=strlabel)

            # Properties
            ax.grid()
            ax.set_xlabel("Number of elements $nb_{el}$")
            ax.set_ylabel("Relative error (%)")
            ax.legend()
            fig.tight_layout()
            fig.savefig(folder + 'Error_constructionI_' + geoName + '_' + varName + '.png')


# ------------------------
# TEST ERROR SYMMETRY
# ------------------------
if doSymetry:
    for geoName in ['VB', 'CB', 'TR', 'RQA']:
        for varName in ['K', 'C']:

            # Create plot
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,4))


            for degree in range(3, 6):
                
                # Initialize
                norm = []; nbel_list =[]
                
                for cuts in range(1, 5): 
                    print(degree, cuts)
                    
                    # blockPrint()
                    # Define geometry 
                    geometry = {'degree':[degree, degree, degree]}
                    modelGeo = geomdlModel(geoName, **geometry)
                    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                                np.array([cuts, cuts, cuts]))
                
                    # Creation of thermal model object
                    modelPhy1 = cfortran(modelIGA)

                    # Set physical properties
                    material = {'capacity':1, 'conductivity':np.eye(3)}
                    modelPhy1._set_material(material)

                    if varName == "K": var1 = modelPhy1.eval_conductivity_matrix()
                    elif varName == "C": var1 = modelPhy1.eval_capacity_matrix()
                    # enablePrint()

                    # Compare results 
                    error = var1.T - var1
                    norm_temp = sp.linalg.norm(error, np.inf)/sp.linalg.norm(var1, np.inf)*100
                    norm.append(norm_temp)

                    # Set number of elements
                    nbel_list.append(2**cuts/(degree+1))

                # Plot 
                strlabel = 'Degree p = ' + str(degree)
                ax.loglog(nbel_list, norm, label=strlabel)

            # Properties
            ax.set_xlabel("(Parametric support width)" + r"$^{-1}$", fontsize= 16)
            ax.set_ylabel("Relative error (%)", fontsize= 16)
            ax.legend()
            fig.tight_layout()
            fig.savefig(folder + 'Error_symmetry_' + geoName + '_' + varName + '.png')
