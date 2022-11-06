"""
.. Test of assembly and symmetry 
.. We test if the assembly by fortran and python are the same
.. We also test how asymmetric the stiffness and mass matrices are
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_iga import fortran_mf_iga
from lib.fortran_mf_wq import fortran_mf_wq
from lib.python_iga import IGA
from lib.python_wq import WQ
from lib.physics import power_density
from lib.base_functions import relativeError

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test2/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
isIGA       = False
doAssembly  = False
doSymmetry  = True

if isIGA: cfortran = fortran_mf_iga; cpython = IGA
else: cfortran = fortran_mf_wq; cpython = WQ

# ------------------------
# TEST ERROR CONSTRUCTION
# ------------------------
if doAssembly: 
    for geoName in ['VB', 'TR', 'CB']:
        for varName in ['C', 'K', 'F', 'J', 'QP']:

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))

            for degree in range(3, 6):

                norm = []; nbel_list =[]

                for cuts in range(1, 4): 
                    print(degree, cuts, varName)

                    blockPrint()
                    
                    # Define model 
                    geometry = {'degree':[degree, degree, degree]}
                    modelGeo = geomdlModel(geoName, **geometry)
                    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                                np.array([cuts, cuts, cuts]))
                
                    modelPhy1 = cfortran(modelIGA)
                    modelPhy2 = cpython(modelIGA)

                    # Set physical properties
                    material = {'capacity':1, 'conductivity':np.eye(3)}
                    modelPhy1._set_material(material)
                    modelPhy2._set_material(material)
                    funpow = power_density 

                    if varName == 'K': 
                        var1 = modelPhy1.eval_conductivity_matrix()
                        var2 = modelPhy2.eval_conductivity_matrix()

                    elif varName == 'C': 
                        var1 = modelPhy1.eval_capacity_matrix()
                        var2 = modelPhy2.eval_capacity_matrix()
                    
                    elif varName == 'F':
                        var1 = modelPhy1.eval_source_vector(funpow)
                        var2 = modelPhy2.eval_source_vector(funpow)

                    elif varName == 'J':
                        var1 = modelPhy1._detJ
                        var2 = modelPhy2._detJ
                    
                    elif varName == 'QP':
                        var1 = modelPhy1._qp_PS
                        var2 = modelPhy2._qp_PS
                    enablePrint()

                    # Compare results 
                    norm_temp = relativeError(var1, var2)
                    if norm_temp > 1e-5: raise Warning('Fortran and Python give different results')
                    norm.append(norm_temp)
                    nbel_list.append(2**cuts)
                
                label = 'Degree $p =$ ' + str(degree)
                ax.loglog(nbel_list, norm, label=label)

            ax.grid()
            ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
            ax.set_ylabel('Relative error')
            ax.legend()
            fig.tight_layout()
            fig.savefig(folder + 'Error_constructionI_' + geoName + '_' + varName + '.png')

# ------------------------
# TEST ERROR SYMMETRY
# ------------------------
if doSymmetry:
    for geoName in ['VB', 'CB', 'TR', 'RQA']:
        for varName in ['K', 'C']:

            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))

            for degree in range(3, 7):
                
                norm = []; nbel_list =[]
                
                for cuts in range(1, 5): 
                    print(degree, cuts)
                    
                    blockPrint()
                    # Define model 
                    geometry = {'degree':[degree, degree, degree]}
                    modelGeo = geomdlModel(geoName, **geometry)
                    modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                                np.array([cuts, cuts, cuts]))
                
                    modelPhy1 = cfortran(modelIGA)

                    # Set physical properties
                    material = {'capacity':1, 'conductivity':np.eye(3)}
                    modelPhy1._set_material(material)

                    if varName == 'K': var1 = modelPhy1.eval_conductivity_matrix()
                    elif varName == 'C': var1 = modelPhy1.eval_capacity_matrix()
                    enablePrint()

                    # Compare results 
                    norm_temp = relativeError(var1, var1.T)
                    norm.append(norm_temp)
                    nbel_list.append(2**cuts/(degree+1))

                label = 'Degree $p =$ ' + str(degree)
                ax.loglog(nbel_list, norm, label=label)

            ax.set_xlabel('(Parametric support width)' + r'$^{-1}$')
            ax.set_ylabel('Relative error')
            ax.legend()
            fig.tight_layout()
            fig.savefig(folder + 'Error_symmetry_' + geoName + '_' + varName + '.png')
