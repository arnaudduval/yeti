"""
.. Test of fast diagonalization
.. We verify speed of fast diagonalization as direct method to solve Sylvester system
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.create_geomdl import geomdlModel
from lib.fortran_mf_wq import fortran_mf_wq

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test5/'
if not os.path.isdir(folder): os.mkdir(folder)
folder_data = os.path.dirname(full_path) + '/data/'

# Set global variables
dataExist     = True
withReference = False
degree_list   = range(2, 10)
cuts          = 6
extension     = '.pdf'

if not dataExist:

    timeMF_Mass_matrix = np.zeros((len(degree_list), 2))
    timeMF_Mass_matrix[:, 0] = degree_list
    timeMF_Stiff_matrix = np.zeros((len(degree_list), 2))
    timeMF_Stiff_matrix[:, 0] = degree_list
    timeMF_SM_matrix = np.zeros((len(degree_list), 2))
    timeMF_SM_matrix[:, 0] = degree_list
    timePython = np.zeros((len(degree_list), 2))
    timePython[:, 0] = degree_list

    for i, degree in enumerate(degree_list):
        
        blockPrint()

        nbel = 2**cuts
        nb_ctrlpts = degree + nbel - 2

        # Create model
        geometry = {'degree':[degree, degree, degree]}
        modelGeo = geomdlModel('CB', **geometry)
        modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                    np.array([cuts, cuts, cuts]))
        modelPhy = fortran_mf_wq(modelIGA)

        # Add material 
        conductivity = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
        material = {'capacity':1.0, 'conductivity':conductivity}
        modelPhy._set_material(material)

        # Block boundaries
        Dirichlet = np.array([[1, 1], [1, 1], [1, 1]])
        Dirichlet_dict = {'thermal':Dirichlet}
        modelPhy._set_dirichlet_boundaries(Dirichlet_dict)

        V = np.random.random(nb_ctrlpts**3)

        # --------------
        # Compute matrix
        # --------------
        dof = modelPhy._thermal_dof
        CC = modelPhy.eval_capacity_matrix()[:, dof][dof, :]
        
        start = time.process_time()
        R = CC.dot(V)
        stop = time.process_time()

        time_python = stop - start
        timePython[i, 1] = time_python
        np.savetxt(folder_data+'matvec_Python_'+str(cuts)+'.dat', timePython)

        # ------------------
        # Compute MF product
        # ------------------
        MFtime1 = modelPhy.eval_Cu(V, table=Dirichlet)[-1]
        timeMF_Mass_matrix[i, 1] = MFtime1
        np.savetxt(folder_data+'matvec_MF_Mass_'+str(cuts)+'.dat', timeMF_Mass_matrix)

        MFtime2 = modelPhy.eval_Ku(V, table=Dirichlet)[-1]
        timeMF_Stiff_matrix[i, 1] = MFtime2
        np.savetxt(folder_data+'matvec_MF_Stiff_'+str(cuts)+'.dat', timeMF_Stiff_matrix)

        MFtime3 = modelPhy.eval_KCu(V, table=Dirichlet)[-1]
        timeMF_SM_matrix[i, 1] = MFtime3
        np.savetxt(folder_data+'matvec_MF_SM_'+str(cuts)+'.dat', timeMF_SM_matrix)
        
        enablePrint()
        print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, time_python))
        print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, MFtime1))
        print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, MFtime2))
        print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, MFtime3))

else:

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))
    markers = ['o', 'v', 's', 'X', '+', 'p']

    # Load data
    file_P = pd.read_table(folder_data + 'matvec_Python_'+str(cuts)+'.dat', sep=' ', names=['degree', 'P1', 'P2']) 
    file_M = pd.read_table(folder_data + 'matvec_MF_Mass_'+str(cuts)+'.dat', sep=' ', names=['degree', 'P1', 'P2']) 
    file_K = pd.read_table(folder_data + 'matvec_MF_Stiff_'+str(cuts)+'.dat', sep=' ', names=['degree', 'P1', 'P2']) 
    file_KM = pd.read_table(folder_data + 'matvec_MF_SM_'+str(cuts)+'.dat', sep=' ', names=['degree', 'P1', 'P2']) 
    degree = file_M.degree
    arrays = [file_M.P1, file_K.P1, file_KM.P1]
    labels = [  r'$\mathsf{M}x$' + ', ' + r'$h^{-1}=$ ' + str(2**cuts), 
                r'$\mathsf{K}x$' + ', ' + r'$h^{-1}=$ ' + str(2**cuts), 
                r'$\mathsf{A}x$' + ', ' + r'$h^{-1}=$ ' + str(2**cuts)]

    ax.semilogy(file_P.degree, file_P.P1, '--', label='Python'+', '+r'$h^{-1}=$ ' + str(2**cuts), marker=markers[0])
    for i, [array, label] in enumerate(zip(arrays, labels)):
        ax.semilogy(degree, array, '--', label=label, marker=markers[i+1])

    ax.legend(loc='best')
    ax.set_xlabel('Degree ' + r'$p$')
    ax.set_ylabel('CPU time (s)')
    ax.set_xlim([1, 11])
    ax.set_ylim([0.01, 100])

    fig.tight_layout()
    fig.savefig(folder + 'ProductMF' + extension)

    if withReference:

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))

        # Load data
        file = pd.read_table(folder_data + 'matvec_MF_lit.dat', sep='\t', names=['degree', 'Ku64', 'Cu64']) 
        arrays = [file.Cu64, file.Ku64]
        labels = ['MF-WQ Mass', 'MF-WQ Stiffness']

        for array, label in zip(arrays, labels):
            ax.semilogy(file.degree, array, 'o--', label=label)

        # Set properties
        ax.legend(loc='best')
        ax.set_xlabel('Degree')
        ax.set_ylabel('CPU time (s)')
        ax.set_xlim([1, 11])
        ax.set_ylim([0.001, 10])
        fig.tight_layout()
        fig.savefig(folder + 'ProductMF_lit' + '.png')