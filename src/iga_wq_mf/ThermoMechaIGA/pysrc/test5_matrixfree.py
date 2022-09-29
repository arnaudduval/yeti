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
folder_data = os.path.dirname(full_path) + '/data/'
folder_figure = os.path.dirname(full_path) + '/results/test5/'
if not os.path.isdir(folder_figure): os.mkdir(folder_figure)

# Set global variable
dataExist = True
withReference = True
degree_list = range(7, 8)
cut_list = range(5, 6)

# Initialize
timeMF_Mass_matrix = np.zeros((len(degree_list), len(cut_list)+1))
timeMF_Mass_matrix[:, 0] = degree_list
timeMF_Stiff_matrix = np.zeros((len(degree_list), len(cut_list)+1))
timeMF_Stiff_matrix[:, 0] = degree_list
timePython = np.zeros((len(degree_list), len(cut_list)+1))
timePython[:, 0] = degree_list

if not dataExist:

    for j, cuts in enumerate(cut_list):
        for i, degree in enumerate(degree_list):
            
            # blockPrint()
            # Set number of elements
            nbel = 2**cuts
            nb_ctrlpts = degree + nbel - 2

            # Create geometry
            geometry = {'degree':[degree, degree, degree]}
            modelGeo = geomdlModel('CB', **geometry)
            modelIGA = modelGeo.export_IGAparametrization(nb_refinementByDirection=
                                                        np.array([cuts, cuts, cuts]))
            # Create model
            modelPhy = fortran_mf_wq(modelIGA)

            # Add material 
            conductivity = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
            material = {'capacity':1.0, 'conductivity':conductivity}
            modelPhy._set_material(material)

            # Block boundaries
            Dirichlet = np.array([[1, 1], [1, 1], [1, 1]])
            Dirichlet_dict = {'thermal':Dirichlet}
            modelPhy._set_dirichlet_boundaries(Dirichlet_dict)

            # Create random vector 
            V = np.random.random(nb_ctrlpts**3)

            # # Compute matrix C
            # dof = modelPhy._thermal_dof
            # CC = modelPhy.eval_capacity_matrix()[:, dof][dof, :]
            # start = time.time()
            # R = CC.dot(V)
            # stop = time.time()
            # time_python = stop - start
            # timePython[i, j+1] = time_python
            # # np.savetxt(folder_data+'matvec_Python.dat', timePython)

            # Compute Matrix free product
            MFtime1 = modelPhy.eval_Cu(V, table=Dirichlet)[-1]
            timeMF_Mass_matrix[i, j+1] = MFtime1
            np.savetxt(folder_data+'matvec_MF_Mass.dat', timeMF_Mass_matrix)

            MFtime2 = modelPhy.eval_Ku(V, table=Dirichlet)[-1]
            timeMF_Stiff_matrix[i, j+1] = MFtime2
            np.savetxt(folder_data+'matvec_MF_Stiff.dat', timeMF_Stiff_matrix)
            
            enablePrint()
            # print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, time_python))
            print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, MFtime1))
            print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, MFtime2))
            print('')

else:

    # Create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))

    # Extract litterature data
    file_P = pd.read_table(folder_data + 'matvec_Python.dat', sep=' ', names=['degree', 'Cu32', 'Cu64']) 
    file_M = pd.read_table(folder_data + 'matvec_MF_Mass.dat', sep=' ', names=['degree', 'Cu32', 'Cu64']) 
    file_K = pd.read_table(folder_data + 'matvec_MF_Stiff.dat', sep=' ', names=['degree', 'Ku32', 'Ku64']) 
    degree = file_M.degree
    arrays = [file_M.Cu64, file_K.Ku64]
    labels = ['MF-WQ Mass', 'MF-WQ Stiffness']

    for array, label in zip(arrays, labels):
        ax.semilogy(degree, array, 'o--', label=label)

    ax.semilogy(file_P.degree, file_P.Cu64, 'o--', label='Python built-in SpMDV')

    # Set properties
    ax.legend(loc='best')
    ax.set_xlabel("Polynomial degree")
    ax.set_ylabel("CPU time (s)")
    ax.set_xlim([1, 11])
    ax.set_ylim([0.01, 100])
    fig.tight_layout()
    fig.savefig(folder_figure + 'ProductMF' + '.pdf')

    if withReference:
        # Create plot
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))

        # Extract litterature data
        file = pd.read_table(folder_data + 'matvec_MF_lit.dat', sep='\t', names=['degree', 'Ku64', 'Cu64']) 
        degree = file.degree
        arrays = [file.Cu64, file.Ku64]
        labels = ['MF-WQ Mass', 'MF-WQ Stiffness']

        for array, label in zip(arrays, labels):
            ax.semilogy(degree, array, 'o--', label=label)

        # Set properties
        ax.legend(loc='best')
        ax.set_xlabel("Degree")
        ax.set_ylabel("CPU time (s)")
        ax.set_xlim([1, 11])
        ax.set_ylim([0.001, 10])
        fig.tight_layout()
        fig.savefig(folder_figure + 'ProductMF_lit' + '.png')