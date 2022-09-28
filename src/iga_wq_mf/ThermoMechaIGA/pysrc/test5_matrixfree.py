"""
.. Test of fast diagonalization
.. We verify speed of fast diagonalization as direct method to solve Sylvester system
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.base_functions import (create_knotvector, 
                                wq_find_basis_weights_fortran, 
                                fast_diagonalization, 
                                erase_rows_csr, 
                                eigen_decomposition, 
                                compute_eig_diag
)

# Select folder
full_path = os.path.realpath(__file__)
folder_data = os.path.dirname(full_path) + '/data/'
folder_figure = os.path.dirname(full_path) + '/results/test5/'
if not os.path.isdir(folder_figure): os.mkdir(folder_figure)

# Set filename
filename_data = folder_data + 'FD_time.dat' 

# Set global variable
dataExist = True
withReference = True
degree_list = range(2, 7)
cut_list = range(5, 10)

# Initialize
timeFD_matrix = np.zeros((len(cut_list), len(degree_list)+1))
timeFD_matrix[:, 0] = np.array([2**i for i in cut_list])


if not dataExist:

    for i, cuts in enumerate(cut_list):
        for j, degree in enumerate(degree_list):
        
            # Set number of elements
            nbel = 2**cuts
            nb_ctrlpts = degree + nbel - 2

            # Basis and weights
            knotvector = create_knotvector(degree, nbel)
            B, W, indi, indj = wq_find_basis_weights_fortran(degree, knotvector)[2:]

            # Erase data
            rows2erase = [0, -1]
            indi_t, indj_t, [B_t, W_t] = erase_rows_csr(rows2erase, indi, indj, [B, W])
            data_t = [B_t[:, 0], B_t[:, 1], W_t[:, 0], W_t[:, -1]]

            # Compute eigen decomposition
            eig_t, U_t = eigen_decomposition(indi_t, indj_t, data_t)
            eig_diag = compute_eig_diag(eig_t, eig_t, eig_t)

            # Create random vector 
            V = np.random.random(nb_ctrlpts**3)

            # Compute fast diagonalization
            start = time.time()
            fast_diagonalization(U_t, U_t, U_t, eig_diag, V, fdtype='steady')
            stop = time.time()
            FDtime = stop - start
            print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, FDtime))
            timeFD_matrix[i, j+1] = FDtime
            # np.savetxt(filename_data, timeFD_matrix)

else:

    # Create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))

    # Extract litterature data
    file = pd.read_table(folder_data + 'matvec_MF.dat', sep='\t', names=['degree', 'Ku64', 'Cu64']) 
    degree = file.degree
    arrays = [file.Cu64, file.Ku64]
    labels = ['MF-WQ Mass', 'MF-WQ Stiffness']

    for array, label in zip(arrays, labels):
        ax.semilogy(degree, array, 'o--', label=label)

    # Set properties
    ax.legend(loc='best')
    ax.set_xlabel("Degree")
    ax.set_ylabel("CPU time (s)")
    ax.set_xlim([1, 10])
    ax.set_ylim([0.001, 10])
    fig.tight_layout()
    fig.savefig(folder_figure + 'ProductMF_lit' + '.png')


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
        ax.set_xlim([1, 10])
        ax.set_ylim([0.001, 10])
        fig.tight_layout()
        fig.savefig(folder_figure + 'ProductMF_lit' + '.png')