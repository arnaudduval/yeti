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
withReference = False
degree_list = range(2, 7)
cut_list = range(6, 10)

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

            # Create random vector 
            V = np.random.random(nb_ctrlpts**3)

            # Compute fast diagonalization
            start = time.process_time()
            eig_t, U_t = eigen_decomposition(indi_t, indj_t, data_t)
            eig_diag = solver.compute_diagonal_py(eig_t, eig_t, eig_t, np.arange(4))
            fast_diagonalization(U_t, U_t, U_t, eig_diag, V, fdtype='steady')
            stop = time.process_time()
            FDtime = stop - start
            print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, FDtime))
            timeFD_matrix[i, j+1] = FDtime
            np.savetxt(filename_data, timeFD_matrix)

else:

    # Create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
    markers = ['o', 'v', 's', 'X', '+', 'p']

    # Extract litterature data
    file = pd.read_table(folder_data + 'FD_time.dat', sep=' ', names=['nbel', 'p2', 'p3', 'p4', 'p5', 'p6'])
    nbel = file.nbel
    times = np.column_stack((file.p2, file.p3, file.p4, file.p5, file.p6))

    # Plot data
    for i in range(5):
        nb_ctrlpts = nbel
        ax.loglog(nb_ctrlpts, times[:, i], '--', label='degree ' + r'$p=$' + str(i+2), marker=markers[i])

    # Compute slope
    slope = np.polyfit(np.log10(nb_ctrlpts),np.log10(times[:, 2]), 1)[0]
    slope = round(slope, 3)
    annotation.slope_marker((nb_ctrlpts[1], times[1, 2]), slope, 
                            poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xticks([64, 128, 256, 512])

    # Set properties
    ax.legend(loc='best')
    ax.set_xlabel("Discretization level " + r'$h^{-1}$')
    ax.set_ylabel("CPU time (s)")
    ax.set_ylim([0.01, 100])

    fig.tight_layout()
    fig.savefig(folder_figure + 'FastDiag' + '.pdf')

    if withReference:
        # Create plot
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))

        # Extract litterature data
        file = pd.read_table(folder_data + 'FD_time_lit.dat', sep='\t', names=['nbel', 'p2', 'p3', 'p4', 'p5', 'p6'])
        nbel = file.nbel[:-1]
        times = np.column_stack((file.p2, file.p3, file.p4, file.p5, file.p6))

        # Plot data
        for i in range(np.shape(times)[1]):
            nb_ctrlpts = (nbel + i)**3
            ax.loglog(nb_ctrlpts, times[:-1, i], 'o--', label='degree ' + r'$p=$' + str(i+2))

        # Set properties
        ax.legend(loc='best')
        ax.set_xlabel("Total number of DOF")
        ax.set_ylabel("CPU time (s)")
        
        fig.tight_layout()
        fig.savefig(folder_figure + 'FastDiag_lit' + '.png')