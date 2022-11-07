"""
.. Test of setup time
.. We test how much time it takes to compute K and C matrices
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.base_functions import wq_find_basis_weights_fortran, create_knotvector

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test3/'
if not os.path.isdir(folder): os.mkdir(folder)
folder_data = os.path.dirname(full_path) + '/data/'

# Set global variables
dataExist   = True
nbel        = 64 # or 40
degree_list = np.arange(2, 8)

# Set filename
filename_data = folder_data + 'setup_time_' + str(nbel) + '.dat'
filename_figure = folder + 'setup_time_' + str(nbel) + '.pdf' 

if not dataExist:

    time_matrix = np.zeros((len(degree_list), 3))
    time_matrix[:, 0] = degree_list
    
    for i, degree in enumerate(degree_list):

        knotvector = create_knotvector(degree, nbel)
        nnz_I, qp_position, B, W, indi, indj = wq_find_basis_weights_fortran(degree, knotvector)
        nb_qp_list, indices, dB, dW, size_I = [], [], [], [], []
        for dim in range(3):
            nb_qp_list.append(len(qp_position))
            indices.append(indi); indices.append(indj) 
            dB.append(B); dW.append(W); size_I.append(nnz_I)

        # ----------------
        # Capacity matrix
        # ----------------
        coefs = np.ones(len(qp_position)**3)
        inputs = [coefs, *nb_qp_list, *indices, *dB, *dW, *size_I]
        
        start = time.process_time()
        assembly.wq_get_capacity_3d(*inputs)
        stop = time.process_time()
        time_matrix[i, 1] = stop - start         

        # --------------------
        # Conductivity matrix
        # --------------------
        coefs = np.ones((3, 3, len(qp_position)**3))
        inputs = [coefs, *nb_qp_list, *indices, *dB, *dW, *size_I]

        start = time.process_time()
        assembly.wq_get_conductivity_3d(*inputs)
        stop = time.process_time()
        time_matrix[i, 2] = stop - start 
        del inputs, coefs

        print(time_matrix[i, :])
        np.savetxt(filename_data, time_matrix)

else :

    if nbel == 40:
        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))

        # Load data
        mydata = np.loadtxt(filename_data)
        litterature = pd.read_table(folder_data + 'setup_time_lit_40.dat', sep='\t', names=['degree', 'time'])

        ax1.loglog(litterature.degree, litterature.time, 'o--', label='Literature')
        ax1.loglog(mydata[:,0], mydata[:,1], 'x--', label='This work')

        # Compute slope
        slope = np.polyfit(np.log10(mydata[:,0]), np.log10(mydata[:,1]), 1)[0]
        slope = round(slope, 3)
        annotation.slope_marker((mydata[2,0], mydata[2,1]), slope, 
                                poly_kwargs={'facecolor': (0.73, 0.8, 1)})

        # Reference
        c1 = 0.125
        xx = litterature.degree
        yy = c1*xx**3
        ax1.loglog(xx, yy, '-', label=r'$O(p^3)$')

        ax1.legend()
        ax1.set_xlim([2, 10])
        ax1.set_ylim([1e-1, 1e3])
        ax1.set_ylabel('CPU setup time (s)')
        ax1.set_xlabel('Polynomial degree $p$')

        fig.tight_layout()
        fig.savefig(filename_figure)

    if nbel == 64:
        fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
        
        # Load data
        mydata = np.loadtxt(filename_data)
        litterature = pd.read_table(folder_data + 'setup_time_lit_64.dat', sep='\t', names=['degree', 'C', 'K'])

        ax1.loglog(litterature.degree, litterature.C, 'o--', label='Literature')
        ax1.loglog(mydata[:,0], mydata[:,1], 'x--', label='This work')

        ax2.loglog(litterature.degree, litterature.K, 'o--')
        ax2.loglog(mydata[:,0], mydata[:,2], 'x--')

        # Reference
        c1, c2 = 0.6, 0.125
        xx = litterature.degree
        yy = c1*xx**3
        yy2 = c2*xx**3
        ax1.loglog(xx, yy, '-', label=r'$O(p^3)$')
        ax2.loglog(xx, yy2, '-')
        ax1.legend()

        for ax in [ax1, ax2]:
            ax.set_xlim([2, 10])
            ax.set_ylim([1, 1e4])
            ax.set_ylabel('CPU setup time (s)')
            ax.set_xlabel('Polynomial degree $p$')

        fig.tight_layout()
        fig.savefig(filename_figure)

        


