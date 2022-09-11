"""
.. Test of setup time
.. We test how much time it takes to compute K and C matrices
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.base_functions import wq_find_basis_weights_fortran, create_knotvector
from iga_wq_mf import assembly

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test3/'
if not os.path.isdir(folder): os.mkdir(folder)

# Simulation setup
dataExist = False
nbel = 64
degree_list = np.arange(2, 11)

# Set filename
filename_WQ = folder + 'setup_time_WQ_'  + 'nbel_' + str(nbel) 

# Output
time_matrix = np.zeros((len(degree_list), 3))
time_matrix[:, 0] = degree_list

if not dataExist:
    for i, degree in enumerate(degree_list):

        # Define basis (the same for all dimensions)
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
        
        start = time.time()
        assembly.wq_get_capacity_3d(*inputs)
        stop = time.time()
        time_matrix[i, 1] = stop - start 
        del inputs, coefs

        # --------------------
        # Conductivity matrix
        # --------------------
        coefs = np.ones((3, 3, len(qp_position)**3))
        inputs = [coefs, *nb_qp_list, *indices, *dB, *dW, *size_I]

        start = time.time()
        assembly.wq_get_conductivity_3d(*inputs)
        stop = time.time()
        time_matrix[i, 2] = stop - start 
        del inputs, coefs

        np.savetxt(filename_WQ, time_matrix)

else :

    if nbel == 40:
        # Load data
        mydata = np.loadtxt(filename_WQ)
        litterature = pd.read_table(folder + 'lit_WQ.dat', sep='\t', names=['degree', 'time'])

        # Plot data
        fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
        ax1.loglog(litterature.degree, litterature.time, 'o--', label='Litterature')
        ax1.loglog(mydata[:,0], mydata[:,1], 'x--', label='My algorithm')

        # Compute slope
        slope = np.polyfit(np.log10(mydata[:,0]), np.log10(mydata[:,1]), 1)[0]
        slope = round(slope, 3)
        annotation.slope_marker((mydata[2,0], mydata[2,1]), slope, 
                                text_kwargs={'fontsize': 14},
                                poly_kwargs={'facecolor': (0.73, 0.8, 1)})

        # Plot reference
        c1 = 0.125
        xx = litterature.degree
        yy = c1*xx**3
        ax1.loglog(xx, yy, '-', label=r'$O(C p^3)$')

        # Properties
        for ax in [ax1]:
            ax.legend()
            ax.set_xlim([2, 10])
            ax.set_ylim([1e-1, 1e3])
            ax.set_ylabel('CPU setup time (s)')
            ax.set_xlabel('Polynomial degree p')

        fig.tight_layout()
        fig.savefig(filename_WQ)

    if nbel == 64:
        # Load data
        mydata = np.loadtxt(filename_WQ)
        litterature = pd.read_table(folder + 'lit_MF.dat', sep='\t', names=['degree', 'C', 'K'])

        # Plot data
        fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(12,4))
        ax1.loglog(litterature.degree, litterature.C, 'o--', label='Litterature')
        ax1.loglog(mydata[:,0], mydata[:,1], 'x--', label='My algorithm')

        ax2.loglog(litterature.degree, litterature.K, 'o--', label='Litterature')
        ax2.loglog(mydata[:,0], mydata[:,2], 'x--', label='My algorithm')

        # Plot reference
        c1, c2 = 0.125, 0.125
        xx = litterature.degree
        yy = c1*xx**3
        yy2 = c2*xx**3
        ax1.loglog(xx, yy, '-', label=r'$O(C p^3)$')
        ax2.loglog(xx, yy2, '-', label=r'$O(C p^3)$')

        # Properties
        for ax in [ax1, ax2]:
            ax.legend()
            ax.set_xlim([2, 10])
            ax.set_ylim([1, 1e4])
            ax.set_ylabel('CPU setup time (s)')
            ax.set_xlabel('Polynomial degree p')

        fig.tight_layout()
        fig.savefig(filename_WQ)

        


