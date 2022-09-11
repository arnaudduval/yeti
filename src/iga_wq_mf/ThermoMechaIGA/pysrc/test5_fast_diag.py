"""
.. Test of fast diagonalization
.. We verify speed of fast diagonalization as direct method to solve Sylvester system
.. Joaquin Cornejo 
"""

from lib.__init__ import *
from lib.base_functions import erase_rows_csr, create_knotvector, wq_find_basis_weights_fortran
from iga_wq_mf import solver

# Select folder
full_path = os.path.realpath(__file__)
folder_file = os.path.dirname(full_path) + '/data/'
folder_figure = os.path.dirname(full_path) + '/results/test5/'
if not os.path.isdir(folder_figure): os.mkdir(folder_figure)

# Set global variable
dataExist = True
degree = 5

if not dataExist:
    
    for nbel in range(60, 80, 50):         
        # Define basis 
        knotvector = create_knotvector(degree, nbel)
        qp_wq, B, W, indi, indj = wq_find_basis_weights_fortran(degree, knotvector)[1:]

        # Erase data
        rows2erase = [0, -1]
        indi_t, indj_t, data_t = erase_rows_csr(rows2erase, indi, indj, [B, W])
        [B_t, W_t] = data_t

        # Initialize
        shape_matrices, indices, data_B, data_W = [], [], [], []
        for dim in range(3):
            shape_matrices.append(len(qp_wq))
            indices.append(indi_t); indices.append(indj_t) 
            data_B.append(B_t); data_W.append(W_t)

        # Solve sylvester equation P s = r
        inputs = [*shape_matrices, *indices, *data_B, *data_W]
        start = time.time()
        solver.test_precondfd(*inputs) # !!!!!! Maybe it not longer used
        stop = time.time()
        print(stop-start)

else:

    # Create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))

    # Import data
    file1 = pd.read_table(folder_file + 'TenP_new.dat', sep='\t', names=['nbel', 'time'])
    file2 = pd.read_table(folder_file + 'TenP_former.dat', sep='\t', names=['nbel', 'time'])
    file3 = pd.read_table(folder_file + 'TenP_lit.dat', sep='\t', names=['nbel', 'time'])

    # Post treatment
    files = [file1, file2, file3]
    labels = ['New algorithm', 'Former algorithm', 'Litterature']
    for file, label in zip(files, labels):

        # Extract data
        nbel = file.nbel**3
        times = file.time
        nbdata = len(nbel)

        # Plot data
        ax.loglog(nbel, times, 'o--', label=label)

        # Compute slope
        slope = np.polyfit(np.log10(nbel),np.log10(times), 1)[0]
        slope = round(slope, 3)
        annotation.slope_marker((nbel[round((nbdata-1)/2)], times[round((nbdata-1)/2)]), slope, 
                                text_kwargs={'fontsize': 14},
                                poly_kwargs={'facecolor': (0.73, 0.8, 1)})

    # Set properties
    ax.legend(loc='best')
    ax.set_xlabel("Total DOF")
    ax.set_ylabel("CPU time (s)")
    fig.tight_layout()
    fig.savefig(folder_figure + 'TensorProd' + '.png')

    # --------------------------

    # Create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))
    colors = ['#377eb8', '#ff7f00']

    # Import data
    file1 = pd.read_table(folder_file + 'MFProd.dat', sep='\t', names=['degree', 'Cu64', 'Cu128', 'Ku64', 'Ku128'])
    file2 = pd.read_table(folder_file + 'MFProd_lit.dat', sep='\t', names=['degree', 'Ku64', 'Cu64']) 

    # Litterature
    degree = file2.degree
    arrays = [file2.Cu64, file2.Ku64]
    labels = ['MF-WQ Mass Lit.', 'MF-WQ Stiffness Lit.']

    for array, label, color in zip(arrays, labels, colors):
        ax.semilogy(degree, array, 'x--', label=label, color=color)

    # My data
    degree = file1.degree
    arrays = [file1.Cu64, file1.Ku64]
    labels = ['MF-WQ Mass', 'MF-WQ Stiffness']

    for array, label, color in zip(arrays, labels, colors):
        ax.semilogy(degree, array, 'o--', label=label, color=color)

    # Set properties
    ax.legend(loc='best')
    ax.set_xlabel("Degree")
    ax.set_ylabel("CPU time (s)")
    ax.set_xlim([1, 10])
    ax.set_ylim([0.001, 10])
    fig.tight_layout()
    fig.savefig(folder_figure + 'ProductMF' + '.png')
