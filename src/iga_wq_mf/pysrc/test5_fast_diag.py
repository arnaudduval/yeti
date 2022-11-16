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
								eigen_decomposition
)

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/test5/'
if not os.path.isdir(folder): os.mkdir(folder)
folder_data = os.path.dirname(full_path) + '/data/'

# Set global variables
dataExist     = True
withReference = True
degree_list   = range(2, 7)
cut_list      = range(6, 10)

# Set filename
filename_data = folder_data + 'FD_time.dat' 
extension     = '.png'

if not dataExist:

	time_matrix = np.zeros((len(cut_list), len(degree_list)+1))
	time_matrix[:, 0] = np.array([2**i for i in cut_list])

	for i, cuts in enumerate(cut_list):
		for j, degree in enumerate(degree_list):
		
			nbel = 2**cuts
			nb_ctrlpts = degree + nbel - 2
			knotvector = create_knotvector(degree, nbel)
			B, W, indi, indj = wq_find_basis_weights_fortran(degree, knotvector)[2:]

			# Erase data
			rows2erase = [0, -1]
			indi_t, indj_t, [B_t, W_t] = erase_rows_csr(rows2erase, indi, indj, [B, W])
			data_t = [B_t[:, 0], B_t[:, 1], W_t[:, 0], W_t[:, -1]]

			# Compute fast diagonalization
			V = np.random.random(nb_ctrlpts**3)
			start = time.process_time()

			eig_t, U_t = eigen_decomposition(indi_t, indj_t, data_t)
			eig_diag = np.random.random(nb_ctrlpts**3)
			fast_diagonalization(U_t, U_t, U_t, eig_diag, V, fdtype='steady')
			
			stop = time.process_time()
			time_t = stop - start

			print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, time_t))
			time_matrix[i, j+1] = time_t
			# np.savetxt(filename_data, time_matrix)

else:

	fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))

	# Load data
	file = pd.read_table(folder_data + 'FD_time.dat', sep=' ', names=['nbel', 'p2', 'p3', 'p4', 'p5', 'p6'])
	times = np.column_stack((file.p2, file.p3, file.p4, file.p5, file.p6))

	for i in range(5):
		ax.loglog(file.nbel, times[:, i], '--', label='degree ' + r'$p=$' + str(i+2), marker=markerSet[i])

	# Compute slope
	slope = np.polyfit(np.log10(file.nbel),np.log10(times[:, 2]), 1)[0]
	slope = round(slope, 3)
	annotation.slope_marker((file.nbel[1], times[1, 2]), slope, 
							poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

	ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
	ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
	ax.set_xticks([64, 128, 256, 512])

	ax.legend(loc='best')
	ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
	ax.set_ylabel('CPU time (s)')
	ax.set_ylim([0.01, 100])
	
	fig.tight_layout()
	fig.savefig(folder + 'FastDiag' + extension)

	if withReference:
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,4))

		# Load data
		file = pd.read_table(folder_data + 'FD_time_lit.dat', sep='\t', names=['nbel', 'p2', 'p3', 'p4', 'p5', 'p6'])
		times = np.column_stack((file.p2, file.p3, file.p4, file.p5, file.p6))

		# Plot data
		for i in range(np.shape(times)[1]):
			ax.loglog(file.nbel, times[:, i], 'o--', label='degree ' + r'$p=$' + str(i+2))

		ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
		ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		ax.set_xticks([128, 256, 512, 1024])

		ax.legend(loc='best')
		ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
		ax.set_ylabel('CPU time (s)')
		
		fig.tight_layout()
		fig.savefig(folder + 'FastDiag_lit' + '.png')