"""
.. Test of setup time
.. We test how much time it takes to compute K and C matrices
.. TO DO:
	rebuild the functions to assemble the matrices
.. Joaquin Cornejo 
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_quadrules import WeightedQuadrature

# Select folder
full_path = os.path.realpath(__file__)
folder2save = os.path.dirname(full_path) + '/results/biblio/'
folder2find = os.path.dirname(full_path) + '/data/'

# Set global variables
dataExist   = True
nbel        = 64 # or 40 or 64
degree_list = np.arange(2, 8)

# Set filename
filename_data   = folder2find + 'setup_time_' + str(nbel) + '.dat'
filename_figure = folder2save + 'setup_time_' + str(nbel) + '.pdf' 

if not dataExist:

	time_matrix = np.zeros((len(degree_list), 3))
	time_matrix[:, 0] = degree_list
	
	for i, degree in enumerate(degree_list):

		knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel)
		weightedQuad = WeightedQuadrature(degree, knotvector)
		info = weightedQuad.getQuadratureRulesInfo()
		qp, [indi, indj], dersbasis, dersweights = info

		nbqp, indices, basis, wgt = [], [], [], []
		for dim in range(3):
			nbqp.append(len(qp))
			indices.append(indi); indices.append(indj) 
			basis.append(dersbasis); wgt.append(dersweights)

		# ----------------
		# Capacity matrix
		# ----------------
		nnz_I_list, nnz = np.array([-1, -1, -1], dtype=np.int32), 1
		coefs  = np.ones(len(qp)**3)
		inputs = [coefs, *nbqp, *indices, *basis, *wgt]

		start = time.process_time()
		assembly.wq_get_capacity_3d(*inputs, nnz_I_list, nnz)
		nnz = np.prod(nnz_I_list)
		info = assembly.wq_get_capacity_3d(*inputs, nnz_I_list, nnz)
		stop = time.process_time()
		time_matrix[i, 1] = stop - start         

		# --------------------
		# Conductivity matrix
		# --------------------
		nnz_I_list, nnz = np.array([-1, -1, -1], dtype=np.int32), 1
		coefs = np.ones((3, 3, len(qp)**3))
		inputs = [coefs, *nbqp, *indices, *basis, *wgt]

		start = time.process_time()
		assembly.wq_get_conductivity_3d(*inputs, nnz_I_list, nnz)
		nnz = np.prod(nnz_I_list)
		info = assembly.wq_get_conductivity_3d(*inputs, nnz_I_list, nnz)
		stop = time.process_time()
		time_matrix[i, 2] = stop - start 

		print(time_matrix[i, :])
		np.savetxt(filename_data, time_matrix)

else :

	if nbel == 40:
		fig, ax = plt.subplots()

		# Load data
		thiswork    = np.loadtxt(filename_data)
		litterature = pd.read_table(folder2find + 'setup_time_lit_40.dat', sep='\t', names=['degree', 'time'])

		ax.loglog(litterature.degree, litterature.time, 'o--', label='Literature')
		ax.loglog(thiswork[:, 0], thiswork[:, 1], 'x--', label='This work')

		# Compute slope
		slope = np.polyfit(np.log10(thiswork[:, 0]), np.log10(thiswork[:, 1]), 1)[0]
		slope = round(slope, 3)
		annotation.slope_marker((thiswork[2, 0], thiswork[2, 1]), slope, poly_kwargs={'facecolor': (0.73, 0.8, 1)})

		ax.legend()
		ax.set_xlim([2, 10])
		ax.set_ylim([1e-1, 1e3])
		ax.set_ylabel('CPU setup time (s)')
		ax.set_xlabel('Polynomial degree $p$')
		fig.tight_layout()
		fig.savefig(filename_figure)

	if nbel == 64:
		fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
		
		# Load data
		thiswork = np.loadtxt(filename_data)
		litterature = pd.read_table(folder2find + 'setup_time_lit_64.dat', sep='\t', names=['degree', 'C', 'K'])

		axs[0].loglog(litterature.degree, litterature.C, 'o--', label='Literature')
		axs[0].loglog(thiswork[:,0], thiswork[:,1], 'x--', label='This work')

		axs[1].loglog(litterature.degree, litterature.K, 'o--')
		axs[1].loglog(thiswork[:,0], thiswork[:,2], 'x--')

		for ax in axs:
			ax.set_xlim([2, 10])
			ax.set_ylim([1, 1e4])
			ax.set_ylabel('CPU setup time (s)')
			ax.set_xlabel('Polynomial degree $p$')

		fig.tight_layout()
		fig.savefig(filename_figure)