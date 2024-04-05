"""
.. Test of fast diagonalization
.. We verify speed of fast diagonalization as direct method to solve Sylvester system
.. TO DO:
	rebuild the functions to compute fast diagonalization 
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
quadArgs    = {'quadrule': 'iga', 'type': 'leg'}
dataExist     = False
withReference = True
degree_list   = range(2, 6)
cuts_list     = range(6, 10)

if not dataExist:

	time_matrix = np.zeros((len(cuts_list), len(degree_list)+1))
	time_matrix[:, 0] = np.array([2**i for i in cuts_list])

	for i, cuts in enumerate(cuts_list):
		for j, degree in enumerate(degree_list):
			nbel = 2**cuts
			knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel)
			weightedQuad = WeightedQuadrature(degree, knotvector, quadArgs={})
			info = weightedQuad.getQuadratureRulesInfo()
			qp, indices, dersbasis, dersweights = info; nbqp = len(qp)
			nb_ctrlpts = weightedQuad.nbctrlpts
			V = np.random.random(nb_ctrlpts**3)

			inpts = [nbqp, nbqp, nbqp, *indices, *indices, *indices, 
					dersbasis, dersbasis, dersbasis, dersweights, dersweights, dersweights]
			
			start = time.time()
			eigensolver.fastdiagonalization_3d(*inpts, V)
			stop  = time.time()
			time_t = stop - start

			print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, time_t))
			time_matrix[i, j+1] = time_t
			# np.savetxt(folder2find + 'FD_time.dat', time_matrix)

else:

	fig, ax = plt.subplots()

	# Load data
	file = pd.read_table(folder2find + 'FD_time.dat', sep=' ', names=['nbel', 'p2', 'p3', 'p4', 'p5', 'p6'])
	CPUtime = np.column_stack((file.p2, file.p3, file.p4, file.p5, file.p6))
	
	for i in range(5): ax.loglog(file.nbel, CPUtime[:, i], '--', label='degree ' + r'$p=$' + str(i+2), marker=MARKERLIST[i])
	slope = np.polyfit(np.log10(file.nbel),np.log10(CPUtime[:, 2]), 1)[0]
	slope = round(slope, 3)
	annotation.slope_marker((file.nbel[1], CPUtime[1, 2]), slope, 
							poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

	ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
	ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
	ax.set_xticks([64, 128, 256, 512])

	ax.legend(loc='best')
	ax.set_xlabel('Discretization level ' + r'$h^{-1}$')
	ax.set_ylabel('CPU time (s)')
	ax.set_ylim([0.01, 100])
	fig.tight_layout()
	fig.savefig(folder2save + 'FastDiag' + '.pdf')

	if withReference:
		fig, ax = plt.subplots()

		# Load data
		file = pd.read_table(folder2find + 'FD_time_lit.dat', sep='\t', names=['nbel', 'p2', 'p3', 'p4', 'p5', 'p6'])
		CPUtime = np.column_stack((file.p2, file.p3, file.p4, file.p5, file.p6))

		for i in range(np.shape(CPUtime)[1]): ax.loglog(file.nbel, CPUtime[:, i], 'o--', label='degree ' + r'$p=$' + str(i+2))

		ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
		ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		ax.set_xticks([128, 256, 512, 1024])
		ax.set_xticklabels([r'$128^3$', r'$256^3$', r'$512^3$', r'$1024^3$'])

		ax.legend(loc='best')
		ax.set_xlabel('Total number of elements')
		ax.set_ylabel('CPU time (s)')
		fig.tight_layout()
		fig.savefig(folder2save + 'FastDiag_lit' + '.pdf')