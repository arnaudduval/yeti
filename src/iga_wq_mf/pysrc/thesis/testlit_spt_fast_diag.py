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
quadArgs      = {'quadrule': 'iga', 'type': 'leg'}
dataExist     = True
degree_list   = range(2, 6)
cuts_list     = range(4, 8)

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
			V = np.random.random(nb_ctrlpts**4)

			inpts = [nbqp, nbqp, nbqp, nbqp, *indices, *indices, *indices, *indices, 
					dersbasis, dersbasis, dersbasis, dersbasis, dersweights, dersweights, dersweights, dersweights]
			
			start = time.time()
			eigensolver.sptfastdiagonalization_3d(*inpts, V)
			stop  = time.time()
			time_t = stop - start

			print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, time_t))
			time_matrix[i, j+1] = time_t
			# np.savetxt(folder2find + 'spt_FD_time.dat', time_matrix)

else:

	fig, ax = plt.subplots()

	# Load data
	file = pd.read_table(folder2find + 'spt_FD_time.dat', sep=' ', names=['nbel', 'p2', 'p3', 'p4', 'p5', 'p6'])
	CPUtime = np.column_stack((file.p2, file.p3, file.p4, file.p5))
	
	for i in range(4): ax.loglog(file.nbel**4, CPUtime[:, i], '--', label='degree ' + r'$p=$' + str(i+2), marker=MARKERLIST[i])
	# slope = np.polyfit(np.log10(file.nbel),np.log10(CPUtime[:, 2]), 1)[0]
	# slope = round(slope, 3)
	# annotation.slope_marker((file.nbel[1], CPUtime[1, 2]), slope, 
	# 						poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

	ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
	ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
	ax.set_xticks([1e4, 16**4, 32**4, 64**4, 128**4, 1e9])
	ax.set_xticklabels([r'$10^4$', r'$16^4$', r'$32^4$', r'$64^4$', r'$128^4$', r'$10^9$'])
	ax.minorticks_off()

	# ax.set_xticks([16**4, 32**4, 64**4, 128**4])
	# ax.set_xticklabels([r'$16^4$', r'$32^4$', r'$64^4$', r'$128^4$'])

	ax.legend(loc='best')
	ax.set_xlabel('Total number of elements')
	ax.set_ylabel('Wall time (s)')
	ax.set_ylim([0.01, 500])
	ax.set_xlim([1e4, 1e9])
	fig.tight_layout()
	fig.savefig(folder2save + 'sptFastDiag' + '.pdf')
