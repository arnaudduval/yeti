"""
.. Test of fast diagonalization
.. We verify speed of fast diagonalization as direct method to solve Sylvester system
.. TO DO: rebuild the functions to compute fast diagonalization 
.. Joaquin Cornejo 
"""

from thesis.IncrementalHeat.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_quadrules import WeightedQuadrature

# Set global variables
quadArgs = {'quadrule': 'iga', 'type': 'leg'}
degList  = range(2, 6)
cutsList = range(5, 10)
RUNSIMU  = False
filename = FOLDER2SAVE + 'FD_CPUtime.dat'

if RUNSIMU:

	time_matrix = np.zeros((len(cutsList), len(degList)+1))
	time_matrix[:, 0] = np.array([2**i for i in cutsList])

	for i, cuts in enumerate(cutsList):
		for j, degree in enumerate(degList):
			nbel = 2**cuts
			knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel)
			weightedQuad = WeightedQuadrature(degree, knotvector, quadArgs={})
			info = weightedQuad.getQuadratureRulesInfo()
			qp, indices, dersbasis, dersweights = info
			nbqp = weightedQuad.nbqp
			nb_ctrlpts = weightedQuad.nbctrlpts
			v_in = np.random.random(nb_ctrlpts**3)

			inpts = [nbqp, nbqp, nbqp, 
					*indices, *indices, *indices, 
					dersbasis, dersbasis, dersbasis, 
					dersweights, dersweights, dersweights]
			
			start = time.process_time()
			eigensolver.fastdiagonalization_3d(*inpts, v_in)
			stop  = time.process_time()
			time_t = stop - start

			print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, time_t))
			time_matrix[i, j+1] = time_t
			np.savetxt(filename, time_matrix)

else:

	fig, ax = plt.subplots()
	time_matrix = np.loadtxt(filename)
	
	for i in range(4): 
		ax.loglog(time_matrix[:, 0]**3, time_matrix[:, i+1], '-', marker=MARKERLIST[i],
				label='degree ' + r'$p=$' + str(i+2))
	slope = round(np.polyfit(np.log10(time_matrix[:, 0]**3),np.log10(time_matrix[:, 2]), 1)[0], 2)
	annotation.slope_marker((time_matrix[2, 0]**3, time_matrix[2, -1]), slope, 
							poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

	ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
	ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
	ax.set_xticks([1e4, 32**3, 64**3, 128**3, 256**3, 512**3, 1e9])
	ax.set_xticklabels([r'$10^4$', r'$32^3$', r'$64^3$', r'$128^3$', r'$256^3$', r'$512^3$', r'$10^9$'])
	ax.minorticks_off()

	ax.legend(loc='best')
	ax.set_xlabel('Total number of elements')
	ax.set_ylabel('CPU time (s)')
	ax.set_ylim([0.01, 1e2])
	ax.set_xlim([1e4, 1e9])
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE + 'FastDiag' + '.pdf')