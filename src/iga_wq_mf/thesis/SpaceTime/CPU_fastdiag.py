from thesis.SpaceTime.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_quadrules import WeightedQuadrature

# Set global variables
RUNSIMU = False
degList = np.arange(1, 6)
cutList = np.arange(4, 8)
NDIM = 4
filename = FOLDER2SAVE + 'sptFD_time'

if RUNSIMU:

	timeMatrix = np.zeros((len(cutList)+1, len(degList)+1))
	timeMatrix[1:, 0] = np.array([2**i for i in cutList])
	timeMatrix[0, 1:] = degList

	for i, cuts in enumerate(cutList):
		for j, degree in enumerate(degList):
			nbel = 2**cuts
			knotvector = createUniformKnotvector_Rmultiplicity(degree, nbel)
			weightedQuad = WeightedQuadrature(degree, knotvector, quadArgs={})
			info = weightedQuad.getQuadratureRulesInfo()
			qp, indices, dersbasis, dersweights = info; nbqp = len(qp)
			nb_ctrlpts = weightedQuad.nbctrlpts
			vector_in = np.random.random(nb_ctrlpts**NDIM)

			inpts = [*[nbqp for i in range(NDIM)], *[indices[j] for i in range(NDIM) for j in range(2)], 
					*[dersbasis for i in range(NDIM)], *[dersweights for i in range(NDIM)]]

			start = time.process_time()
			if NDIM == 3: eigensolver.sptfastdiagonalization_2d(*inpts, vector_in)
			if NDIM == 4: eigensolver.sptfastdiagonalization_3d(*inpts, vector_in)
			finish = time.process_time()
			timeElapsed = finish - start

			print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, timeElapsed))
			timeMatrix[i+1, j+1] = timeElapsed
			np.savetxt(filename+'.dat', timeMatrix)

# Load data
file = np.loadtxt(filename+'.dat')
degList = file[0, 1:]; nbelList = file[1:, 0]; timeElapsed = file[1:, 1:]

fig, ax = plt.subplots(figsize=(5, 4))
for i, degree in enumerate(degList): 
	ax.loglog(nbelList**4, timeElapsed[:, i], label=r'$p_s=p_t=$' + str(int(degree)), marker='s', color='k', alpha=(i+1)/len(degList))

slope = round(np.polyfit(np.log(nbelList**4), np.log(timeElapsed[:, 2]), 1)[0], 1)
annotation.slope_marker((nbelList[-2]**4,  timeElapsed[-2, 2]), slope, 
				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

if NDIM == 3: ...
if NDIM == 4: 
	ax.set_xticks([1e4, 16**4, 32**4, 64**4, 128**4, 1e9])
	ax.set_xticklabels([r'$10^4$', r'$16^4$', r'$32^4$', r'$64^4$', r'$128^4$', r'$10^9$'])

ax.minorticks_off()
ax.legend(loc='best')
ax.set_xlabel('Total number of elements')
ax.set_ylabel('CPU time (s)')
ax.set_ylim([1e-2, 1e3])
ax.set_xlim([1e4, 1e9])
fig.tight_layout()
fig.savefig(filename+'.pdf')