from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_quadrules import WeightedQuadrature

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/solver/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
TODOSIMU = False
degList = np.arange(1, 6)
cutList = np.arange(4, 8)
NDIM = 4
filename = folder + 'sptFD_time'

if TODOSIMU:

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

			start = time.time()
			if NDIM == 3: eigensolver.sptfastdiagonalization_2d(*inpts, vector_in)
			if NDIM == 4: eigensolver.sptfastdiagonalization_3d(*inpts, vector_in)
			finish = time.time()
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

ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

if NDIM == 3: ...
if NDIM == 4: 
	ax.set_xticks([1e4, 16**4, 32**4, 64**4, 128**4, 1e9])
	ax.set_xticklabels([r'$10^4$', r'$16^4$', r'$32^4$', r'$64^4$', r'$128^4$', r'$10^9$'])

ax.minorticks_off()
ax.legend(loc='best')
ax.set_xlabel('Total number of elements')
ax.set_ylabel('Wall time (s)')
ax.set_ylim([1e-2, 1e3])
ax.set_xlim([1e4, 1e9])
fig.tight_layout()
fig.savefig(filename+'.pdf')
