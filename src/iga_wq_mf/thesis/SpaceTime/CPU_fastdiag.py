from thesis.SpaceTime.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_quadrules import WeightedQuadrature

# Set global variables
RUNSIMU = False
degList = np.arange(1, 6)
cutList = np.arange(4, 8)
NDIM = 4
filename = 'sptFD_time'

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
			np.savetxt(FOLDER2DATA+filename+'.dat', timeMatrix)

# Make dummie mappable
fig, ax = plt.subplots()
cmap = mpl.colors.ListedColormap(COLORLIST[:len(degList)])
c = np.arange(1,len(degList)+1, dtype=int)
dummie_cax = ax.scatter(c, c, c=c, cmap=cmap)
ax.cla()

# Load data
file = np.loadtxt(FOLDER2DATA+filename+'.dat')
degList = file[0, 1:]; nbelList = file[1:, 0]; timeElapsed = file[1:, 1:]

for i, degree in enumerate(degList): 
	ax.loglog(nbelList**4, timeElapsed[:, i], marker='s')

slope = round(np.polyfit(np.log(nbelList**4), np.log(timeElapsed[:, 2]), 1)[0], 1)
annotation.slope_marker((nbelList[-2]**4,  timeElapsed[-2, 2]), slope, 
				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

cbar = plt.colorbar(dummie_cax)
cbar.set_label('Degree '+r'$p_s=p_t$')
tick_locs = 1+(np.arange(len(degList)) + 0.5)*(len(degList)-1)/len(degList)
cbar.set_ticks(tick_locs)
cbar.set_ticklabels(np.array(degList, dtype=int))

ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

if NDIM == 3: ...
if NDIM == 4: 
	ax.set_xticks([1e4, 16**4, 32**4, 64**4, 128**4, 1e9])
	ax.set_xticklabels([r'$10^4$', r'$16^4$', r'$32^4$', r'$64^4$', r'$128^4$', r'$10^9$'])

ax.minorticks_off()
ax.set_xlabel('Total number of elements')
ax.set_ylabel('CPU time (s)')
ax.set_ylim([1e-2, 1e3])
ax.set_xlim([1e4, 1e9])
fig.tight_layout()
fig.savefig(FOLDER2SAVE+filename+'.pdf')