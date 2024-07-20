from thesis.Elliptic.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_quadrules import WeightedQuadrature

# Set global variables
TYPESIMU = 'meca'
RUNSIMU = False
degList = np.arange(1, 6)
cutList = np.arange(5, 9)
NDIM = 3
filename = FOLDER2SAVE + 'FD_time'

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

			inpts = [*[nbqp for i in range(NDIM)], *[indices[j] for i in range(NDIM) for j in range(2)], 
					*[dersbasis for i in range(NDIM)], *[dersweights for i in range(NDIM)]]

			start = time.process_time()
			if TYPESIMU == 'heat':
				vector_in = np.random.random(nb_ctrlpts**NDIM)
				if NDIM == 2: eigensolver.fastdiagonalization_2d(*inpts, vector_in)
				if NDIM == 3: eigensolver.fastdiagonalization_3d(*inpts, vector_in)
			elif TYPESIMU == 'meca':
				vector_in = np.random.random((NDIM, nb_ctrlpts**NDIM))
				if NDIM == 2: eigensolver.elfastdiagonalization_2d(*inpts, vector_in)
				if NDIM == 3: eigensolver.elfastdiagonalization_3d(*inpts, vector_in)

			finish = time.process_time()
			timeElapsed = finish - start

			print('For p = %s, nbel = %s, time: %.4f' %(degree, nbel, timeElapsed))
			timeMatrix[i+1, j+1] = timeElapsed
			np.savetxt(filename+TYPESIMU+'.dat', timeMatrix)

# Load data
fig, ax = plt.subplots(figsize=(6, 4))

file = np.loadtxt(filename+'heat'+'.dat')
degList = file[0, 1:]; nbelList = file[1:, 0]; timeElapsed = file[1:, 1:]
for i, degree in enumerate(degList): 
	ax.loglog(nbelList**3, timeElapsed[:, i], marker='s', color=COLORLIST[0], alpha=(i+1)/len(degList))
ax.loglog([], [], marker='s', color=COLORLIST[0], label='Steady heat')
slope = round(np.polyfit(np.log(nbelList**3), np.log(timeElapsed[:, 2]), 1)[0], 1)
annotation.slope_marker((nbelList[-2]**3,  timeElapsed[-2, 2]), slope, 
				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

file = np.loadtxt(filename+'meca'+'.dat')
degList = file[0, 1:]; nbelList = file[1:, 0]; timeElapsed = file[1:, 1:]
for i, degree in enumerate(degList): 
	ax.loglog(nbelList**3, timeElapsed[:, i], marker='o', color=COLORLIST[1], alpha=(i+1)/len(degList))
ax.loglog([], [], marker='o', color=COLORLIST[1], label='Elasticity')

ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

if NDIM == 2: ...
if NDIM == 3: 
	ax.set_xticks([1e4, 32**3, 64**3, 128**3, 256**3, 1e8])
	ax.set_xticklabels([r'$10^3$', r'$32^3$', r'$64^3$', r'$128^3$', r'$256^3$', r'$10^8$'])

ax.minorticks_off()
ax.legend(loc='best')
ax.set_xlabel('Total number of elements')
ax.set_ylabel('CPU time (s)')
ax.set_ylim([1e-3, 1e2])
fig.tight_layout()
fig.savefig(filename+'.pdf')
