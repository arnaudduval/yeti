from thesis.Elliptic.__init__ import *
from pysrc.lib.lib_base import createUniformKnotvector_Rmultiplicity
from pysrc.lib.lib_quadrules import WeightedQuadrature, GaussQuadrature

# Set global variables
RUNSIMU = False
degList = range(1, 10)
NDIM = 3
nbel = 64

if RUNSIMU:

	timeMF_conductivity = np.zeros((len(degList), 2))
	timeMF_conductivity[:, 0] = degList

	timeMF_stiffness = np.zeros((len(degList), 2))
	timeMF_stiffness[:, 0] = degList

	timeMF_capacity = np.zeros((len(degList), 2))
	timeMF_capacity[:, 0] = degList

	timeMF_conductivityPy = np.zeros((len(degList), 2))
	timeMF_conductivityPy[:, 0] = degList

	timeMF_stiffnessPy = np.zeros((len(degList), 2))
	timeMF_stiffnessPy[:, 0] = degList

	timeMF_capacityPy = np.zeros((len(degList), 2))
	timeMF_capacityPy[:, 0] = degList

	for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
		
		quadArgs = {'quadrule':quadrule, 'type':quadtype}

		for i, deg in enumerate(degList):
			
			knotvector = createUniformKnotvector_Rmultiplicity(deg, nbel)
			if quadArgs['quadrule'] == 'wq':
				quadraturerule = WeightedQuadrature(deg, knotvector, quadArgs=quadArgs)
			if quadArgs['quadrule'] == 'iga':
				quadraturerule = GaussQuadrature(deg, knotvector, quadArgs=quadArgs)
			info = quadraturerule.getQuadratureRulesInfo()
			quadpts, indices, dersbasis, dersweights = info; nbqp = len(quadpts)
			nb_ctrlpts = quadraturerule.nbctrlpts

			inpts = [*[nbqp for i in range(NDIM)], *[indices[j] for i in range(NDIM) for j in range(2)], 
					*[dersbasis for i in range(NDIM)], *[dersweights for i in range(NDIM)]]

			nbctrlpts_total = np.product(np.array([nb_ctrlpts for i in range(NDIM)]))
			nbqp_total = np.product(np.array([nbqp for i in range(NDIM)]))
			invJ, detJ = np.ones((NDIM, NDIM, nbqp_total)), np.ones(nbqp_total)

			# ------------------
			# Compute MF product
			# ------------------
			enablePrint()
			print('******')

			if quadrule != 'iga' or deg < 9:
				prop = np.ones(nbqp_total)
				inpts_solver = [*inpts, False, invJ, detJ, prop, np.random.random(nbctrlpts_total)]
				start = time.process_time()
				if NDIM == 2: array_out = heatsolver.mf_capacity_2d(*inpts_solver)
				if NDIM == 3: array_out = heatsolver.mf_capacity_3d(*inpts_solver)
				finish = time.process_time()
				print('Time Capacity:%.2e' %(finish-start))
				timeMF_capacity[i, 1] = finish - start

				np.savetxt(FOLDER2DATA+'MF_capacity_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_capacity)

			if quadrule != 'iga' or deg < 6: 

				prop = np.ones((NDIM, NDIM, nbqp_total))
				inpts_solver = [*inpts, invJ, detJ, prop, np.random.random(nbctrlpts_total)]
				start = time.process_time()
				if NDIM == 2: array_out = heatsolver.mf_conductivity_2d(*inpts_solver)
				if NDIM == 3: array_out = heatsolver.mf_conductivity_3d(*inpts_solver)
				finish = time.process_time()
				print('Time Conductivity:%.2e' %(finish-start))
				timeMF_conductivity[i, 1] = finish - start
				del prop

				np.savetxt(FOLDER2DATA+'MF_conductivity_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_conductivity)

				prop = np.ones((2, nbqp_total))
				inpts_solver = [*inpts, invJ, detJ, prop, np.random.random((NDIM, nbctrlpts_total))]
				start = time.process_time()
				if NDIM == 2: array_out = plasticitysolver.mf_stiffness_2d(*inpts_solver)
				if NDIM == 3: array_out = plasticitysolver.mf_stiffness_3d(*inpts_solver)
				finish = time.process_time()
				print('Time Stiffness:%.2e' %(finish-start))
				timeMF_stiffness[i, 1] = finish - start
				del prop

				np.savetxt(FOLDER2DATA+'MF_stiffness_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_stiffness)

			if quadArgs['type'] != 1: continue

			if deg > 7: continue
			matrix = buildpseudomatrix_ht3d(quadraturerules=[quadraturerule for i in range(NDIM)], nbctrlpts_total=nbctrlpts_total)
			start = time.process_time()
			matrix @ np.random.random(nbctrlpts_total)
			finish = time.process_time()
			print('Time Conductivity Python:%.2e' %(finish-start))
			timeMF_conductivityPy[i, 1] = finish - start
			del matrix
			np.savetxt(FOLDER2DATA+'MF_conductivity_Py'+'.dat', timeMF_conductivityPy)

			if deg > 4: continue
			matrix = buildpseudomatrix_el3d(quadraturerules=[quadraturerule for i in range(NDIM)], nbctrlpts_total=nbctrlpts_total)
			start = time.process_time()
			matrix @ np.random.random((NDIM*nbctrlpts_total))
			finish = time.process_time()
			print('Time Stiffness Python:%.2e' %(finish-start))
			timeMF_stiffnessPy[i, 1] = finish - start
			del matrix

			np.savetxt(FOLDER2DATA+'MF_stiffness_Py'+'.dat', timeMF_stiffnessPy)

filenamelist = ['MF_conductivity_', 'MF_stiffness_']
sufixList = ['iga_leg', 'wq_1', 'wq_2']
labelList = ['MF-GL', 'MF-WQ 1', 'MF-WQ 2']
plotoptions = [CONFIGLINE0, CONFIGLINE1, CONFIGLINE2]

# Load data
for filename in filenamelist:
	fig, ax = plt.subplots()

	for label, sufix, plotops in zip(labelList, sufixList, plotoptions):
		file = np.loadtxt(FOLDER2DATA+filename+sufix+'.dat') 
		degList = file[:, 0]; timeElapsed = file[:, 1]
		ax.semilogy(degList, timeElapsed, label=label, marker=plotops['marker'],
					markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])

	file = np.loadtxt(FOLDER2DATA+filename+'Py'+'.dat') 
	degList = file[:, 0]; timeElapsed = file[:, 1]
	ax.semilogy(degList, timeElapsed, label='Built-in SpMdV')

	ax.minorticks_off()
	ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.2), loc='upper center')
	# ax.legend()
	ax.set_xlabel('Degree ' + r'$p$')
	ax.set_ylabel('CPU time (s)')
	ax.set_xlim([0, 10])
	ax.set_ylim([1e-1, 1e3])
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE+filename+'.png')