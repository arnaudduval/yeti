from thesis.Elliptic.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat, mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import heatproblem, mechaproblem

# Set global variables
RUNSIMU = False
degList = range(1, 10)
cuts = 6
quadArgs = {'quadrule':'wq', 'type':1}

if RUNSIMU:

	timeMF_conductivity = np.zeros((len(degList), 2))
	timeMF_conductivity[:, 0] = degList

	timeMF_stiffness = np.zeros((len(degList), 2))
	timeMF_stiffness[:, 0] = degList

	timeMF_conductivityPy = np.zeros((len(degList), 2))
	timeMF_conductivityPy[:, 0] = degList

	timeMF_stiffnessPy = np.zeros((len(degList), 2))
	timeMF_stiffnessPy[:, 0] = degList

	for i, deg in enumerate(degList):
		
		geoArgs = {'name': 'RQA', 'degree': deg*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)
		}
		blockPrint()			
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)
		dim = modelPhy.dim

		heatmaterial = heatmat()
		heatmaterial.addConductivity(inpt=1.0, isIsotropic=True, shape=dim)

		mecamaterial = mechamat({'elastic_modulus':1e3, 'elastic_limit':1e10, 
						'poisson_ratio':0.3, 'isoHardLaw': {'name':'none'}})

		# Set Dirichlet boundaries	
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		boundary.add_DirichletConstTemperature(table=np.ones((dim, 2), dtype=bool))
		boundary.add_DirichletDisplacement(table=np.ones((dim, 2, dim), dtype=bool))
		enablePrint()

		# Solve elastic problem
		hproblem = heatproblem(heatmaterial, modelPhy, boundary)
		mproblem = mechaproblem(mecamaterial, modelPhy, boundary)

		# ------------------
		# Compute MF product
		# ------------------
		enablePrint()
		print('******')

		start = time.process_time()
		hproblem.compute_mfConductivity(np.random.random(boundary._nbctrlpts_total))
		finish = time.process_time()
		print('Time Conductivity:%.2e' %(finish-start))
		timeMF_conductivity[i, 1] = finish - start

		start = time.process_time()
		mproblem.compute_mfStiffness(np.random.random((dim, boundary._nbctrlpts_total)))
		finish = time.process_time()
		print('Time Stiffness:%.2e' %(finish-start))
		timeMF_stiffness[i, 1] = finish - start

		np.savetxt(FOLDER2SAVE+'MF_conductivity_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_conductivity)
		np.savetxt(FOLDER2SAVE+'MF_stiffness_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_stiffness)

		if deg > 6: continue
		matrix = buildpseudomatrix_ht3d(hproblem)
		start = time.process_time()
		matrix @ np.random.random(boundary._nbctrlpts_total)
		finish = time.process_time()
		print('Time Conductivity Python:%.2e' %(finish-start))
		timeMF_conductivityPy[i, 1] = finish - start
		del matrix

		matrix = buildpseudomatrix_el3d(mproblem)
		start = time.process_time()
		matrix @ np.random.random((dim*boundary._nbctrlpts_total))
		finish = time.process_time()
		print('Time Stiffness Python:%.2e' %(finish-start))
		timeMF_stiffnessPy[i, 1] = finish - start
		del matrix

		np.savetxt(FOLDER2DATA+'MF_conductivity_Py'+'.dat', timeMF_conductivityPy)
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
	fig.savefig(FOLDER2SAVE+filename+'.pdf')