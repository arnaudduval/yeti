from thesis.Elliptic.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat, mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import heatproblem, mechaproblem

# Set global variables
RUNSIMU = False
degList = range(1, 10)
cuts = 5
quadArgs = {'quadrule':'iga', 'type':'leg'}

if RUNSIMU:

	timeMF_conductivity = np.zeros((len(degList), 2))
	timeMF_conductivity[:, 0] = degList

	timeMF_stiffness = np.zeros((len(degList), 2))
	timeMF_stiffness[:, 0] = degList

	for i, degList in enumerate(degList):
		
		geoArgs = {'name': 'VB', 'degree': degList*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)
		}
		blockPrint()			
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)

		heatmaterial = heatmat()
		heatmaterial.addCapacity(inpt=1.0, isIsotropic=True)
		heatmaterial.addConductivity(inpt=1.0, isIsotropic=True, shape=3)

		mecamaterial = mechamat({'elastic_modulus':1e3, 'elastic_limit':1e10, 
						'poisson_ratio':0.3, 'isoHardLaw': {'name':'none'}})

		# Set Dirichlet boundaries	
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		boundary.add_DirichletConstTemperature(table=np.ones((3, 2), dtype=bool))
		boundary.add_DirichletDisplacement(table=np.ones((3, 2, 3), dtype=bool))
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
		mproblem.compute_mfStiffness(np.random.random((3, boundary._nbctrlpts_total)))
		finish = time.process_time()
		print('Time Stiffness:%.2e' %(finish-start))
		timeMF_stiffness[i, 1] = finish - start

		np.savetxt(FOLDER2SAVE+'MF_conductivity_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_conductivity)
		np.savetxt(FOLDER2SAVE+'MF_stiffness_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_stiffness)

fig, ax = plt.subplots(figsize=(6, 4))
IgaPlot = {'marker': 's', 'linestyle': '-', 'markersize': 10}
WQ1Plot = {'marker': 'x', 'linestyle': '--', 'markersize': 6}
WQ2Plot = {'marker': 'o', 'linestyle': ':', 'markersize': 6}
plotoptions = [IgaPlot, WQ1Plot, WQ2Plot]
sufixList = ['iga_leg', 'wq_1', 'wq_2']
labels = ['Steady heat', 'Elasticity']

# Load data
for sufix, plotops in zip(sufixList, plotoptions):
	file_K1 = np.loadtxt(FOLDER2SAVE+'MF_conductivity_'+sufix+'.dat') 
	file_S1 = np.loadtxt(FOLDER2SAVE+'MF_stiffness_'+sufix+'.dat') 

	degList = file_K1[:, 0]
	timeElapsedList = [file_K1[:, 1], file_S1[:, 1]]
	quadrule = sufix.split('_')[0]

	for i, [timeElapsed, label] in enumerate(zip(timeElapsedList, labels)):
		color = COLORLIST[i]
		if quadrule == 'iga':
			ax.semilogy(degList, timeElapsed, label='IGA-GL '+label, color=color, marker=plotops['marker'],
						markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])
		else:
			ax.semilogy(degList, timeElapsed, color=color, marker=plotops['marker'],
						markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])

ax.semilogy([], [], color='k', marker=WQ1Plot['marker'], markerfacecolor='w',
				markersize=WQ1Plot['markersize'], linestyle=WQ1Plot['linestyle'], label='IGA-WQ 1')
ax.semilogy([], [], color='k', marker=WQ2Plot['marker'], markerfacecolor='w',
				markersize=WQ2Plot['markersize'], linestyle=WQ2Plot['linestyle'], label='IGA-WQ 2')

ax.minorticks_off()
ax.legend(ncol=2, loc='upper center')
ax.set_xlabel('Degree ' + r'$p$')
ax.set_ylabel('CPU time (s)')
ax.set_xlim([0, 10])
ax.set_ylim([1e-1, 1e4])
fig.tight_layout()
fig.savefig(FOLDER2SAVE + 'MF_time' + '.pdf')