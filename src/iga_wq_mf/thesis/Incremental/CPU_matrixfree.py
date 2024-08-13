from thesis.Incremental.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat, mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import heatproblem, mechaproblem

# Set global variables
RUNSIMU = False
degList = range(9, 10)
cuts = 6
quadArgs = {'quadrule':'wq', 'type':1}
filename = FOLDER2SAVE + 'MF_time' 

if RUNSIMU:

	timeMF_conductivity = np.zeros((len(degList), 2))
	timeMF_conductivity[:, 0] = degList

	timeMF_capacity = np.zeros((len(degList), 2))
	timeMF_capacity[:, 0] = degList

	for i, degList in enumerate(degList):
		
		geoArgs = {'name': 'VB', 'degree': degList*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)
		}
		blockPrint()			
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)
		dim = modelPhy.dim

		heatmaterial = heatmat()
		heatmaterial.addCapacity(inpt=1.0, isIsotropic=True)
		heatmaterial.addConductivity(inpt=1.0, isIsotropic=True, shape=dim)

		# Set Dirichlet boundaries	
		boundary = boundaryCondition(modelPhy.nbctrlpts)
		boundary.add_DirichletConstTemperature(table=np.ones((dim, 2), dtype=bool))
		enablePrint()

		# Solve elastic problem
		hproblem = heatproblem(heatmaterial, modelPhy, boundary)

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
		# np.savetxt(FOLDER2SAVE+'MF_conductivity_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_conductivity)

		# start = time.process_time()
		# hproblem.compute_mfCapacity(np.random.random(boundary._nbctrlpts_total))
		# finish = time.process_time()
		# print('Time Capacity:%.2e' %(finish-start))
		# timeMF_capacity[i, 1] = finish - start
		# np.savetxt(FOLDER2SAVE+'MF_capacity_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_capacity)

# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(5.5, 5.5))
# plotoptions = [CONFIGLINE1, CONFIGLINE2]
# sufixList = ['iga_leg', 'wq_1', 'wq_2']
# labels = ['MF-GL', 'MF-WQ 1 conductivity', 'MF-WQ 1 capacity']

# # Load data
# for j, [sufix, plotops] in enumerate(zip(sufixList, plotoptions)):
# 	file_K1 = np.loadtxt(FOLDER2SAVE+'MF_conductivity_'+sufix+'.dat') 
# 	file_C1 = np.loadtxt(FOLDER2SAVE+'MF_capacity_'+sufix+'.dat') 

# 	degList = [file_K1[:, 0], file_C1[:, 0]]
# 	timeElapsedList = [file_K1[:, 1], file_C1[:, 1]]
# 	quadrule = sufix.split('_')[1]

# 	for i, [deg, timeElapsed, label] in enumerate(zip(degList, timeElapsedList, labels)):
# 		color = COLORLIST[i]
# 		if quadrule == '1':
# 			ax.semilogy(deg, timeElapsed, label = label, color=color, marker=plotops['marker'],
# 						markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])
# 		else: 
# 			ax.semilogy(deg, timeElapsed, color=color, marker=plotops['marker'],
# 						markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])

# ax.semilogy([], [], color='k', marker=CONFIGLINE2['marker'], markerfacecolor='w',
# 				markersize=CONFIGLINE2['markersize'], linestyle=CONFIGLINE2['linestyle'], label='MF-WQ 2')

# ax.minorticks_off()
# ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.2), loc='upper center')
# ax.set_xlabel('Degree ' + r'$p$')
# ax.set_ylabel('CPU time (s)')
# ax.set_xlim([0, 10])
# ax.set_ylim([1e-2, 50])
# fig.tight_layout()
# fig.savefig(filename+'.pdf')