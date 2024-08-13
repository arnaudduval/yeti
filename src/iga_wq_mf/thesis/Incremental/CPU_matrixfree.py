from thesis.Incremental.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import heatmat, mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import heatproblem, mechaproblem

# Set global variables
RUNSIMU = False
degList = range(1, 6)
cuts = 6
quadArgs = {'quadrule':'iga', 'type':'leg'}
filename = 'MF_time' 

if RUNSIMU:

	timeMF_conductivity = np.zeros((len(degList), 2))
	timeMF_conductivity[:, 0] = degList

	timeMF_capacity = np.zeros((len(degList), 2))
	timeMF_capacity[:, 0] = degList

	for i, deg in enumerate(degList):
		
		geoArgs = {'name': 'VB', 'degree': deg*np.ones(3, dtype=int), 
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
		np.savetxt(FOLDER2DATA+'MF_conductivity_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_conductivity)

		start = time.process_time()
		hproblem.compute_mfCapacity(np.random.random(boundary._nbctrlpts_total))
		finish = time.process_time()
		print('Time Capacity:%.2e' %(finish-start))
		timeMF_capacity[i, 1] = finish - start
		np.savetxt(FOLDER2DATA+'MF_capacity_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_capacity)

fig, ax = plt.subplots()
plotoptions = [CONFIGLINE0, CONFIGLINE1, CONFIGLINE2]
sufixList = ['iga_leg', 'wq_1', 'wq_2']
labels = ['MF-GL', 'MF-WQ 1', 'MF-WQ 2']

# Load data
for j, [sufix, plotops, label] in enumerate(zip(sufixList, plotoptions, labels)):
	file_K1 = np.loadtxt(FOLDER2DATA+'MF_conductivity_'+sufix+'.dat') 
	file_C1 = np.loadtxt(FOLDER2DATA+'MF_capacity_'+sufix+'.dat') 

	degs = file_K1[:, 0]
	timeElapsed = file_K1[:, 1] + file_C1[:, 1]
	
	ax.semilogy(degs, timeElapsed, marker=plotops['marker'], label=label,
				markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])

ax.minorticks_off()
ax.legend(ncol=3, bbox_to_anchor=(0.5, 1.2), loc='upper center')
ax.set_xlabel('Degree ' + r'$p$')
ax.set_ylabel('CPU time (s)')
ax.set_xlim([0, 10])
ax.set_ylim([1e-1, 1e2])
fig.tight_layout()
fig.savefig(FOLDER2SAVE+filename+'.pdf')