from thesis.SpaceTime.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job3d import stheatproblem

# Set global variables
TODOSIMU = False
degList = range(1, 7)
cuts = 4
# quadArgs = {'quadrule':'wq', 'type':2}
quadArgs = {'quadrule':'iga', 'type':'leg'}

if TODOSIMU:

	timeMF_Mass = np.zeros((len(degList), 2))
	timeMF_Mass[:, 0] = degList
	timeMF_Stiff = np.zeros((len(degList), 2))
	timeMF_Stiff[:, 0] = degList
	timeMF_dersMass = np.zeros((len(degList), 2))
	timeMF_dersMass[:, 0] = degList
	timeMF_dersStiff = np.zeros((len(degList), 2))
	timeMF_dersStiff[:, 0] = degList

	for i, degList in enumerate(degList):
		
		geoArgs = {'name': 'VB', 'degree': degList*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)
		}
		blockPrint()			
		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)
		timecrv  = part1D(createUniformOpenCurve(degList, int(2**cuts), 1.0), {'quadArgs':quadArgs})

		heatmaterial = heatmat()
		heatmaterial.addCapacity(inpt=1.0, isIsotropic=True)
		heatmaterial.addConductivity(inpt=1.0, isIsotropic=True, shape=3)
		heatmaterial.addCapacityDers(inpt=1.0, isIsotropic=True)
		heatmaterial.addConductivityDers(inpt=1.0, isIsotropic=True, shape=3)
	
		# Set Dirichlet boundaries	
		sptnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], timecrv.nbctrlpts_total])
		boundary = boundaryCondition(sptnbctrlpts)
		dirichlet_table = np.ones((4, 2)); dirichlet_table[-1, -1] = 0
		boundary.add_DirichletConstTemperature(table=dirichlet_table)
		enablePrint()

		# Solve elastic problem
		problem = stheatproblem(heatmaterial, modelPhy, timecrv, boundary)
		vector_in = np.random.random(boundary._nbctrlpts_total)

		# ------------------
		# Compute MF product
		# ------------------
		enablePrint()
		print('******')
		start = time.process_time()
		problem.compute_mfSTCapacity(vector_in)
		finish = time.process_time()
		print('Time Capacity:%.2e' %(finish-start))
		timeMF_Mass[i, 1] = finish - start

		start = time.process_time()
		problem.compute_mfSTConductivity(vector_in)
		finish = time.process_time()
		print('Time Conductivity:%.2e' %(finish-start))
		timeMF_Stiff[i, 1] = finish - start

		start = time.process_time()
		problem.compute_mfSTCapacityDers(vector_in)
		finish = time.process_time()
		print('Time Ders Capacity:%.2e' %(finish-start))
		timeMF_dersMass[i, 1] = finish - start

		start = time.process_time()
		problem.compute_mfSTConductivityDers(vector_in)
		finish = time.process_time()
		print('Time Ders Conductivity:%.2e' %(finish-start))
		timeMF_dersStiff[i, 1] = finish - start

		np.savetxt(FOLDER2SAVE+'sptMF_Mass_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_Mass)
		np.savetxt(FOLDER2SAVE+'sptMF_Stiff_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_Stiff)
		np.savetxt(FOLDER2SAVE+'sptMF_MassDers_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_dersMass)
		np.savetxt(FOLDER2SAVE+'sptMF_StiffDers_'+quadArgs['quadrule']+'_'+str(quadArgs['type'])+'.dat', timeMF_dersStiff)

fig, ax = plt.subplots(figsize=(5, 4))
IgaPlot = {'marker': 's', 'linestyle': '-', 'markersize': 10}
WQ1Plot = {'marker': 'x', 'linestyle': '--', 'markersize': 6}
WQ2Plot = {'marker': 'o', 'linestyle': ':', 'markersize': 6}
plotoptions = [IgaPlot, WQ2Plot]
sufixList = ['iga_leg', 'wq_2']
labels = [r'$\mathsf{A}v_{in}$', r'$\mathsf{B}v_{in}$']

# Load data
for sufix, plotops in zip(sufixList, plotoptions):
	file_M1 = np.loadtxt(FOLDER2SAVE+'sptMF_Mass_'+sufix+'.dat') 
	file_K1 = np.loadtxt(FOLDER2SAVE+'sptMF_Stiff_'+sufix+'.dat') 
	file_M2 = np.loadtxt(FOLDER2SAVE+'sptMF_MassDers_'+sufix+'.dat') 
	file_K2 = np.loadtxt(FOLDER2SAVE+'sptMF_StiffDers_'+sufix+'.dat') 

	degList = file_M1[:, 0]
	timeElapsedList = [file_M1[:, 1]+file_K1[:, 1], file_M2[:, 1]+file_K2[:, 1]]
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
				markersize=WQ1Plot['markersize'], linestyle=WQ1Plot['linestyle'], label='IGA-WQ 4')
ax.semilogy([], [], color='k', marker=WQ2Plot['marker'], markerfacecolor='w',
				markersize=WQ2Plot['markersize'], linestyle=WQ2Plot['linestyle'], label='IGA-WQ 2')

ax.minorticks_off()
ax.legend(ncol=2, loc='upper center')
ax.set_xlabel('Degree ' + r'$p$')
ax.set_ylabel('CPU time (s)')
ax.set_xlim([0, 10])
ax.set_ylim([1e-1, 1e4])
fig.tight_layout()
fig.savefig(FOLDER2SAVE + 'sptMF_time' + '.pdf')