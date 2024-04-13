from pysrc.sptmpaper.input_data import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/1D/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
TODOSIMU = True
FIG_CASE = 0

IgaPlot = {'marker': 's', 'linestyle': '-', 'markersize': 10}
WQ1Plot = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
WQ2Plot = {'marker': 'x', 'linestyle': ':', 'markersize': 6}
IncPlot = {'marker': 'd', 'linestyle': '-.', 'markersize': 6}

lastsufix = 'linear' if ISLINEAR else 'nonlin'

if FIG_CASE == 0:
	if TODOSIMU:
		degList = np.array([1, 2, 3, 4, 5])
		cutList = np.arange(1, 6)
		for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
			sufix = '_' + quadrule + '_' + str(quadtype) + '_' + lastsufix
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			AbserrorTable = np.zeros((len(degList)+1, len(cutList)+1))
			relerrorTable = np.zeros((len(degList)+1, len(cutList)+1))
			AbserrorTable[0, 1:] = cutList; relerrorTable[0, 1:] = cutList
			AbserrorTable[1:, 0] = degList; relerrorTable[1:, 0] = degList
			filenameA1 = folder+'L2abserror_meshpar'+sufix+'.dat'
			filenameR1 = folder+'L2relerror_meshpar'+sufix+'.dat'
			for j, cuts in enumerate(cutList):
				for i, degree in enumerate(degList):
					blockPrint()
					dirichlet_table = np.ones((2, 2)); dirichlet_table[-1, 1] = 0
					problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensitySquare_spt, 
													dirichlet_table=dirichlet_table, quadArgs=quadArgs, 
													degree_time=degree, cuts_time=cuts, is1dim=IS1DIM)
					
					enablePrint()
					AbserrorTable[i+1, j+1], relerrorTable[i+1, j+1] = problem_spt.normOfError(temp_spt, 
																	normArgs={'type':'L2', 
																	'exactFunction':exactTemperatureSquare_spt},)

					np.savetxt(filenameA1, AbserrorTable)
					np.savetxt(filenameR1, relerrorTable)

	plotoptions = [IgaPlot, WQ1Plot, WQ2Plot]
	lastsufix = 'linear' if ISLINEAR else 'nonlin'
	figname = folder + 'SPTNLL2Convergence'+lastsufix+'.pdf'
	filenames = ['L2abserror_meshpar_iga_leg_', 'L2abserror_meshpar_wq_1_', 'L2abserror_meshpar_wq_2_']

	fig, ax = plt.subplots(figsize=(8, 6))
	for filename, plotops in zip(filenames, plotoptions):
		quadrule = filename.split('_')[2]
		table = np.loadtxt(folder+filename+lastsufix+'.dat')	
		nbels = 2**(table[0, 1:])
		degList = table[1:, 0]
		errList  = table[1:, 1:]
		for i, degree in enumerate(degList):
			color = COLORLIST[i]
			if quadrule == 'iga': 
				ax.loglog(nbels, errList[i, :], label='IGA-GL deg. '+str(int(degree)), color=color, marker=plotops['marker'],
							markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])		
			else: 
				ax.loglog(nbels, errList[i, :], color=color, marker=plotops['marker'], markerfacecolor='w',
						markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
			fig.savefig(figname)

	ax.loglog([], [], color='k', marker=WQ1Plot['marker'], markerfacecolor='w',
					markersize=WQ1Plot['markersize'], linestyle=WQ1Plot['linestyle'], label='IGA-WQ 4')
	ax.loglog([], [], color='k', marker=WQ2Plot['marker'], markerfacecolor='w',
			markersize=WQ2Plot['markersize'], linestyle=WQ2Plot['linestyle'], label='IGA-WQ 2')

	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}$')
	ax.set_xlabel('Number of elements by space-time direction')
	ax.set_xlim(left=1, right=100)
	ax.set_ylim(top=1e1, bottom=1e-10)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(figname)

elif FIG_CASE == 1:

	filenameA1 = folder + '1incheatAbs'+lastsufix
	filenameR1 = folder + '1incheatRel'+lastsufix
	filenameT1 = folder + '1incheatTim'+lastsufix

	filenameA2 = folder + '1sptheatAbs'+lastsufix
	filenameR2 = folder + '1sptheatRel'+lastsufix
	filenameT2 = folder + '1sptheatTim'+lastsufix

	degList = np.array([1, 2, 3, 4])
	cutList = np.arange(1, 6)

	if TODOSIMU:

		A1errorList = np.ones((len(degList), len(cutList)))
		R1errorList = np.ones((len(degList), len(cutList)))
		T1timeList = np.ones((len(degList), len(cutList)))

		A2errorList = np.ones((len(degList), len(cutList)))
		R2errorList = np.ones((len(degList), len(cutList)))
		T2timeList = np.ones((len(degList), len(cutList)))
		
		for j, cuts in enumerate(cutList):
			for i, degree in enumerate(degList):

				blockPrint()
				start = time.process_time()
				dirichlet_table = np.ones((2, 2))
				problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensitySquare_inc, 
													dirichlet_table=dirichlet_table, cuts_time=3, is1dim=IS1DIM)
				finish = time.process_time()
				T1timeList[i, j] = finish - start
				
				# Error of last step
				A1errorList[i, j], R1errorList[i, j] = problem_inc.normOfError(temp_inc[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperatureSquare_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})
				

				start = time.process_time()
				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensitySquare_spt, 
													dirichlet_table=dirichlet_table, 
													quadArgs={'quadrule':'wq', 'type':2},
													degree_time=2, cuts_time=3, is1dim=IS1DIM)
					
				finish = time.process_time()
				T2timeList[i, j] = finish - start

				# Error of last "step"
				newtemp_spt = np.reshape(temp_spt, newshape=(problem_spt.part.nbctrlpts_total, problem_spt.time.nbctrlpts_total), order='F')
				A2errorList[i, j], R2errorList[i, j] = problem_inc.normOfError(newtemp_spt[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperatureSquare_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})

				enablePrint()

			np.savetxt(filenameA1+'.dat', A1errorList)
			np.savetxt(filenameR1+'.dat', R1errorList)
			np.savetxt(filenameT1+'.dat', T1timeList)
			np.savetxt(filenameA2+'.dat', A2errorList)
			np.savetxt(filenameR2+'.dat', R2errorList)
			np.savetxt(filenameT2+'.dat', T2timeList)

	errorList1 = np.loadtxt(filenameA1+'.dat')
	errorList2 = np.loadtxt(filenameA2+'.dat')
	fig, ax = plt.subplots(figsize=(8, 6))
	for i, degree in enumerate(degList):
		color = COLORLIST[i]
		ax.loglog(2**cutList, errorList2[i, :], color=color, marker=IgaPlot['marker'], markerfacecolor='w',
					markersize=IgaPlot['markersize'], linestyle=IgaPlot['linestyle'], )
		
		ax.loglog(2**cutList, errorList1[i, :], color=color, marker=IncPlot['marker'], markerfacecolor='w',
					markersize=IncPlot['markersize'], linestyle=IncPlot['linestyle'], label='INC-IGA deg. '+str(degree))
		
	ax.loglog([], [], color='k', marker=IgaPlot['marker'], markerfacecolor='w', 
					markersize=IgaPlot['markersize'], linestyle=IgaPlot['linestyle'], label='ST-IGA-GL')

	
	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Omega)}$')
	ax.set_xlabel('Number of elements by spatial direction')
	ax.set_xlim(left=1, right=100)
	ax.set_ylim(top=2e1, bottom=1e-4)
	ax.legend(loc='upper right')
	fig.tight_layout()
	fig.savefig(folder + 'INCNLL2Convergence' +  '.pdf')
	plt.close(fig)

elif FIG_CASE == 2:

	filenameA1 = folder + 'incheatAbs'+lastsufix
	filenameR1 = folder + 'incheatRel'+lastsufix
	filenameT1 = folder + 'incheatTim'+lastsufix
	degree, cuts = 8, 6
	quadArgs = {'quadrule':'wq', 'type':2}
	cutsincList = np.arange(4, 7)
	degsptList = np.arange(1, 5)
	abserrorInc, relerrorInc = np.ones(len(cutsincList)), np.ones(len(cutsincList))
	abserrorSpt, relerrorSpt = np.ones((len(degsptList), len(cutsincList))), np.ones((len(degsptList), len(cutsincList)))

	if TODOSIMU:
		for i, cutsinc in enumerate(cutsincList):

			blockPrint()
			# Incremental
			dirichlet_table = np.ones((2, 2))
			problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensitySquare_inc, dirichlet_table=dirichlet_table,
														cuts_time=cutsinc, quadArgs=quadArgs, is1dim=IS1DIM)
			
			abserrorInc[i], relerrorInc[i] = problem_inc.normOfError(temp_inc[:, -1], 
										normArgs={'type':'L2',
													'exactFunction':exactTemperatureSquare_inc,
													'exactExtraArgs':{'time':time_inc[-1]}})
			
			np.savetxt(folder+'2abserrorstag_inc'+'.dat', abserrorInc)
			np.savetxt(folder+'2relerrorstag_inc'+'.dat', relerrorInc)

			# Space time
			for j, degspt in enumerate(degsptList):
				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensitySquare_spt, dirichlet_table=dirichlet_table,
													degree_time=degspt, cuts_time=cutsinc, quadArgs=quadArgs, is1dim=IS1DIM)
					
				newtemp_spt = np.reshape(temp_spt, newshape=(problem_spt.part.nbctrlpts_total, problem_spt.time.nbctrlpts_total), order='F')
				abserrorSpt[j, i], relerrorSpt[j, i] = problem_inc.normOfError(newtemp_spt[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperatureSquare_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})

				np.savetxt(folder+'2abserrorstag_spt'+'.dat', abserrorSpt)
				np.savetxt(folder+'2relerrorstag_spt'+'.dat', relerrorSpt)
			enablePrint()

	fig, ax = plt.subplots(figsize=(8, 6))
	errorList1 = np.loadtxt(folder+'2abserrorstag_spt'+'.dat')
	for i, deg in enumerate(degsptList):
		ax.loglog(2**cutsincList+1, errorList1[i, :], color=COLORLIST[i], marker=IgaPlot['marker'], markerfacecolor='w',
					markersize=IgaPlot['markersize'], linestyle=IgaPlot['linestyle'], label='SPT-IGA deg. '+str(int(deg)))

	errorList1 = np.loadtxt(folder+'2abserrorstag_inc'+'.dat')
	ax.loglog(2**cutsincList+1, errorList1, marker=IncPlot['marker'], markerfacecolor='w', color='k',
					markersize=IncPlot['markersize'], linestyle=IncPlot['linestyle'], label='INC-IGA')
	

	ax.set_ylabel('Stagnation error')
	ax.set_xlabel('Number of control points on time')
	# ax.set_xlim(left=1, right=100)
	# ax.set_ylim(top=1e1, bottom=1e-10)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(folder+'stagnationError'+'.pdf')