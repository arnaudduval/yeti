from pysrc.sptmpaper.input_data import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/' + GEONAME + '2/'
if not os.path.isdir(folder): os.mkdir(folder)
blockPrint()

# Set global variables
TODOSIMU = True
FIG_CASE = 4
extension = '.dat'

IgaPlot = {'marker': 's', 'linestyle': '-', 'markersize': 10}
WQ1Plot = {'marker': 'x', 'linestyle': '--', 'markersize': 6}
WQ2Plot = {'marker': 'o', 'linestyle': ':', 'markersize': 6}
IncPlot = {'marker': 'd', 'linestyle': '-.', 'markersize': 6}
IncPlot2 = {'marker': '*', 'linestyle': '-.', 'markersize': 6}

def exactTemperature_inc(args):
	func = None
	if GEONAME == 'QA': func = exactTemperatureRing_inc(args)
	elif GEONAME == 'TP': func = exactTemperatureTrap_inc(args)
	else: raise Warning('Not possible')
	return func

def exactTemperature_spt(args):
	func = None
	if GEONAME == 'QA': func = exactTemperatureRing_spt(args)
	elif GEONAME == 'TP': func = exactTemperatureTrap_spt(args)
	else: raise Warning('Not possible')
	return func

def powerDensity_inc(args):
	func = None
	if GEONAME == 'QA': func = powerDensityRing_inc(args)
	elif GEONAME == 'TP': func = powerDensityTrap_inc(args)
	else: raise Warning('Not possible')
	return func

def powerDensity_spt(args):
	func = None
	if GEONAME == 'QA': func = powerDensityRing_spt(args)
	elif GEONAME == 'TP': func = powerDensityTrap_spt(args)
	else: raise Warning('Not possible')
	return func

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
			filenameA1 = folder+'0L2abserror'+sufix+extension
			filenameR1 = folder+'0L2relerror'+sufix+extension
			for j, cuts in enumerate(cutList):
				for i, degree in enumerate(degList):
					geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
							'nb_refinementByDirection': np.array([cuts, cuts, 1])}
					dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
					problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, 
													dirichlet_table=dirichlet_table, quadArgs=quadArgs, 
													degree_time=degree, nbel_time=2**cuts, geoArgs=geoArgs)
					AbserrorTable[i+1, j+1], relerrorTable[i+1, j+1] = problem_spt.normOfError(temp_spt, 
																	normArgs={'type':'L2', 
																	'exactFunction':exactTemperature_spt},)

					np.savetxt(filenameA1, AbserrorTable)
					np.savetxt(filenameR1, relerrorTable)

	plotoptions = [IgaPlot, WQ1Plot, WQ2Plot]
	lastsufix = 'linear' if ISLINEAR else 'nonlin'
	figname = folder + 'Convergence'+lastsufix+'.pdf'
	filenames = ['0L2abserror_iga_leg_', '0L2abserror_wq_1_', '0L2abserror_wq_2_']

	fig, ax = plt.subplots(figsize=(7, 6))
	for filename, plotops in zip(filenames, plotoptions):
		quadrule = filename.split('_')[1]
		table = np.loadtxt(folder+filename+lastsufix+extension)	
		nbels = 2**(table[0, 1:])
		degList = table[1:, 0]
		errList  = table[1:, 1:]
		for i, degree in enumerate(degList):
			color = COLORLIST[i]
			if quadrule == 'iga': 
				ax.loglog(nbels, errList[i, :], label='ST-IGA-GL '+r'$p_s=p_t=$'+str(int(degree)), color=color, marker=plotops['marker'],
							markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])		
				slope = np.polyfit(np.log10(nbels[2:]),np.log10(errList[i, 2:]), 1)[0]
				slope = round(slope, 1)
				annotation.slope_marker((nbels[-2], errList[i, -2]), slope, 
								poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
			
			else: 
				ax.loglog(nbels, errList[i, :], color=color, marker=plotops['marker'], markerfacecolor='w',
						markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
			fig.savefig(figname)

	ax.loglog([], [], color='k', marker=WQ2Plot['marker'], markerfacecolor='w',
			markersize=WQ2Plot['markersize'], linestyle=WQ2Plot['linestyle'], label='ST-IGA-WQ 2')

	ax.loglog([], [], color='k', marker=WQ1Plot['marker'], markerfacecolor='w',
					markersize=WQ1Plot['markersize'], linestyle=WQ1Plot['linestyle'], label='ST-IGA-WQ 4')
	
	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}$')
	ax.set_xlabel('Number of elements by space-time direction')
	ax.set_xlim(left=1, right=50)
	ax.set_ylim(top=1e1, bottom=1e-8)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(figname)

elif FIG_CASE == 1:

	filenameA1 = folder + '1incheatAbs2'+lastsufix
	filenameR1 = folder + '1incheatRel2'+lastsufix
	filenameT1 = folder + '1incheatTim2'+lastsufix

	filenameA2 = folder + '1sptheatAbs2'+lastsufix
	filenameR2 = folder + '1sptheatRel2'+lastsufix
	filenameT2 = folder + '1sptheatTim2'+lastsufix

	degsptList = np.arange(1, 3, dtype=int)
	degList = np.array([1, 2, 3, 4, 5])
	cutList = np.arange(1, 7)

	if TODOSIMU:

		A1errorList = np.ones((len(degList), len(cutList)))
		R1errorList = np.ones((len(degList), len(cutList)))
		T1timeList = np.ones((len(degList), len(cutList)))

		A2errorList = np.ones((len(degList), len(cutList), len(degsptList)))
		R2errorList = np.ones((len(degList), len(cutList), len(degsptList)))
		T2timeList = np.ones((len(degList), len(cutList), len(degsptList)))
		
		for j, cuts in enumerate(cutList):
			for i, degree in enumerate(degList):
				geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, cuts, 1])}

				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt_inc = simulate_spacetime(degree, cuts, powerDensity_spt, 
													dirichlet_table=dirichlet_table, geoArgs=geoArgs, 
													quadArgs={'quadrule':'iga', 'type':'leg'},
													degree_time=1, nbel_time=2**5, solveSystem=False)[0]
				
				start = time.process_time()
				dirichlet_table = np.ones((2, 2))
				problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensity_inc, 
													dirichlet_table=dirichlet_table, geoArgs=geoArgs,
													nbel_time=2**5)
				finish = time.process_time()
				T1timeList[i, j] = finish - start
				
				A1errorList[i, j], R1errorList[i, j] = problem_spt_inc.normOfError(
														np.ravel(temp_inc, order='F'), 
														normArgs={'type':'L2',
																'exactFunction':exactTemperature_spt})
				
				np.savetxt(filenameA1+extension, A1errorList)
				np.savetxt(filenameR1+extension, R1errorList)
				np.savetxt(filenameT1+extension, T1timeList)

				for k, degspt in enumerate(degsptList):
					start = time.process_time()
					dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
					problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, 
														dirichlet_table=dirichlet_table, geoArgs=geoArgs, 
														quadArgs={'quadrule':'wq', 'type':2},
														degree_time=degspt, nbel_time=2**5)
					finish = time.process_time()
					T2timeList[i, j, k] = finish - start
					
					A2errorList[i, j, k], R2errorList[i, j, k] = problem_spt.normOfError(temp_spt, 
															normArgs={'type':'L2',
																	'exactFunction':exactTemperature_spt,})

					np.save(filenameA2, A2errorList)
					np.save(filenameR2, R2errorList)
					np.save(filenameT2, T2timeList)

	errorList1 = np.loadtxt(filenameA1+extension)
	errorList2 = np.load(filenameA2+'.npy')
	# fig, ax = plt.subplots(figsize=(7, 6))
	# for i, degree in enumerate(degList):
	# 	color = COLORLIST[i]
	# 	ax.loglog(2**cutList, errorList2[i, :, 0], color=color, marker=IgaPlot['marker'], markerfacecolor='w', 
	# 				markersize=IgaPlot['markersize'], linestyle=IgaPlot['linestyle'], 
	# 				label='ST-IGA-GL '+r'$p_s=$'+str(degree))
		
	# 	ax.loglog(2**cutList, errorList1[i, :], color=color, marker=IncPlot['marker'], markerfacecolor='w',
	# 				markersize=IncPlot['markersize'], linestyle=IncPlot['linestyle'])
		
	# ax.loglog([], [], color='k', marker=IncPlot['marker'], markerfacecolor='w',
	# 			markersize=IncPlot['markersize'], linestyle=IncPlot['linestyle'], label='INC-IGA-GL')

	# ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}$')
	# ax.set_xlabel('Number of elements by spatial direction')
	# ax.set_xlim(left=1, right=100)
	# ax.set_ylim(top=1e1, bottom=1e-5)
	# ax.legend(loc='lower left')
	# fig.tight_layout()
	# fig.savefig(folder + 'IncSptConv1'+  '.pdf')

	fig, ax = plt.subplots(figsize=(7, 6))
	for i, degree in enumerate(degList):
		color = COLORLIST[i]
		ax.loglog(2**cutList, errorList2[i, :, 1], color=color, marker=IgaPlot['marker'], markerfacecolor='w', 
					markersize=IgaPlot['markersize'], linestyle=IgaPlot['linestyle'], 
					label='ST-IGA-GL '+r'$p_s=$'+str(degree))
		
		ax.loglog(2**cutList, errorList1[i, :], color=color, marker=IncPlot['marker'], markerfacecolor='w',
					markersize=IncPlot['markersize'], linestyle=IncPlot['linestyle'])
	
	ax.loglog([], [], color='k', marker=IncPlot['marker'], markerfacecolor='w',
				markersize=IncPlot['markersize'], linestyle=IncPlot['linestyle'], label='INC-IGA-GL')
	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}$')
	ax.set_xlabel('Number of elements by spatial direction')
	ax.set_xlim(left=1, right=100)
	ax.set_ylim(top=1e1, bottom=1e-5)
	fig.tight_layout()
	fig.savefig(folder + 'IncSptConv2'+  '.pdf')

elif FIG_CASE == 2:

	degree, cuts = 8, 6
	quadArgs = {'quadrule':'wq', 'type':2}
	nbelincList = np.arange(2, 42, 4)
	degsptList = np.arange(1, 5)
	abserrorInc, relerrorInc = np.ones(len(nbelincList)), np.ones(len(nbelincList))
	abserrorSpt, relerrorSpt = np.ones((len(degsptList), len(nbelincList))), np.ones((len(degsptList), len(nbelincList)))

	if TODOSIMU:
		for i, nbelinc in enumerate(nbelincList):
			geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': np.array([cuts, cuts, 1])}

			dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
			problem_spt_inc = simulate_spacetime(degree, cuts, powerDensity_spt, 
												dirichlet_table=dirichlet_table, geoArgs=geoArgs, 
												quadArgs={'quadrule':'iga', 'type':'leg'},
												degree_time=1, nbel_time=nbelinc, solveSystem=False)[0]
			
			dirichlet_table = np.ones((2, 2))
			problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensity_inc, dirichlet_table=dirichlet_table,
														geoArgs=geoArgs, nbel_time=nbelinc, quadArgs=quadArgs, alpha=1.)
			
			abserrorInc[i], relerrorInc[i] = problem_spt_inc.normOfError(np.ravel(temp_inc, order='F'), 
										normArgs={'type':'L2',
												'exactFunction':exactTemperature_spt})
			
			np.savetxt(folder+'2abserrorstag_inc2'+extension, abserrorInc)
			np.savetxt(folder+'2relerrorstag_inc2'+extension, relerrorInc)

			for j, degspt in enumerate(degsptList):
				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, dirichlet_table=dirichlet_table,
													geoArgs=geoArgs, degree_time=degspt, nbel_time=nbelinc, quadArgs=quadArgs)
					
				abserrorSpt[j, i], relerrorSpt[j, i] = problem_spt.normOfError(temp_spt, 
														normArgs={'type':'L2',
																'exactFunction':exactTemperature_spt,})

				np.savetxt(folder+'2abserrorstag_spt'+extension, abserrorSpt)
				np.savetxt(folder+'2relerrorstag_spt'+extension, relerrorSpt)

	fig, ax = plt.subplots(figsize=(6, 5))
	errorList1 = np.loadtxt(folder+'2abserrorstag_spt'+extension)
	for i, deg in enumerate(degsptList):
		nbctrlpts = nbelincList+deg
		ax.loglog(nbctrlpts, errorList1[i, :], color=COLORLIST[i], marker=IgaPlot['marker'], markerfacecolor='w',
				markersize=IgaPlot['markersize'], linestyle=IgaPlot['linestyle'], label='ST-IGA-GL '+r'$p_t=$'+str(int(deg)))
		slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[i, 3:]), 1)[0]
		slope = round(slope, 1)
		annotation.slope_marker((nbctrlpts[-5], errorList1[i, -5]), slope, 
						poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

	errorList1 = np.loadtxt(folder+'2abserrorstag_inc'+extension)
	nbctrlpts = nbelincList+1
	ax.loglog(nbctrlpts, errorList1, marker=IncPlot['marker'], markerfacecolor='w', color='k',
					markersize=IncPlot['markersize'], linestyle=IncPlot['linestyle'], 
					label='INC-IGA-GL')
	slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[3:]), 1)[0]
	slope = round(slope, 1)
	annotation.slope_marker((nbctrlpts[-5], errorList1[-5]), slope, 
					poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
	
	# errorList1 = np.loadtxt(folder+'2abserrorstag_inc2'+extension)
	# ax.loglog(nbctrlpts, errorList1, marker=IncPlot2['marker'], markerfacecolor='w', color='k',
	# 				markersize=IncPlot2['markersize'], linestyle=IncPlot2['linestyle'], 
	# 				label='INC-IGA-GL '+r'$\alpha=1$')
	# slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[3:]), 1)[0]
	# slope = round(slope, 1)
	# annotation.slope_marker((nbctrlpts[-5], errorList1[-5]), slope, 
	# 				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
	
	ax.set_ylabel('Min. error in '+r'$L^2(\Pi)$' +'-norm')
	ax.set_xlabel('Number of control points in time')
	ax.set_xlim(left=2, right=50)
	ax.set_ylim(top=1e1, bottom=1e-8)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(folder+'StagnationError'+'.pdf')

elif FIG_CASE == 3:

	filenameA2 = folder + '32incheatAbs'
	filenameR2 = folder + '32incheatRel'
	filenameT2 = folder + '32incheatTim'

	filenameA3 = folder + '32sptheatAbs'
	filenameR3 = folder + '32sptheatRel'
	filenameT3 = folder + '32sptheatTim'

	degList = np.array([1, 2, 3, 4, 5])
	cutList = np.arange(4, 7)

	if TODOSIMU:
		# A2errorList = np.ones((len(degList), len(cutList)))
		# R2errorList = np.ones((len(degList), len(cutList)))
		# T2timeList = np.ones((len(degList), len(cutList)))
		
		# quadArgs = {'quadrule': 'wq', 'type': 2}
		# sufix = '_wq_2_' + lastsufix + extension
		# for j, cuts in enumerate(cutList):
		# 	for i, degree in enumerate(degList):
		# 		geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
		# 		'nb_refinementByDirection': np.array([cuts, cuts, 1])}

		# 		dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
		# 		problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, 
		# 											dirichlet_table=dirichlet_table, geoArgs=geoArgs, 
		# 											degree_time=1, nbel_time=2**cuts, quadArgs=quadArgs,
		# 											isadaptive=False, solveSystem=False)

		# 		start = time.process_time()
		# 		dirichlet_table = np.ones((2, 2))
		# 		problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensity_inc, 
		# 											dirichlet_table=dirichlet_table, geoArgs=geoArgs, nbel_time=2**cuts)
		# 		finish = time.process_time()
		# 		T2timeList[i, j] = finish - start
				
		# 		A2errorList[i, j], R2errorList[i, j] = problem_spt.normOfError(np.ravel(temp_inc, order='F'), 
		# 												normArgs={'type':'L2',
		# 														'exactFunction':exactTemperature_spt,})

		# 		np.savetxt(filenameA2+sufix, A2errorList)
		# 		np.savetxt(filenameR2+sufix, R2errorList)
		# 		np.savetxt(filenameT2+sufix, T2timeList)

		A3errorList = np.ones((len(degList), len(cutList)))
		R3errorList = np.ones((len(degList), len(cutList)))
		T3timeList = np.ones((len(degList), len(cutList)))
		
		for quadrule, quadtype in zip(['wq'], [2]):
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			sufix = '_' + quadrule + '_' + str(quadtype) + '_' + lastsufix + extension
			for j, cuts in enumerate(cutList):
				for i, degree in enumerate(degList):
					geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': np.array([cuts, cuts, 1])}

					start = time.process_time()
					dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
					problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, 
														dirichlet_table=dirichlet_table, geoArgs=geoArgs, 
														degree_time=degree, nbel_time=2**cuts, quadArgs=quadArgs,
														isadaptive=False)
						
					end = time.process_time()
					T3timeList[i, j] = end - start

					A3errorList[i, j], R3errorList[i, j] = problem_spt.normOfError(temp_spt, 
															normArgs={'type':'L2',
																	'exactFunction':exactTemperature_spt, })

					np.savetxt(filenameA3+sufix, A3errorList)
					np.savetxt(filenameR3+sufix, R3errorList)
					np.savetxt(filenameT3+sufix, T3timeList)

	position = 2
	assert position in [1, 2], 'Must be one or 2'
	if position==1: fig, ax = plt.subplots(figsize=(6, 4))
	if position==2: fig, ax = plt.subplots(figsize=(5.5, 4))
	# cmap = plt.get_cmap('jet', len(degList))
	cmap = mpl.colors.ListedColormap(COLORLIST[:len(degList)])

	Elist = np.loadtxt(filenameA2+'_wq_2_'+lastsufix+extension)
	Tlist = np.loadtxt(filenameT2+'_wq_2_'+lastsufix+extension)
	im = ax.scatter(Tlist[:len(degList), position], Elist[:len(degList), position], c=degList,
			cmap=cmap, marker=IncPlot['marker'], s=10*IncPlot['markersize'])
	ax.loglog(Tlist[:len(degList), position], Elist[:len(degList), position], 
			color='k', marker='', linestyle=IncPlot['linestyle'])
	
	if position==1:
		cbar = fig.colorbar(im); cbar.set_label('Degree')
		tick_locs = 1+(np.arange(len(degList)) + 0.5)*(len(degList)-1)/len(degList)
		cbar.set_ticks(tick_locs)
		cbar.set_ticklabels(degList)

	for quadrule, quadtype, plotvars in zip(['iga', 'wq'], ['leg', 2], [IgaPlot, WQ2Plot]):
		sufix = '_' + quadrule + '_' + str(quadtype) + '_' + lastsufix + extension
		Elist = np.loadtxt(filenameA3+sufix)
		Tlist = np.loadtxt(filenameT3+sufix)
		ax.scatter(Tlist[:len(degList), position], Elist[:len(degList), position], c=degList,
						cmap=cmap, marker=plotvars['marker'], s=10*plotvars['markersize'])
			
		ax.loglog(Tlist[:len(degList), position], Elist[:len(degList), position], 
				color='k', marker='', linestyle=plotvars['linestyle'])

	ax.loglog([], [], color='k', marker=IncPlot['marker'], alpha=0.5,
			markersize=IncPlot['markersize'], linestyle=IncPlot['linestyle'], label='INC-IGA-WQ 2')
	
	ax.loglog([], [], color='k', marker=IgaPlot['marker'], alpha=0.5,
		markersize=IgaPlot['markersize'], linestyle=IgaPlot['linestyle'], label='ST-IGA-GL')
		
	ax.loglog([], [], color='k', marker=WQ2Plot['marker'], alpha=0.5,
		markersize=WQ2Plot['markersize'], linestyle=WQ2Plot['linestyle'], label='ST-IGA-WQ 2')

	ax.set_xlabel('Wall time (s)')
	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}$')
	if position==1: ax.legend(loc='lower left')
	if position==1: ax.set_xlim(left=1e0, right=1e3)
	if position==2: ax.set_xlim(left=1e1, right=1e4)
	ax.set_ylim(top=1e-1, bottom=1e-9)
	fig.tight_layout()
	fig.savefig(folder + 'SPTINC_CPUError' + str(position) +  '.pdf')

elif FIG_CASE == 4:

	degree, cuts = 3, 4
	subfolderfolder = folder + str(degree) + '_' + str(cuts) + '/' 
	if not os.path.isdir(subfolderfolder): os.mkdir(subfolderfolder)

	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, cuts, 1])}

	if TODOSIMU:
		for [i, isadaptive], prefix1 in zip(enumerate([False, True]), ['exact', 'inexact']):
			for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
				prefix = prefix1 + '_' + prefix2 + '_'

				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, _, _, output = simulate_spacetime(degree, cuts, powerDensity_spt, 
												dirichlet_table=dirichlet_table,
												degree_time=degree, nbel_time=2**cuts, geoArgs=geoArgs,
												isadaptive=isadaptive, isfull=isfull,
												getOthers=True)
				
				sol_list  = output['Solution']; resKrylov_list = output['KrylovRes']
				resNewton_list = output['NewtonRes']; threshold_list = output['Threshold']
				L2error, L2relerror = [], []

				for solution in sol_list:
					err, relerr  = problem_spt.normOfError(solution, normArgs={'type':'L2', 
																	'exactFunction':exactTemperature_spt})
					L2error.append(err); L2relerror.append(relerr)
				resKrylovclean = np.array([]); counter_list = [0]
				for _ in resKrylov_list: 
					resKrylovclean = np.append(resKrylovclean, _[np.nonzero(_)])
					counter_list.append(counter_list[-1] + len(_[np.nonzero(_)]))

				np.savetxt(subfolderfolder+prefix+'CumulKrylovRes'+extension, resKrylovclean)
				np.savetxt(subfolderfolder+prefix+'Inner_loops'+extension, counter_list)
				np.savetxt(subfolderfolder+prefix+'NewtonRes'+extension, resNewton_list)
				np.savetxt(subfolderfolder+prefix+'L2error'+extension, L2error)
				np.savetxt(subfolderfolder+prefix+'L2relerror'+extension, L2relerror)
				np.savetxt(subfolderfolder+prefix+'Threshold'+extension, np.array(threshold_list))

	fig1, ax1 = plt.subplots(figsize=(6, 4))
	fig2, ax2 = plt.subplots(figsize=(6, 4))
	fig3, ax3 = plt.subplots(figsize=(6, 4))
	fig4, ax4 = plt.subplots(figsize=(6, 4))
	figs = [fig1, fig2, fig3, fig4]; axs  = [ax1, ax2, ax3, ax4]
	linestyle_list = ['-', '--', '-', '--']
	marker_list = ['o', 'o', 's', 's']

	for [i, isadaptive], prefix1 in zip(enumerate([True, False]), ['inexact', 'exact']):
		for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
			l = j + i*2
			legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
			prefix = prefix1 + '_' + prefix2 + '_'
			nbInnerLoops = np.loadtxt(subfolderfolder+prefix+'Inner_loops'+extension)
			newtonRes = np.loadtxt(subfolderfolder+prefix+'NewtonRes'+extension)
			L2relerror = np.loadtxt(subfolderfolder+prefix+'L2relerror'+extension)
			newtonRes = newtonRes/newtonRes[0]
			
			ylim1 = np.power(10, np.floor(np.log10(np.min(L2relerror))))
			ylim2 = np.power(10, np.floor(np.log10(np.min(newtonRes))))
			ylim = np.min([ylim1, ylim2])
			for caseplot, fig, ax in zip(range(1, 5), figs, axs):
				if caseplot == 1:
					yy = L2relerror; xx = nbInnerLoops[:len(L2relerror)]
					xlim = 10*np.ceil(np.max(nbInnerLoops)/10)
					ylabel = 'Relative '+r'$L^2(\Pi)$' + ' norm of error'
					xlabel = 'Number of inner iterations'
				elif caseplot == 2:
					yy = newtonRes; xx = nbInnerLoops[:len(newtonRes)]
					xlim = 10*np.ceil(np.max(nbInnerLoops)/10)
					ylabel = 'Relative norm of outer residue'
					xlabel = 'Number of inner iterations'
				elif caseplot == 3:
					yy = newtonRes; xx = np.arange(0, len(newtonRes))
					xlim = 8
					ylabel = 'Relative norm of outer residue'
					xlabel = 'Number of outer iterations'
				elif caseplot == 4:
					yy = L2relerror; xx = np.arange(0, len(L2relerror))
					xlim = 8
					ylabel = 'Relative '+r'$L^2(\Pi)$' + ' norm of error'
					xlabel = 'Number of outer iterations'

				ax.semilogy(xx, yy, label=legendname, marker=marker_list[l], linestyle=linestyle_list[l])
				ax.set_xlim(right=xlim, left=0)
				ax.set_ylim(top=1e1, bottom=ylim)
				ax.set_xlabel(xlabel)
				ax.set_ylabel(ylabel)
				ax.legend()
				fig.tight_layout()
				fig.savefig(folder+'NLConvergence_iters'+'_'+str(degree)+str(cuts)+str(caseplot)+'.pdf')

	fig, ax = plt.subplots(figsize=(6, 4))
	for [i, isadaptive], prefix1 in zip(enumerate([True]), ['inexact']):
		for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
			l = j + i*2
			legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
			prefix = prefix1 + '_' + prefix2 + '_'
			threshold = np.loadtxt(subfolderfolder+prefix+'Threshold'+extension)

			yy = threshold; xx = np.arange(0, len(threshold))
			xlim = 8
			ylabel = 'Forcing term'
			xlabel = 'Number of outer iterations'

			ax.semilogy(xx, yy, label=legendname, marker=marker_list[l], linestyle=linestyle_list[l])
			ax.set_xlim(right=xlim, left=0)
			# ax.set_ylim(top=ylim[0], bottom=ylim[1])
			ax.set_ylim(top=1, bottom=1e-4)
			ax.set_xlabel(xlabel)
			ax.set_ylabel(ylabel)
			ax.legend()
			fig.tight_layout()
			fig.savefig(folder+'NLTolerance'+'_'+str(degree)+str(cuts)+'.pdf')