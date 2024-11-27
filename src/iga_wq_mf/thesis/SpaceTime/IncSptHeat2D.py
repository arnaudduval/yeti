from thesis.SpaceTime.__init__ import *
from thesis.SpaceTime.input_data import *

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

# Set global variables
SUFIX = ('lin' if ISLINEAR else 'nonlin') + GEONAME
PLOTRELATIVE = True
RUNSIMU = False
FIG_CASE = 0
EXTENSION = '.dat'

if RUNSIMU: assert (not IS1DIM), 'Try 2D methods'

if FIG_CASE == 0:
	if RUNSIMU:
		degList = np.array([1, 2, 3, 4, 5])
		cutList = np.arange(1, 6)
		for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
			sufix = '_' + quadrule + '_' + str(quadtype) + '_' + SUFIX
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			AbserrorTable = np.zeros((len(degList)+1, len(cutList)+1))
			relerrorTable = np.zeros((len(degList)+1, len(cutList)+1))
			AbserrorTable[0, 1:] = cutList; relerrorTable[0, 1:] = cutList
			AbserrorTable[1:, 0] = degList; relerrorTable[1:, 0] = degList
			filenameA1 = FOLDER2DATA+'01L2abserror'+sufix+EXTENSION
			filenameR1 = FOLDER2DATA+'01L2relerror'+sufix+EXTENSION
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

	plotoptions = [CONFIGLINE0, CONFIGLINE1, CONFIGLINE2]
	figname = FOLDER2SAVE+'L2Convergence'+SUFIX+'.pdf'
	if PLOTRELATIVE: filenames = ['0L2relerror_iga_leg_', '0L2relerror_wq_1_', '0L2relerror_wq_2_']
	else: filenames = ['0L2abserror_iga_leg_', '0L2abserror_wq_1_', '0L2abserror_wq_2_']

	fig, ax = plt.subplots(figsize=(5, 5))
	for filename, plotops in zip(filenames, plotoptions):
		quadrule = filename.split('_')[1]
		table = np.loadtxt(FOLDER2DATA+filename+SUFIX+EXTENSION)	
		nbels = 2**(table[0, 1:])
		degList = table[1:, 0]
		errList  = table[1:, 1:]
		for i, degree in enumerate(degList):
			color = COLORLIST[i]
			if quadrule == 'iga': 
				ax.loglog(nbels, errList[i, :], label='ST-IGA-GL $p=$ '+str(int(degree)), color=color, marker=plotops['marker'],
							markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])		
				slope = np.polyfit(np.log10(nbels[2:]),np.log10(errList[i, 2:]), 1)[0]
				slope = round(slope, 1)
				annotation.slope_marker((nbels[-2], errList[i, -2]), slope, 
								poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
			
			else: 
				ax.loglog(nbels, errList[i, :], color=color, marker=plotops['marker'], markerfacecolor='w',
						markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
			fig.savefig(figname)

	ax.loglog([], [], color='k', marker=CONFIGLINE1['marker'], markerfacecolor='w',
			markersize=CONFIGLINE1['markersize'], linestyle=CONFIGLINE1['linestyle'], label='ST-IGA-WQ-1')

	ax.loglog([], [], color='k', marker=CONFIGLINE2['marker'], markerfacecolor='w',
					markersize=CONFIGLINE2['markersize'], linestyle=CONFIGLINE2['linestyle'], label='ST-IGA-WQ-2')
	
	if PLOTRELATIVE:
		ax.set_ylabel('Relative ' + r'$L^2(\Pi)$' + ' error')
		ax.set_ylim(top=1e0, bottom=1e-10)
	else:
		ax.set_ylabel(r'$L^2(\Pi)$' + ' error')
		ax.set_ylim(top=1e1, bottom=1e-7)

	ax.set_xlabel('Number of elements by space-time direction')
	ax.set_xlim(left=1, right=40)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(figname)

elif FIG_CASE == 1:

	filenameA1 = FOLDER2DATA + '1incheatAbs2'+SUFIX
	filenameR1 = FOLDER2DATA + '1incheatRel2'+SUFIX
	filenameT1 = FOLDER2DATA + '1incheatTim2'+SUFIX

	filenameA2 = FOLDER2DATA + '1sptheatAbs2'+SUFIX
	filenameR2 = FOLDER2DATA + '1sptheatRel2'+SUFIX
	filenameT2 = FOLDER2DATA + '1sptheatTim2'+SUFIX

	degsptList = np.arange(1, 3, dtype=int)
	degList = np.array([1, 2, 3, 4, 5])
	cutList = np.arange(1, 7)
	nbel_time = 2**4

	if RUNSIMU:

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
													dirichlet_table=dirichlet_table, 
													geoArgs=geoArgs, 
													quadArgs={'quadrule':'iga', 'type':'leg'},
													degree_time=1, 
													nbel_time=nbel_time, 
													solveSystem=False)[0]
				
				start = time.process_time()
				dirichlet_table = np.ones((2, 2))
				problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensity_inc, 
													dirichlet_table=dirichlet_table, 
													geoArgs=geoArgs,
													nbel_time=nbel_time)
				finish = time.process_time()
				T1timeList[i, j] = finish - start
				
				A1errorList[i, j], R1errorList[i, j] = problem_spt_inc.normOfError(
														np.ravel(temp_inc, order='F'), 
														normArgs={'type':'L2',
																'exactFunction':exactTemperature_spt})
				
				np.savetxt(filenameA1+EXTENSION, A1errorList)
				np.savetxt(filenameR1+EXTENSION, R1errorList)
				np.savetxt(filenameT1+EXTENSION, T1timeList)

				for k, degspt in enumerate(degsptList):
					start = time.process_time()
					dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
					problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, 
														dirichlet_table=dirichlet_table, 
														geoArgs=geoArgs, 
														quadArgs={'quadrule':'wq', 'type':2},
														degree_time=degspt,
														nbel_time=nbel_time)
					finish = time.process_time()
					T2timeList[i, j, k] = finish - start
					
					A2errorList[i, j, k], R2errorList[i, j, k] = problem_spt.normOfError(temp_spt, 
															normArgs={'type':'L2',
																	'exactFunction':exactTemperature_spt,})

					np.save(filenameA2, A2errorList)
					np.save(filenameR2, R2errorList)
					np.save(filenameT2, T2timeList)

	if PLOTRELATIVE:
		filename1 = FOLDER2DATA + '1incheatRel2'+SUFIX
		filename2 = FOLDER2DATA + '1sptheatRel2'+SUFIX
	else:
		filename1 = FOLDER2DATA + '1incheatAbs2'+SUFIX
		filename2 = FOLDER2DATA + '1sptheatAbs2'+SUFIX
	
	errorList1 = np.loadtxt(filename1+EXTENSION)
	errorList2 = np.load(filename2+'.npy')
	fig, ax = plt.subplots(figsize=(5, 5))
	for i, degree in enumerate(degList):
		color = COLORLIST[i]
		ax.loglog(2**cutList, errorList2[i, :, 0], color=color, marker=CONFIGLINE0['marker'], markerfacecolor='w', 
					markersize=CONFIGLINE0['markersize'], linestyle=CONFIGLINE0['linestyle'], 
					label='ST-IGA-GL $p_s=$ '+str(degree))
		
		ax.loglog(2**cutList, errorList1[i, :], color=color, marker=CONFIGLINE4['marker'], markerfacecolor='w',
					markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'])
		
	ax.loglog([], [], color='k', marker=CONFIGLINE4['marker'], markerfacecolor='w',
				markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'], label='INC-IGA-GL '+ r'$\theta=0.5$')

	if PLOTRELATIVE:
		ax.set_ylabel('Relative ' + r'$L^2(\Pi)$' + ' error')
		ax.set_ylim(top=1e0, bottom=1e-6)
	else:
		ax.set_ylabel(r'$L^2(\Pi)$' + ' error')
		ax.set_ylim(top=1e1, bottom=1e-5)

	ax.set_xlabel('Number of elements by spatial direction')
	ax.set_xlim(left=1, right=100)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE + 'IncSptConv1'+ SUFIX + '.pdf')

	if PLOTRELATIVE:
		filename1 = FOLDER2DATA + '1incheatRel'+SUFIX
		filename2 = FOLDER2DATA + '1sptheatRel'+SUFIX
	else:
		filename1 = FOLDER2DATA + '1incheatAbs'+SUFIX
		filename2 = FOLDER2DATA + '1sptheatAbs'+SUFIX

	errorList1 = np.loadtxt(filename1+EXTENSION)
	errorList2 = np.load(filename2+'.npy')
	fig, ax = plt.subplots(figsize=(5, 5))
	for i, degree in enumerate(degList):
		color = COLORLIST[i]
		ax.loglog(2**cutList, errorList2[i, :, 1], color=color, marker=CONFIGLINE0['marker'], markerfacecolor='w', 
					markersize=CONFIGLINE0['markersize'], linestyle=CONFIGLINE0['linestyle'], 
					label='ST-IGA-GL $p_s=$ '+str(degree))
		
		ax.loglog(2**cutList, errorList1[i, :], color=color, marker=CONFIGLINE4['marker'], markerfacecolor='w',
					markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'])
	
	ax.loglog([], [], color='k', marker=CONFIGLINE4['marker'], markerfacecolor='w',
				markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'], label='INC-IGA-GL ' + r'$\theta=0.5$')
	
	if PLOTRELATIVE:
		ax.set_ylabel('Relative ' + r'$L^2(\Pi)$' + ' error')
		ax.set_ylim(top=1e0, bottom=1e-6)
	else:
		ax.set_ylabel(r'$L^2(\Pi)$' + ' error')
		ax.set_ylim(top=1e1, bottom=1e-5)

	ax.set_xlabel('Number of elements by spatial direction')
	ax.set_xlim(left=1, right=100)
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE + 'IncSptConv2'+ SUFIX + '.pdf')

elif FIG_CASE == 2:

	degree, cuts = 8, 6
	quadArgs = {'quadrule':'wq', 'type':2}
	nbelincList = np.array([2**cuts for cuts in range(2, 7)], dtype=int)
	degsptList = np.arange(1, 5)
	abserrorInc, relerrorInc = np.ones(len(nbelincList)), np.ones(len(nbelincList))
	abserrorSpt, relerrorSpt = np.ones((len(degsptList), len(nbelincList))), np.ones((len(degsptList), len(nbelincList)))

	if RUNSIMU:
		for i, nbelinc in enumerate(nbelincList):
			geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
			'nb_refinementByDirection': np.array([cuts, cuts, 1])}

			# dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
			# problem_spt_inc = simulate_spacetime(degree, cuts, powerDensity_spt, 
			# 									dirichlet_table=dirichlet_table, geoArgs=geoArgs, 
			# 									quadArgs={'quadrule':'iga'},
			# 									degree_time=1, nbel_time=nbelinc, solveSystem=False)[0]
			
			# dirichlet_table = np.ones((2, 2))
			# problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensity_inc, dirichlet_table=dirichlet_table,
			# 											geoArgs=geoArgs, nbel_time=nbelinc, quadArgs=quadArgs)
			
			# abserrorInc[i], relerrorInc[i] = problem_spt_inc.normOfError(np.ravel(temp_inc, order='F'), 
			# 							normArgs={'type':'L2',
			# 									'exactFunction':exactTemperature_spt})

			# # abserrorInc[i], relerrorInc[i] = problem_inc.normOfError(temp_inc[:, -1], 
			# # 							normArgs={'type':'L2',
			# # 									'exactFunction':exactTemperature_inc, 
			# # 									'exactExtraArgs':{'time':time_inc[-1]}})
			
			# np.savetxt(FOLDER2DATA+'2abserrorstag_inc'+SUFIX+EXTENSION, abserrorInc)
			# np.savetxt(FOLDER2DATA+'2relerrorstag_inc'+SUFIX+EXTENSION, relerrorInc)

			for j, degspt in enumerate(degsptList):
				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, dirichlet_table=dirichlet_table,
													geoArgs=geoArgs, degree_time=degspt, nbel_time=nbelinc, quadArgs=quadArgs, isfull=True)
					
				abserrorSpt[j, i], relerrorSpt[j, i] = problem_spt.normOfError(temp_spt, 
														normArgs={'type':'L2',
																'exactFunction':exactTemperature_spt,})

				# abserrorSpt[j, i], relerrorSpt[j, i] = problem_inc.normOfError(np.reshape(temp_spt, order='F', 
				# 										newshape=(problem_spt.part.nbctrlpts_total, time_spt.nbctrlpts_total))[:, -1], 
				# 										normArgs={'type':'L2',
				# 												'exactFunction':exactTemperature_inc,
				# 												'exactExtraArgs':{'time':time_inc[-1]}})

				np.savetxt(FOLDER2DATA+'2abserrorstag_spt'+SUFIX+EXTENSION, abserrorSpt)
				np.savetxt(FOLDER2DATA+'2relerrorstag_spt'+SUFIX+EXTENSION, relerrorSpt)

	fig, ax = plt.subplots(figsize=(5,5))

	if PLOTRELATIVE: errorList1 = np.loadtxt(FOLDER2DATA+'2relerrorstag_spt'+SUFIX+EXTENSION)
	else: errorList1 = np.loadtxt(FOLDER2DATA+'2abserrorstag_spt'+SUFIX+EXTENSION)
	for i, deg in enumerate(degsptList):
		nbctrlpts = nbelincList+deg
		ax.loglog(nbctrlpts, errorList1[i, :], color=COLORLIST[i], marker=CONFIGLINE0['marker'], markerfacecolor='w',
				markersize=CONFIGLINE0['markersize'], linestyle=CONFIGLINE0['linestyle'], label='ST-IGA-GL '+r'$p_t=$'+str(int(deg)))
		slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[i, 3:]), 1)[0]
		slope = round(slope, 1)
		annotation.slope_marker((nbctrlpts[-2], errorList1[i, -2]), slope, 
						poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

	if PLOTRELATIVE: errorList1 = np.loadtxt(FOLDER2DATA+'2relerrorstag_inc'+SUFIX+EXTENSION)
	else: errorList1 = np.loadtxt(FOLDER2DATA+'2abserrorstag_inc'+SUFIX+EXTENSION)
	
	nbctrlpts = nbelincList+1
	ax.loglog(nbctrlpts, errorList1, marker=CONFIGLINE5['marker'], markerfacecolor='w', color='k',
					markersize=CONFIGLINE5['markersize'], linestyle=CONFIGLINE5['linestyle'], 
					label='Crank-Nicolson')
	slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[3:]), 1)[0]
	slope = round(slope, 1)
	annotation.slope_marker((nbctrlpts[-2], errorList1[-2]), slope, 
					poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
	
	if PLOTRELATIVE: 
		ax.set_ylabel('Relative '+r'$L^2(\Pi)$' +' error')
		ax.set_ylim(top=1e0, bottom=1e-12)
	else: 
		ax.set_ylabel(r'$L^2(\Pi)$' +' error')
		ax.set_ylim(top=1e1, bottom=1e-8)
	
	ax.set_xlabel('Number of control points in time \n(or number of time-steps)')
	ax.set_xlim(left=2, right=100)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE+'StagnationError'+SUFIX+'.pdf')
