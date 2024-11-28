from article.__init__ import *
from article.input_data import *
from pysrc.lib.lib_base import vtk2png

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
FIG_CASE = 1
EXTENSION = '.dat'

if RUNSIMU: assert (not IS1DIM), 'Try 2D methods'

if FIG_CASE == 0:
	if RUNSIMU:
		degList = np.array([4, 5])
		cutList = np.arange(6, 7)
		# for quadrule, quadtype in zip(['wq', 'wq', 'iga'], [1, 2, 'leg']):
		for quadrule, quadtype in zip(['iga'], ['leg']):
			sufix = '_' + quadrule + '_' + str(quadtype) + '_' + SUFIX
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			AbserrorTable = np.zeros((len(degList)+1, len(cutList)+1))
			relerrorTable = np.zeros((len(degList)+1, len(cutList)+1))
			AbserrorTable[0, 1:] = cutList; relerrorTable[0, 1:] = cutList
			AbserrorTable[1:, 0] = degList; relerrorTable[1:, 0] = degList
			filenameA1 = FOLDER2DATA+'0L2abserror'+sufix+EXTENSION
			filenameR1 = FOLDER2DATA+'0L2relerror'+sufix+EXTENSION
			for j, cuts in enumerate(cutList):
				for i, degree in enumerate(degList):
					geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
							'nb_refinementByDirection': np.array([cuts, cuts, 1])}
					problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, quadArgs=quadArgs, 
																		degree_time=degree, nbel_time=2**cuts, geoArgs=geoArgs)
					AbserrorTable[i+1, j+1], relerrorTable[i+1, j+1] = problem_spt.normOfError(temp_spt, 
																		normArgs={'type':'L2',
																				'exactFunction':exactTemperature_spt},)

					print(relerrorTable[i+1, j+1])
					# np.savetxt(filenameA1, AbserrorTable)
					# np.savetxt(filenameR1, relerrorTable)

	plotoptions = [CONFIGLINE0, CONFIGLINE1, CONFIGLINE2]
	figname = FOLDER2SAVE+'L2Convergence'+SUFIX+'.pdf'
	if PLOTRELATIVE: filenames = ['0L2relerror_iga_leg_', '0L2relerror_wq_1_', '0L2relerror_wq_2_']
	else: filenames = ['0L2abserror_iga_leg_', '0L2abserror_wq_1_', '0L2abserror_wq_2_']

	fig, ax = plt.subplots(figsize=(5, 5))
	for filename, plotops in zip(filenames, plotoptions):
		quadrule = filename.split('_')[1]
		table = np.loadtxt(FOLDER2DATA+filename+SUFIX+EXTENSION)	
		nbels = 2**(table[0, 1:])
		# degList = table[1:-1, 0]
		# errList  = table[1:, 1:]
		degList = table[1:, 0]
		errList  = table[1:, 1:]
		for i, degree in enumerate(degList):
			color = COLORLIST[i]
			if quadrule == 'iga': 
				ax.loglog(nbels, errList[i, :], label='ST-IGA-GL $p=$ '+str(int(degree)), color=color, marker=plotops['marker'],
							markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])		
				# slope = np.polyfit(np.log10(nbels[2:]),np.log10(errList[i, 2:]), 1)[0]
				# slope = round(slope, 1)
				# annotation.slope_marker((nbels[-2], errList[i, -2]), slope, 
				# 				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
			
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
		ax.set_ylim(top=1e2, bottom=1e-8)
	else:
		ax.set_ylabel(r'$L^2(\Pi)$' + ' error')
		ax.set_ylim(top=1e1, bottom=1e-7)

	ax.set_xlabel('Number of elements by space-time direction')
	ax.set_xlim(left=1, right=100)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(figname)


	# ================

	for i, [filename, plotops] in enumerate(zip(filenames, plotoptions)):
		fig, ax = plt.subplots(figsize=(5, 5))
		figname = FOLDER2SAVE+'L2Convergence'+str(i)+SUFIX+'.pdf'
	
		quadrule = filename.split('_')[1]
		table = np.loadtxt(FOLDER2DATA+filename+SUFIX+EXTENSION)	
		nbels = 2**(table[0, 1:])

		degList = table[1:, 0]
		errList = table[1:, 1:]
		
		for i, degree in enumerate(degList):
			color = COLORLIST[i]
			tmpyy = errList[i, :]
			yy = tmpyy[tmpyy>0]
			xx =nbels[tmpyy>0]

			if quadrule == 'iga': 
				ax.loglog(xx, yy, label='ST-IGA-GL $p=$ '+str(int(degree)), color=color, marker=plotops['marker'],
							markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])		
				if degree < 4:
					slope = np.polyfit(np.log10(xx[-2:]),np.log10(yy[-2:]), 1)[0]
					slope = round(slope, 1)
					annotation.slope_marker((xx[-2], yy[-2]), slope, 
									poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
			
			else: 
				ax.loglog(xx, yy, color=color, marker=plotops['marker'], markerfacecolor='w',
						markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
			fig.savefig(figname)

		# ax.loglog([], [], color='k', marker=CONFIGLINE1['marker'], markerfacecolor='w',
		# 		markersize=CONFIGLINE1['markersize'], linestyle=CONFIGLINE1['linestyle'], label='ST-IGA-WQ-1')

		# ax.loglog([], [], color='k', marker=CONFIGLINE2['marker'], markerfacecolor='w',
		# 				markersize=CONFIGLINE2['markersize'], linestyle=CONFIGLINE2['linestyle'], label='ST-IGA-WQ-2')
		
		if PLOTRELATIVE:
			ax.set_ylabel('Relative ' + r'$L^2(\Pi)$' + ' error')
			ax.set_ylim(top=1e1, bottom=1e-7)
		else:
			ax.set_ylabel(r'$L^2(\Pi)$' + ' error')
			ax.set_ylim(top=1e1, bottom=1e-7)

		ax.set_xlabel('Number of elements by space-time direction')
		ax.set_xlim(left=1, right=100)
		if quadrule == 'iga': ax.legend(loc='lower left')
		fig.tight_layout()
		fig.savefig(figname)

elif FIG_CASE == 1:

	degree, cuts = 3, 4
	subfolder = FOLDER2DATA + GEONAME + '_' + str(degree) + '_' + str(cuts) + '/' 
	if not os.path.isdir(subfolder): os.mkdir(subfolder)

	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, cuts, 1])}

	if RUNSIMU:
		for [i, isadaptive], prefix1 in zip(enumerate([False, True]), ['exact', 'inexact']):
			for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
				prefix = prefix1 + '_' + prefix2 + '_'

				problem_spt, _, _, output = simulate_spacetime(degree, cuts, powerDensity_spt, geoArgs=geoArgs,
																degree_time=degree, nbel_time=2**cuts, 
																isadaptive=isadaptive, isfull=isfull,
																getOthers=True)
				
				solution_list  = output['Solution']; resLin_list = output['KrylovRes']
				resNonLin_list = output['NewtonRes']; adapThres_list = output['Threshold']
				L2error_list, L2relerror_list = [], []

				for solution in solution_list:
					L2error, L2relerror  = problem_spt.normOfError(solution, normArgs={'type':'L2', 
																	'exactFunction':exactTemperature_spt})
					L2error_list.append(L2error); L2relerror_list.append(L2relerror)
				resLinclean = np.array([]); counter_list = [0]
				for _ in resLin_list: 
					resLinclean = np.append(resLinclean, _[np.nonzero(_)])
					counter_list.append(counter_list[-1] + len(_[np.nonzero(_)]))

				np.savetxt(subfolder+prefix+'CumulKrylovRes'+EXTENSION, resLinclean)
				np.savetxt(subfolder+prefix+'Inner_loops'+EXTENSION, counter_list)
				np.savetxt(subfolder+prefix+'NewtonRes'+EXTENSION, resNonLin_list)
				np.savetxt(subfolder+prefix+'L2error'+EXTENSION, L2error_list)
				np.savetxt(subfolder+prefix+'L2relerror'+EXTENSION, L2relerror_list)
				np.savetxt(subfolder+prefix+'Threshold'+EXTENSION, np.array(adapThres_list))

	fig1, ax1 = plt.subplots(figsize=(5,4))
	fig2, ax2 = plt.subplots(figsize=(5,4))
	fig3, ax3 = plt.subplots(figsize=(5,4))
	fig4, ax4 = plt.subplots(figsize=(5,4))
	figs = [fig1, fig2, fig3, fig4]; axs  = [ax1, ax2, ax3, ax4]
	linestyle_list = ['-', '--', '-', '--']
	marker_list = ['o', 'o', 's', 's']

	for [i, isadaptive], prefix1 in zip(enumerate([True, False]), ['inexact', 'exact']):
		for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
			l = j + i*2
			legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
			prefix = prefix1 + '_' + prefix2 + '_'
			nbInnerLoops = np.loadtxt(subfolder+prefix+'Inner_loops'+EXTENSION)
			newtonRes = np.loadtxt(subfolder+prefix+'NewtonRes'+EXTENSION)
			L2relerror = np.loadtxt(subfolder+prefix+'L2relerror'+EXTENSION)
			newtonRes = newtonRes/newtonRes[0]
			
			ylim1 = np.power(10, np.floor(np.log10(np.min(L2relerror))))
			ylim2 = np.power(10, np.floor(np.log10(np.min(newtonRes))))
			# ylim = np.min([ylim1, ylim2])
			for caseplot, fig, ax in zip(range(1, 5), figs, axs):
				ylim = 1e-11
				if caseplot == 1:
					yy = L2relerror; xx = nbInnerLoops[:len(L2relerror)]
					xlim = 800; ylim=1e-5
					ylabel = 'Relative '+r'$L^2(\Pi)$' + ' norm of error'
					xlabel = 'Number of matrix-vector products'
				elif caseplot == 2:
					yy = newtonRes; xx = nbInnerLoops[:len(newtonRes)]
					xlim = 800
					ylabel = 'Relative norm of nonlinear residue'
					xlabel = 'Number of matrix-vector products'
				elif caseplot == 3:
					yy = newtonRes; xx = np.arange(0, len(newtonRes)) 
					xlim = 20
					ylabel = 'Relative norm of nonlinear residue'
					xlabel = 'Number of nonlinear iterations'
				elif caseplot == 4:
					yy = L2relerror; xx = np.arange(0, len(L2relerror))
					xlim = 20; ylim=1e-5
					ylabel = 'Relative '+r'$L^2(\Pi)$' + ' norm of error'
					xlabel = 'Number of nonlinear iterations'

				ax.semilogy(xx, yy, label=legendname, marker=marker_list[l], linestyle=linestyle_list[l])
				ax.set_xlim(right=xlim, left=1e0)
				ax.set_ylim(top=1e1, bottom=ylim)
				ax.set_xlabel(xlabel)
				ax.set_ylabel(ylabel)
				if caseplot==3 or caseplot==4: ax.legend()
						
				fig.tight_layout()
				fig.savefig(FOLDER2SAVE+'NLConvergence_iters'+GEONAME+'_'+str(degree)+str(cuts)+str(caseplot)+'.pdf')

	fig, ax = plt.subplots()
	for [i, isadaptive], prefix1 in zip(enumerate([True]), ['inexact']):
		for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
			l = j + i*2
			legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
			prefix = prefix1 + '_' + prefix2 + '_'
			threshold = np.loadtxt(subfolder+prefix+'Threshold'+EXTENSION)

			yy = threshold; xx = np.arange(0, len(threshold))
			xlim = 20
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
			fig.savefig(FOLDER2SAVE+'NLTolerance'+GEONAME+'_'+str(degree)+str(cuts)+'.pdf')

elif FIG_CASE == 2:
	degree, cuts = 6, 4
	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, cuts, 1])}
	problem_spt = simulate_spacetime(degree, cuts, None, geoArgs=geoArgs,
									degree_time=degree, nbel_time=2**cuts, 
									solveSystem=False)[0]
	
	problem_spt.part.postProcessingPrimal(fields={'temp':exactTemperature_inc, 
													'conductivity':nonlinearfunc2,
													'dersconductivity':nonlineardersfunc2,
												}, 
									folder=FOLDER2SAVE, name='spt_'+GEONAME, 
									extraArgs={'time':1.0, 'temperature':exactTemperature_inc})
	vtk2png(folder=FOLDER2SAVE, filename='spt_'+GEONAME, fieldname='temp', cmap='coolwarm', title='Temperature', position_y=0.25)
	vtk2png(folder=FOLDER2SAVE, filename='spt_'+GEONAME, fieldname='conductivity', cmap='viridis', title='Conductivity', position_y=0.25)
	vtk2png(folder=FOLDER2SAVE, filename='spt_'+GEONAME, fieldname='dersconductivity', cmap='viridis', title=' Ders conductivity', position_y=0.25, fmt="%.1e",)