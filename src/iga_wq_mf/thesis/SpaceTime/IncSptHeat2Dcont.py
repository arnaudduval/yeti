from thesis.SpaceTime.__init__ import *
from thesis.SpaceTime.input_data import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
FIG_CASE = 5
EXTENSION = '.dat'

if RUNSIMU: assert (not IS1DIM), 'Try 2D methods'

if FIG_CASE == 3:

	filenameA2 = FOLDER2DATA + '3incheatAbs'
	filenameR2 = FOLDER2DATA + '3incheatRel'
	filenameT2 = FOLDER2DATA + '3incheatTim'

	filenameA3 = FOLDER2DATA + '3sptheatAbs'
	filenameR3 = FOLDER2DATA + '3sptheatRel'
	filenameT3 = FOLDER2DATA + '3sptheatTim'

	degList = np.array([1, 2, 3, 4, 5, 6])
	cutList = np.arange(4, 7)

	if RUNSIMU:
		A2errorList = np.ones((len(degList), len(cutList)))
		R2errorList = np.ones((len(degList), len(cutList)))
		T2timeList = np.ones((len(degList), len(cutList)))
		
		quadArgs = {'quadrule': 'wq', 'type': 2}
		sufix = '_wq_2_' + SUFIX + EXTENSION
		for j, cuts in enumerate(cutList):
			for i, degree in enumerate(degList):
				geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, cuts, 1])}

				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensity_spt, 
													dirichlet_table=dirichlet_table, 
													geoArgs=geoArgs, 
													degree_time=1, 
													nbel_time=2**cuts, 
													quadArgs=quadArgs,
													isadaptive=False, solveSystem=False)

				start = time.process_time()
				dirichlet_table = np.ones((2, 2))
				problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensity_inc, 
													dirichlet_table=dirichlet_table, geoArgs=geoArgs, nbel_time=2**cuts)
				finish = time.process_time()
				T2timeList[i, j] = finish - start
				
				A2errorList[i, j], R2errorList[i, j] = problem_spt.normOfError(np.ravel(temp_inc, order='F'), 
														normArgs={'type':'L2',
																'exactFunction':exactTemperature_spt,})

				np.savetxt(filenameA2+sufix, A2errorList)
				np.savetxt(filenameR2+sufix, R2errorList)
				np.savetxt(filenameT2+sufix, T2timeList)

		A3errorList = np.ones((len(degList), len(cutList)))
		R3errorList = np.ones((len(degList), len(cutList)))
		T3timeList = np.ones((len(degList), len(cutList)))
		
		for quadrule, quadtype in zip(['wq', 'wq', 'iga'], [2, 1, 'leg']):
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			sufix = '_' + quadrule + '_' + str(quadtype) + '_' + SUFIX + EXTENSION
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

	position = 1
	assert position in [1, 2], 'Must be one or 2'
	if position==1: fig, ax = plt.subplots(figsize=(6, 4))
	if position==2: fig, ax = plt.subplots(figsize=(5.5, 4))
	cmap = mpl.colors.ListedColormap(COLORLIST[:len(degList)])

	if PLOTRELATIVE:
		filenameA2 = FOLDER2DATA + '3incheatRel'
		filenameA3 = FOLDER2DATA + '3sptheatRel'

	Elist = np.loadtxt(filenameA2+'_wq_2_'+SUFIX+EXTENSION)
	Tlist = np.loadtxt(filenameT2+'_wq_2_'+SUFIX+EXTENSION)
	im = ax.scatter(Tlist[:len(degList), position], Elist[:len(degList), position], c=degList,
			cmap=cmap, marker=CONFIGLINE4['marker'], s=10*CONFIGLINE4['markersize'])
	ax.loglog(Tlist[:len(degList), position], Elist[:len(degList), position], 
			color='k', marker='', linestyle=CONFIGLINE4['linestyle'])
	
	if position==1:
		cbar = fig.colorbar(im); cbar.set_label('Degree')
		tick_locs = 1+(np.arange(len(degList)) + 0.5)*(len(degList)-1)/len(degList)
		cbar.set_ticks(tick_locs)
		cbar.set_ticklabels(degList)

	for quadrule, quadtype, plotvars in zip(['iga', 'wq'], ['leg', 2], [CONFIGLINE0, CONFIGLINE2]):
		sufix = '_' + quadrule + '_' + str(quadtype) + '_' + SUFIX + EXTENSION
		Elist = np.loadtxt(filenameA3+sufix)
		Tlist = np.loadtxt(filenameT3+sufix)
		ax.scatter(Tlist[:len(degList), position], Elist[:len(degList), position], c=degList,
						cmap=cmap, marker=plotvars['marker'], s=10*plotvars['markersize'])
			
		ax.loglog(Tlist[:len(degList), position], Elist[:len(degList), position], 
				color='k', marker='', linestyle=plotvars['linestyle'])

	ax.loglog([], [], color='k', marker=CONFIGLINE4['marker'], alpha=0.5,
			markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'], label='INC-IGA-WQ')
	
	ax.loglog([], [], color='k', marker=CONFIGLINE0['marker'], alpha=0.5,
		markersize=CONFIGLINE0['markersize'], linestyle=CONFIGLINE0['linestyle'], label='ST-IGA-GL')
		
	ax.loglog([], [], color='k', marker=CONFIGLINE2['marker'], alpha=0.5,
		markersize=CONFIGLINE2['markersize'], linestyle=CONFIGLINE2['linestyle'], label='ST-IGA-WQ 2')

	if PLOTRELATIVE: 
		ax.set_ylabel('Relative ' + r'$L^2$' + ' error')
		ax.set_ylim(top=1e-2, bottom=1e-12)
	else:
		ax.set_ylabel(r'$L^2$' + ' error')
		ax.set_ylim(top=1e-1, bottom=1e-11)

	ax.set_xlabel('CPU time (s)')
	if position==1: ax.legend(loc='lower left')
	if position==1: ax.set_xlim(left=1e0, right=1e3)
	if position==2: ax.set_xlim(left=1e1, right=1e4)
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE + 'SPTINC_CPUError0' + str(position) +  '.pdf')

elif FIG_CASE == 4:

	degree, cuts = 3, 4
	subfolderfolder = FOLDER2DATA + str(degree) + '_' + str(cuts) + '/' 
	if not os.path.isdir(subfolderfolder): os.mkdir(subfolderfolder)

	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, cuts, 1])}

	if RUNSIMU:
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

				np.savetxt(subfolderfolder+prefix+'CumulKrylovRes'+EXTENSION, resKrylovclean)
				np.savetxt(subfolderfolder+prefix+'Inner_loops'+EXTENSION, counter_list)
				np.savetxt(subfolderfolder+prefix+'NewtonRes'+EXTENSION, resNewton_list)
				np.savetxt(subfolderfolder+prefix+'L2error'+EXTENSION, L2error)
				np.savetxt(subfolderfolder+prefix+'L2relerror'+EXTENSION, L2relerror)
				np.savetxt(subfolderfolder+prefix+'Threshold'+EXTENSION, np.array(threshold_list))

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
			nbInnerLoops = np.loadtxt(subfolderfolder+prefix+'Inner_loops'+EXTENSION)
			newtonRes = np.loadtxt(subfolderfolder+prefix+'NewtonRes'+EXTENSION)
			L2relerror = np.loadtxt(subfolderfolder+prefix+'L2relerror'+EXTENSION)
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
				if caseplot==3: ax.legend()
				fig.tight_layout()
				fig.savefig(FOLDER2SAVE+'NLConvergence_iters'+'_'+str(degree)+str(cuts)+str(caseplot)+'.pdf')

	fig, ax = plt.subplots(figsize=(6, 4))
	for [i, isadaptive], prefix1 in zip(enumerate([True]), ['inexact']):
		for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
			l = j + i*2
			legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
			prefix = prefix1 + '_' + prefix2 + '_'
			threshold = np.loadtxt(subfolderfolder+prefix+'Threshold'+EXTENSION)

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
			fig.savefig(FOLDER2SAVE+'NLTolerance'+'_'+str(degree)+str(cuts)+'.pdf')

elif FIG_CASE == 5:
	
	degList = np.array([1, 2, 3, 4, 5, 6])
	cutList = np.arange(4, 7)

	fig, axs = plt.subplots(1, 2, figsize=(8, 3.5))
	cmap = mpl.colors.ListedColormap(COLORLIST[:len(degList)])
	filenameA3 = FOLDER2DATA + '3sptheatRel'
	filenameT3 = FOLDER2DATA + '3sptheatTim'

	for quadrule, quadtype, plotvars, ax in zip(['iga', 'wq'], ['leg', 2], [CONFIGLINE0, CONFIGLINE2], [axs[0], axs[1]]):
		sufix = '_' + quadrule + '_' + str(quadtype) + '_' + SUFIX + EXTENSION
		Elist = np.loadtxt(filenameA3+sufix)
		Tlist = np.loadtxt(filenameT3+sufix)
		for pos in range(np.size(Elist, axis=1)):
			im = ax.scatter(Tlist[:len(degList), pos], Elist[:len(degList), pos], c=degList,
							cmap=cmap, marker=plotvars['marker'], s=10*plotvars['markersize'])
				
			ax.loglog(Tlist[:len(degList), pos], Elist[:len(degList), pos], 
					color='k', marker='', linestyle=plotvars['linestyle'])
			ax.text(Tlist[-1, pos]*1.2, Elist[-1, pos]/5, str(int(2**(pos+4)))+r'$^3$')

	divider1 = make_axes_locatable(axs[0])
	cax1 = divider1.append_axes("right", size="5%", pad=0.1)
	divider2 = make_axes_locatable(axs[1])
	cax2 = divider2.append_axes("right", size="5%", pad=0.1)

	cbar = plt.colorbar(im, cax=cax1)
	fig.delaxes(fig.axes[2])

	cbar = plt.colorbar(im, cax=cax2)
	cbar.set_label('Degree')
	tick_locs = 1+(np.arange(len(degList)) + 0.5)*(len(degList)-1)/len(degList)
	cbar.set_ticks(tick_locs)
	cbar.set_ticklabels(degList)

	axs[0].set_ylabel('Relative ' + r'$L^2$' + ' error')
	axs[0].set_ylim(top=1e-1, bottom=1e-12)
	axs[1].set_ylim(top=1e-1, bottom=1e-12)
	axs[0].set_xlim(left=5e-1, right=5e4)
	axs[1].set_xlim(left=5e-1, right=5e4)
	axs[0].title.set_text('Gauss-Legendre')
	axs[1].title.set_text('Weighted quadrature 2')

	axs[0].set_xlabel('CPU time (s)')
	axs[1].set_xlabel('CPU time (s)')
	plt.tight_layout()
	fig.savefig(FOLDER2SAVE + 'SPTINC_CPUError' +  '.pdf')