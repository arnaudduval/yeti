from pysrc.sptmpaper.input_data import *
from pyevtk.vtk import VtkGroup

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/Other/' + GEONAME + '/'
subfolder = folder +  'steps/'
if not os.path.isdir(folder): os.mkdir(folder)
if not os.path.isdir(subfolder): os.mkdir(subfolder)

def run(folder=None):
	assert folder is not None, 'Folder unknown'
	print("Running group...")
	g = VtkGroup(folder)
	for i in range(33):
		g.addFile(filepath = folder + "out_"+str(i)+".vts", sim_time = i)
	g.save()

# Set global variables
TODOSIMU = True
FIG_CASE = 2

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
		degree_list = np.array([1, 2, 3, 4, 5])
		cuts_list   = np.arange(1, 6)
		for quadrule, quadtype in zip(['iga', 'wq'], ['leg', 1]):
			sufix = '_' + quadrule + '_' + str(quadtype) + '_' + lastsufix
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			L2errorTable = np.zeros((len(degree_list)+1, len(cuts_list)+1))
			L2relerrorTable = np.zeros((len(degree_list)+1, len(cuts_list)+1))
			L2errorTable[0, 1:] = cuts_list; L2relerrorTable[0, 1:] = cuts_list
			L2errorTable[1:, 0] = degree_list; L2relerrorTable[1:, 0] = degree_list
			filenameA1 = folder+'L2error_meshpar'+sufix+'.dat'
			filenameR1 = folder+'L2relerror_meshpar'+sufix+'.dat'
			for j, cuts in enumerate(cuts_list):
				for i, degree in enumerate(degree_list):
					nbels = 2**cuts_list
					geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
							'nb_refinementByDirection': np.array([cuts, cuts, 1])}
					blockPrint()
					dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
					problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, 
													powerdensity=powerDensity_spt, 
													dirichlet_table=dirichlet_table,
													quadArgs=quadArgs, degree_spt=degree, cuts_spt=cuts, geoArgs=geoArgs)
					
					enablePrint()
					L2errorTable[i+1, j+1], L2relerrorTable[i+1, j+1] = problem_spt.normOfError(temp_spt, 
																	normArgs={'type':'L2', 
																	'exactFunction':exactTemperature_spt},)

					np.savetxt(filenameA1, L2errorTable)
					np.savetxt(filenameR1, L2relerrorTable)

	normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
	onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
	onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}
	plotoptions = [normalPlot, onlyMarker1, onlyMarker2]

	lastsufix = 'linear' if ISLINEAR else 'nonlin'
	figname = folder + 'SPTNonLinearConvergenceL2'+lastsufix+'.pdf'
	filenames = ['L2error_meshpar_iga_leg_', 'L2error_meshpar_wq_1_']

	normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
	fig, ax = plt.subplots(figsize=(8, 6))

	for filename, plotops in zip(filenames, plotoptions):
		quadrule = filename.split('_')[2]
		table = np.loadtxt(folder+filename+lastsufix+'.dat')	
		nbels = 2**(table[0, 1:])
		degrees = table[1:, 0]
		errors  = table[1:, 1:]
		for i, degree in enumerate(degrees):
			color = COLORLIST[i]
			if quadrule == 'iga': 
				ax.loglog(nbels, errors[i, :], label='IGA-GL deg. '+str(int(degree)), color=color, marker=plotops['marker'],
							markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])
				
				# slope = np.polyfit(np.log10(nbels[2:]),np.log10(errors[i, 2:]), 1)[0]
				# slope = round(slope, 1)
				# annotation.slope_marker((nbels[-1], errors[i, -1]), slope, 
				# 				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)			
			else: 
				ax.loglog(nbels, errors[i, :], color=color, marker=plotops['marker'], markerfacecolor='w',
						markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
			fig.savefig(figname)

	ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
					markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ")
	# ax.loglog([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
	# 		markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label="IGA-WQ 4")

	# ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$')
	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}$')
	ax.set_xlabel('Number of elements by space-time direction')
	ax.set_xlim(left=1, right=100)
	ax.set_ylim(top=1e1, bottom=1e-10)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(figname)

elif FIG_CASE == 1:

	filenameA1 = folder + 'incheatAbs'+lastsufix
	filenameR1 = folder + 'incheatRel'+lastsufix
	filenameT1 = folder + 'incheatTim'+lastsufix

	filenameA2 = folder + 'sptheatAbs'+lastsufix
	filenameR2 = folder + 'sptheatRel'+lastsufix
	filenameT2 = folder + 'sptheatTim'+lastsufix

	filenameA3 = folder + 'sptheatAbs2'+lastsufix
	filenameR3 = folder + 'sptheatRel2'+lastsufix
	filenameT3 = folder + 'sptheatTim2'+lastsufix

	degree_list = np.array([1, 2, 3, 4, 5])
	cuts_list   = np.arange(1, 7)

	if TODOSIMU:
		exportTimeDependentMaterial(np.linspace(0, 1, 33), 
						temperature=exactTemperature_inc,
						fields={'mat':nonlinearfunc, 'temp':exactTemperature_inc}, 
						geoArgs={'name': GEONAME, 'degree': 5*np.ones(3, dtype=int), 
								'nb_refinementByDirection': 5*np.ones(3, dtype=int)},
						folder=subfolder,)
		run(folder=subfolder)

		Aerror_list1 = np.ones((len(degree_list), len(cuts_list)))
		Rerror_list1 = np.ones((len(degree_list), len(cuts_list)))
		time_list1 = np.ones((len(degree_list), len(cuts_list)))

		Aerror_list2 = np.ones((len(degree_list), len(cuts_list)))
		Rerror_list2 = np.ones((len(degree_list), len(cuts_list)))
		time_list2 = np.ones((len(degree_list), len(cuts_list)))

		Aerror_list3 = np.ones((len(degree_list), len(cuts_list)))
		Rerror_list3 = np.ones((len(degree_list), len(cuts_list)))
		time_list3 = np.ones((len(degree_list), len(cuts_list)))
		
		for j, cuts in enumerate(cuts_list):
			for i, degree in enumerate(degree_list):
				geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, cuts, 1])}

				# Incremental problem
				blockPrint()
				start = time.process_time()
				dirichlet_table = np.ones((2, 2))
				problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, dirichlet_table=dirichlet_table,
												powerdensity=powerDensity_inc, geoArgs=geoArgs)
				end = time.process_time()
				time_list1[i, j] = end - start
				
				# Error of last step
				Aerror_list1[i, j], Rerror_list1[i, j] = problem_inc.normOfError(temp_inc[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperature_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})
				# Space time problem
				start = time.process_time()
				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, dirichlet_table=dirichlet_table,
													powerdensity=powerDensity_spt, geoArgs=geoArgs, degree_spt=degree,)
					
				end = time.process_time()
				time_list2[i, j] = end - start

				# Error of last "step"
				newtemp_spt = np.reshape(temp_spt, newshape=(problem_spt.part.nbctrlpts_total, problem_spt.time.nbctrlpts_total), order='F')
				Aerror_list2[i, j], Rerror_list2[i, j] = problem_inc.normOfError(newtemp_spt[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperature_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})
				
				# Space time problem
				start = time.process_time()
				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, dirichlet_table=dirichlet_table,
													powerdensity=powerDensity_spt, geoArgs=geoArgs, degree_spt=degree, 
													quadArgs={'quadrule':'wq', 'type':1})
					
				end = time.process_time()
				time_list3[i, j] = end - start

				# Error of last "step"
				newtemp_spt = np.reshape(temp_spt, newshape=(problem_spt.part.nbctrlpts_total, problem_spt.time.nbctrlpts_total), order='F')
				Aerror_list3[i, j], Rerror_list3[i, j] = problem_inc.normOfError(newtemp_spt[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperature_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})

				enablePrint()
				print('inc: %.3e, spt: %.3e, spt2: %.3e' %(time_list1[i, j], time_list2[i, j], time_list3[i, j]))
				print('inc: %.3e, spt: %.3e, spt2: %.3e' %(Aerror_list1[i, j], Aerror_list2[i, j], Aerror_list3[i, j], ))
				print('...')
			np.save(filenameA1, Aerror_list1)
			np.save(filenameR1, Rerror_list1)
			np.save(filenameT1, time_list1)
			np.save(filenameA2, Aerror_list2)
			np.save(filenameR2, Rerror_list2)
			np.save(filenameT2, time_list2)
			np.save(filenameA3, Aerror_list3)
			np.save(filenameR3, Rerror_list3)
			np.save(filenameT3, time_list3)

			print('------')

	list1 = np.load(filenameA1+'.npy')
	list2 = np.load(filenameA2+'.npy')
	list3 = np.load(filenameA3+'.npy')
	fig, ax = plt.subplots(figsize=(8, 6))
	for i, degree in enumerate(degree_list):
		if degree%2==0: continue
		color = COLORLIST[i]
		ax.loglog(2**cuts_list, list2[i, :], color=color, marker='s', markerfacecolor='w',
					markersize=10, linestyle='-', label='ST-IGA-GL deg. ' + str(degree))
		
		ax.loglog(2**cuts_list, list3[i, :], color=color, marker='o', markerfacecolor='w',
					markersize=8, linestyle='--')

		ax.loglog(2**cuts_list, list1[i, :], color=color, marker='o',
					markersize=4, linestyle=':')
		
	ax.loglog([], [], color='k', marker='o', markerfacecolor='w',
					markersize=8, linestyle='--', label="ST-IGA-WQ")
	ax.loglog([], [], color='k', marker='o', 
					markersize=4, linestyle=':', label="INC-IGA-GL")
	
	# ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}$')
	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Omega)}$')
	ax.set_xlabel('Number of elements by spatial direction')
	ax.set_xlim(left=1, right=100)
	ax.set_ylim(top=1e1, bottom=1e-10)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(folder + 'FigConvergSptIncHeat' +  '.pdf')
	plt.close(fig)

	# -------------------

	list1 = np.load(filenameT1+'.npy')
	list2 = np.load(filenameT2+'.npy')
	list3 = np.load(filenameT3+'.npy')
	fig, ax = plt.subplots(figsize=(8, 6))
	for i, degree in enumerate(degree_list):
		if degree%2==0: continue
		color = COLORLIST[i]
		ax.loglog(2**cuts_list, list2[i, :], color=color, marker='s', markerfacecolor='w',
					markersize=10, linestyle='-', label='ST-IGA-GL deg. ' + str(degree))
		
		ax.loglog(2**cuts_list, list3[i, :], color=color, marker='o', markerfacecolor='w',
					markersize=8, linestyle='--')

		ax.loglog(2**cuts_list, list1[i, :], color=color, marker='o',
					markersize=4, linestyle=':')
		
	ax.loglog([], [], color='k', marker='o', markerfacecolor='w',
					markersize=8, linestyle='--', label="ST-IGA-WQ")
	ax.loglog([], [], color='k', marker='o', 
					markersize=4, linestyle=':', label="INC-IGA-GL")
	
	ax.set_ylabel('CPU time (s)')
	ax.set_xlabel('Number of elements by spatial direction')
	ax.set_xlim(left=1, right=100)
	ax.set_ylim(top=5e2, bottom=5e-3)
	ax.legend(loc='upper left')
	fig.tight_layout()
	fig.savefig(folder + 'FigTimeSptIncHeat' +  '.pdf')
	plt.close(fig)

	# -------------------
	Elist1 = np.load(filenameA1+'.npy')
	Elist2 = np.load(filenameA2+'.npy')
	Elist3 = np.load(filenameA3+'.npy')
	Tlist1 = np.load(filenameT1+'.npy')
	Tlist2 = np.load(filenameT2+'.npy')
	Tlist3 = np.load(filenameT3+'.npy')
	fig, ax = plt.subplots(figsize=(8, 6))
	for i, degree in enumerate(degree_list):
		if degree%2==0: continue
		color = COLORLIST[i]
		ax.loglog(Tlist2[i, :], Elist2[i, :], color=color, marker='s', markerfacecolor='w',
					markersize=10, linestyle='-', label='ST-IGA-GL deg. ' + str(degree))
		
		ax.loglog(Tlist3[i, :], Elist3[i, :], color=color, marker='o', markerfacecolor='w',
					markersize=8, linestyle='--')

		ax.loglog(Tlist1[i, :], Elist1[i, :], color=color, marker='o',
					markersize=4, linestyle=':')
		
	ax.loglog([], [], color='k', marker='o', markerfacecolor='w',
					markersize=8, linestyle='--', label="ST-IGA-WQ")
	ax.loglog([], [], color='k', marker='o', 
					markersize=4, linestyle=':', label="INC-IGA-GL")
	
	ax.set_xlabel('CPU time (s)')
	ax.set_ylabel(r'$L^2(\Omega)$'+' error')
	ax.set_xlim(left=5e-3, right=5e2)
	ax.set_ylim(top=1e1, bottom=1e-10)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(folder + 'FigProfitSptIncHeat' +  '.pdf')
	plt.close(fig)

elif FIG_CASE == 2:
	extension = '.dat'
	degree, cuts = 4, 4
	geoArgs = {'name': GEONAME, 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, cuts, 1])}

	subfolderfolder = folder + '_' + str(degree) + '_' + str(cuts) + '/' 
	if not os.path.isdir(subfolderfolder): os.mkdir(subfolderfolder)

	if TODOSIMU:
		for [i, isadaptive], prefix1 in zip(enumerate([False, True]), ['exact', 'inexact']):
			for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
				prefix = prefix1 + '_' + prefix2 + '_'
				start=time.time()
				# blockPrint()
				dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, _, _, output = simulate_spacetime(degree, cuts, 
												powerdensity=powerDensity_spt, dirichlet_table=dirichlet_table,
												degree_spt=degree, geoArgs=geoArgs,
												isadaptive=isadaptive, isfull=isfull,
												getOthers=True)
				
				sol_list  = output['Solution']; resKrylov_list = output['KrylovRes']
				resNewton_list = output['NewtonRes']; threshold_list = output['Threshold']
				L2error, L2relerror = [], []
				stop = time.time()

				for solution in sol_list:
					err, relerr  = problem_spt.normOfError(solution, normArgs={'type':'L2', 
																	'exactFunction':exactTemperature_spt})
					L2error.append(err); L2relerror.append(relerr)
				resKrylovclean = np.array([]); counter_list = [0]
				for _ in resKrylov_list: 
					resKrylovclean = np.append(resKrylovclean, _[np.nonzero(_)])
					counter_list.append(counter_list[-1] + len(_[np.nonzero(_)]))
				enablePrint()
				print(err, stop-start)

				np.savetxt(subfolderfolder+prefix+'CumulKrylovRes'+extension, resKrylovclean)
				np.savetxt(subfolderfolder+prefix+'Inner_loops'+extension, counter_list)
				np.savetxt(subfolderfolder+prefix+'NewtonRes'+extension, resNewton_list)
				np.savetxt(subfolderfolder+prefix+'L2error'+extension, L2error)
				np.savetxt(subfolderfolder+prefix+'L2relerror'+extension, L2relerror)
				np.savetxt(subfolderfolder+prefix+'threshold'+extension, np.array(threshold_list))

	fig1, ax1 = plt.subplots(figsize=(8, 6))
	fig2, ax2 = plt.subplots(figsize=(8, 6))
	fig3, ax3 = plt.subplots(figsize=(8, 6))
	fig4, ax4 = plt.subplots(figsize=(8, 6))
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
			
			for caseplot, fig, ax in zip(range(1, 5), figs, axs):
				if caseplot == 1:
					yy = L2relerror; xx = nbInnerLoops[:len(L2relerror)]
					xlim = 10*np.ceil(np.max(nbInnerLoops)/10); ylim = np.power(10, np.floor(np.log10(np.min(L2relerror))))
					ylabel = r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$'
					xlabel = 'Number of inner iterations'
				elif caseplot == 2:
					yy = newtonRes; xx = nbInnerLoops[:len(newtonRes)]
					xlim = 10*np.ceil(np.max(nbInnerLoops)/10); ylim = np.power(10, np.floor(np.log10(np.min(newtonRes))))
					ylabel = 'Relative norm of outer residue'
					xlabel = 'Number of inner iterations'
				elif caseplot == 3:
					yy = newtonRes; xx = np.arange(0, len(newtonRes))
					xlim = 10; ylim = np.power(10, np.floor(np.log10(np.min(newtonRes))))
					ylabel = 'Relative norm of outer residue'
					xlabel = 'Number of outer iterations'
				elif caseplot == 4:
					yy = L2relerror; xx = np.arange(0, len(L2relerror))
					xlim = 10; ylim = np.power(10, np.floor(np.log10(np.min(L2relerror))))
					ylabel = r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$'
					xlabel = 'Number of outer iterations'

				ax.semilogy(xx, yy, label=legendname, marker=marker_list[l], linestyle=linestyle_list[l])
				ax.set_xlim(right=xlim, left=0)
				# ax.set_ylim(top=ylim[0], bottom=ylim[1])
				ax.set_ylim(top=1e1, bottom=ylim)
				ax.set_xlabel(xlabel)
				ax.set_ylabel(ylabel)
				ax.legend()
				fig.tight_layout()
				fig.savefig(folder+'NLConvergence_iters'+'_'+str(degree)+str(cuts)+str(caseplot)+'.pdf')

	fig, ax = plt.subplots(figsize=(8, 6))
	for [i, isadaptive], prefix1 in zip(enumerate([True]), ['inexact']):
		for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
			l = j + i*2
			legendname = prefix1.capitalize() + ' ' + prefix2.capitalize()
			prefix = prefix1 + '_' + prefix2 + '_'
			threshold = np.loadtxt(subfolderfolder+prefix+'threshold'+extension)

			yy = threshold; xx = np.arange(0, len(threshold))
			xlim = 10
			ylabel = 'Forcing term (threshold)'
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

