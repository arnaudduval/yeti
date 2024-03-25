from pysrc.sptmpaper.input_data import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/TR/'
subfolder = folder +  'steps/'
if not os.path.isdir(folder): os.mkdir(folder)
if not os.path.isdir(subfolder): os.mkdir(subfolder)

# Set global variables
TODOSIMU = True
FIG_CASE = 1

if FIG_CASE == 1:

	lastsufix = 'linear' if ISLINEAR else 'nonlin'
	filename1 = folder + 'incrementalheat'+lastsufix
	filename2 = folder + 'spacetimeheat'+lastsufix

	degree_list = np.array([1, 2, 3, 4, 5])
	cuts_list   = np.arange(1, 7)

	if TODOSIMU:
		error_list1 = np.ones((len(degree_list), len(cuts_list), 2**CUTS_TIME))
		error_list2 = np.ones((len(degree_list), len(cuts_list)))
		
		for j, cuts in enumerate(cuts_list):
			for i, degree in enumerate(degree_list):
				geoArgs = {'name': 'TP', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, 1, 1])}

				# Incremental problem
				problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, 
												powerdensity=powerDensity_inc, geoArgs=geoArgs)
				for k, t in enumerate(time_inc[1:-1]):
					error_list1[i, j, k], _ = problem_inc.normOfError(temp_inc[:, k+1], 
																normArgs={'type':'L2',
																		'exactFunction':exactTemperature_inc,
																		'exactExtraArgs':{'time':t}})
				# Space time problem
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, 
													powerdensity=powerDensity_spt,
													degree_spt=degree, geoArgs=geoArgs)

				error_list2[i, j], _ = problem_spt.normOfError(temp_spt, 
													normArgs={'type':'L2',
															'exactFunction':exactTemperature_spt})
		np.save(filename1, error_list1)
		np.save(filename2, error_list2)

	else:
		error_list = np.load(filename1+'.npy')
		for j, k in enumerate(range(0, np.size(error_list, axis=2), 4)):
			fig, ax = plt.subplots(figsize=(9, 6))
		
			for i, degree in enumerate(degree_list):
				color = COLORLIST[i]
				ax.loglog(2**cuts_list, error_list[i, :, k], color=color, marker='o', markerfacecolor='w',
							markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))

			ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L_2(\Omega)}$')
			ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
			ax.set_ylim(top=1e2, bottom=1e-5)
			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
			fig.tight_layout()
			fig.savefig(subfolder + 'FigConvergIncrHeat' + str(j+1) +  '.pdf')
			plt.close(fig)

		error_list = np.load(filename1+'.npy')
		fig, ax = plt.subplots(figsize=(9, 6))
		for i, degree in enumerate(degree_list):
			color = COLORLIST[i]
			ax.loglog(2**cuts_list, error_list[i, :], color=color, marker='o', markerfacecolor='w',
						markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))
			
		ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L_2(\Pi)}$')
		ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
		ax.set_ylim(top=1e2, bottom=1e-5)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(subfolder + 'FigConvergSptHeat' +  '.pdf')
		plt.close(fig)


elif FIG_CASE == 2:
	...
	# if FIG_CASE == 1:
	# 	normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
	# 	onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
	# 	onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}
	# 	plotoptions = [normalPlot, onlyMarker1, onlyMarker2]

	# 	lastsufix = 'linear' if ISLINEAR else 'nonlin'
	# 	figname = folder + 'SPTNonLinearConvergenceL2'+lastsufix+'.pdf'
	# 	filenames = ['L2relerror_meshpar_iga_leg_']

	# 	normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
	# 	fig, ax = plt.subplots(figsize=(8, 6))

	# 	for filename, plotops in zip(filenames, plotoptions):
	# 		quadrule = filename.split('_')[2]
	# 		table = np.loadtxt(folder+filename+lastsufix+extension)	
	# 		nbels   = 2**(table[0, 1:])
	# 		degrees = table[1:, 0]
	# 		errors  = table[1:, 1:]
	# 		for i, degree in enumerate(degrees):
	# 			color = COLORLIST[i]
	# 			if quadrule == 'iga': 
	# 				ax.loglog(nbels, errors[i, :], label='IGA-GL deg. '+str(int(degree)), color=color, marker=plotops['marker'],
	# 							markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
	# 				# slope = np.polyfit(np.log10(nbels[2:]),np.log10(errors[i, 2:]), 1)[0]
	# 				# slope = round(slope, 1)
	# 				# annotation.slope_marker((nbels[-1], errors[i, -1]), slope, 
	# 				# 				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)			
	# 			else: 
	# 				ax.loglog(nbels, errors[i, :], color=color, marker=plotops['marker'], markerfacecolor='w',
	# 						markersize=plotops['markersize'], linestyle=plotops['linestyle'])
						
	# 			fig.savefig(figname)

	# 	# ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
	# 	# 				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ 2")
	# 	# ax.loglog([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
	# 	# 		markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label="IGA-WQ 4")

	# 	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$')
	# 	ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
	# 	ax.set_xlim(left=1, right=100)
	# 	ax.set_ylim(top=1e1, bottom=1e-7)
	# 	ax.legend(loc='lower left')
	# 	fig.tight_layout()
	# 	fig.savefig(figname)