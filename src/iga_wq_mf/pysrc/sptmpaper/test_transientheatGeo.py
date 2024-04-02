from pysrc.sptmpaper.input_data import *
from pyevtk.vtk import VtkGroup

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/TRnonlin4/'
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
FIG_CASE = 1

if FIG_CASE == 1:

	lastsufix = 'linear' if ISLINEAR else 'nonlin'
	filename1 = folder + 'incrementalheat'+lastsufix
	filename2 = folder + 'spacetimeheat'+lastsufix

	degree_list = np.array([1, 2, 3, 4, 5])
	cuts_list   = np.arange(1, 6)

	if TODOSIMU:
		exportTimeDependentMaterial(np.linspace(0, 1, 33), 
						temperature=exactTemperatureQuad_inc,
						fields={'mat':nonlinearfunc, 'temp':exactTemperatureQuad_inc}, 
						geoArgs={'name': 'QA', 'degree': 5*np.ones(3, dtype=int), 
								'nb_refinementByDirection': 5*np.ones(3, dtype=int)},
						folder=subfolder,)
		run(folder=subfolder)

		error_list1 = np.ones((len(degree_list), len(cuts_list), 2**CUTS_TIME))
		error_list2 = np.ones((len(degree_list), len(cuts_list)))
		
		for j, cuts in enumerate(cuts_list):
			for i, degree in enumerate(degree_list):
				geoArgs = {'name': 'TP', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, 1, 1])}

				# Incremental problem
				problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, 
												powerdensity=powerDensitySquare_inc, geoArgs=geoArgs)
				for k, t in enumerate(time_inc[1:-1]):
					error_list1[i, j, k], _ = problem_inc.normOfError(temp_inc[:, k+1], 
																normArgs={'type':'L2',
																		'exactFunction':exactTemperatureQuad_inc,
																		'exactExtraArgs':{'time':t}})
				# Space time problem
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, 
													powerdensity=powerDensitySquare_spt,
													degree_spt=degree, geoArgs=geoArgs)

				error_list2[i, j], _ = problem_spt.normOfError(temp_spt, 
													normArgs={'type':'L2',
															'exactFunction':exactTemperatureQuad_spt})
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
			ax.set_ylim(top=1e5, bottom=1e-2)
			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
			fig.tight_layout()
			fig.savefig(subfolder + 'FigConvergIncrHeat' + str(j+1) +  '.pdf')
			plt.close(fig)

		error_list = np.load(filename2+'.npy')
		fig, ax = plt.subplots(figsize=(9, 6))
		for i, degree in enumerate(degree_list):
			color = COLORLIST[i]
			ax.loglog(2**cuts_list, error_list[i, :], color=color, marker='o', markerfacecolor='w',
						markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))
			
		ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L_2(\Pi)}$')
		ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
		ax.set_ylim(top=1e5, bottom=1e-2)
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig.tight_layout()
		fig.savefig(subfolder + 'FigConvergSptHeat' +  '.pdf')
		plt.close(fig)


elif FIG_CASE == 2:
	extension = '.dat'
	degree, cuts = 4, 5
	geoArgs = {'name': 'TP', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': np.array([cuts, 1, 1])}

	subfolderfolder = folder + '_' + str(degree) + '_' + str(cuts) + '/' 
	if not os.path.isdir(subfolderfolder): os.mkdir(subfolderfolder)

	if TODOSIMU:
		for [i, isadaptive], prefix1 in zip(enumerate([False, True]), ['exact', 'inexact']):
			for [j, isfull], prefix2 in zip(enumerate([True, False]), ['newton', 'picard']):
				prefix = prefix1 + '_' + prefix2 + '_'
				output = {}
				blockPrint()
				problem_spt = simulate_spacetime(degree, cuts, 
												powerdensity=powerDensitySquare_spt,
												degree_spt=degree, geoArgs=geoArgs,
												isadaptive=isadaptive, isfull=isfull,
												outputArgs=output)[0]
				
				displmnt_list  = output['Solution']; resKrylov_list = output['KrylovRes']
				resNewton_list = output['NewtonRes']; threshold_list = output['Threshold']
				L2error, L2relerror = [], []

				for displacement in displmnt_list:
					err, relerr  = problem_spt.normOfError(displacement, normArgs={'type':'L2', 
																	'exactFunction':exactTemperatureSquare_spt})
					L2error.append(err); L2relerror.append(relerr)
				resKrylovclean = np.array([]); counter_list = [0]
				for _ in resKrylov_list: 
					resKrylovclean = np.append(resKrylovclean, _[np.nonzero(_)])
					counter_list.append(counter_list[-1] + len(_[np.nonzero(_)]))
				enablePrint()
				np.savetxt(subfolderfolder+prefix+'CumulKrylovRes'+extension, resKrylovclean)
				np.savetxt(subfolderfolder+prefix+'Inner_loops'+extension, counter_list)
				np.savetxt(subfolderfolder+prefix+'NewtonRes'+extension, resNewton_list)
				np.savetxt(subfolderfolder+prefix+'L2error'+extension, L2error)
				np.savetxt(subfolderfolder+prefix+'L2relerror'+extension, L2relerror)
				np.savetxt(subfolderfolder+prefix+'threshold'+extension, np.array(threshold_list))
	else:
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

