from pysrc.lib.__init__ import *
from pysrc.sptmpaper.input_data import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/spttransient/'
subfolder = folder +  'steps/'
if not os.path.isdir(folder): os.mkdir(folder)
if not os.path.isdir(subfolder): os.mkdir(subfolder)

# Set global variables
TODOSIMU = True

lastsufix = 'linear' if ISLINEAR else 'nonlin'
filename = folder + 'spacetimeheat' + lastsufix

degree_list = np.array([1, 2, 3, 4, 5])
cuts_list   = np.arange(1, 7)

if TODOSIMU:
	error_list = np.ones((len(degree_list), len(cuts_list)))
	for j, cuts in enumerate(cuts_list):
		for i, degree in enumerate(degree_list):
			problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, 
												powerdensity=powerDensity_spt,
												degree_spt=degree)

			error_list[i, j], _ = problem_spt.normOfError(temp_spt, 
												normArgs={'type':'L2',
														'exactFunction':exactTemperature_spt})
	np.save(filename, error_list)

else:
	error_list = np.load(filename+'.npy')
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
