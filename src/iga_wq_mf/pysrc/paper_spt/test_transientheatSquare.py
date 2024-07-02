"""
.. Test of 1D transient heat transfer 
.. Author: Joaquin Cornejo
"""

from pysrc.lib.__init__ import *
from pysrc.sptmpaper.input_data import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/Pseudo2D/lin1d/'
subfolder = folder +  'steps/'
if not os.path.isdir(folder): os.mkdir(folder)
if not os.path.isdir(subfolder): os.mkdir(subfolder)

# Set global variables
TODOSIMU = False
lastsufix = 'linear' if ISLINEAR else 'nonlin'
degree_list = np.array([1, 2, 3, 4, 5])
cuts_list   = np.arange(1, 7)

if IS1DIM:
	# Trace material properties
	XX, TIME = np.meshgrid(np.linspace(0, 1, 201), np.linspace(0, 1, 129))
	TEMPERATURE = exactTemperatureSquare_inc(args={'position':XX, 'time':TIME})

	CONDUCTIVITY = nonlinearfunc(args={'temperature':TEMPERATURE})
	for matfield, figname in zip([CONDUCTIVITY, TEMPERATURE], 
								['conductivity', 'temperature']):
		fig, ax = plt.subplots(figsize=(10, 4))
		im = ax.contourf(XX, TIME, matfield, 21, cmap='viridis')
		cbar = plt.colorbar(im, format='%.1e')
		cbar.set_label(figname.capitalize())

		ax.grid(False)
		ax.set_ylabel('Time')
		ax.set_xlabel('Position')
		fig.tight_layout()
		fig.savefig(folder + 'transHeat_' + figname)

	problem_inc, time_inc, temperature_inc = simulate_incremental(4, 6, 
											powerDensitySquare_inc, is1dim=IS1DIM)
	TEMPERATURE_INTERP = problem_inc.interpolateMeshgridField(temperature_inc, sampleSize=201)[0]
	CONDUCTIVITY_INTERP = nonlinearfunc(args={'temperature':TEMPERATURE_INTERP.T})
	TEMPDIFF = np.abs(TEMPERATURE - TEMPERATURE_INTERP.T)
	TEMPDIFF = np.where(np.abs(TEMPERATURE)<1e-12, 0.0, TEMPDIFF/np.abs(TEMPERATURE))
	CONDDIFF = np.abs(CONDUCTIVITY - CONDUCTIVITY_INTERP)/np.abs(CONDUCTIVITY)

	for matfield, figname in zip([CONDUCTIVITY_INTERP, TEMPERATURE_INTERP.T, CONDDIFF, TEMPDIFF], 
								['conductivity', 'temperature', 'error conductivity', 'error temperature']):
		fig, ax = plt.subplots(figsize=(10, 4))
		im = ax.contourf(XX, TIME, matfield, 21, cmap='viridis')
		cbar = plt.colorbar(im, format='%.1e')
		cbar.set_label(figname.capitalize())

		ax.grid(False)
		ax.set_ylabel('Time')
		ax.set_xlabel('Position')
		fig.tight_layout()
		fig.savefig(folder + 'transHeatinterp_' + figname)

# ------------------------------------------------------------------
filename = folder + 'incrementalheat' + lastsufix
if TODOSIMU:
	error_list = np.ones((len(degree_list), len(cuts_list), 2**CUTS_TIME))
	for j, cuts in enumerate(cuts_list):
		for i, degree in enumerate(degree_list):
			problem_inc, time_inc, TEMPERATURE_INTERP = simulate_incremental(degree, cuts, 
											powerDensitySquare_inc, is1dim=IS1DIM)
			for k, t in enumerate(time_inc[1:-1]):
				_, error_list[i, j, k] = problem_inc.normOfError(TEMPERATURE_INTERP[:, k+1], 
															normArgs={'type':'L2',
																	'exactFunction':exactTemperatureSquare_inc,
																	'exactExtraArgs':{'time':t}})
	
	np.save(filename, error_list)

error_list = np.load(filename+'.npy')
for j, k in enumerate(range(0, np.size(error_list, axis=2), 4)):
	fig, ax = plt.subplots(figsize=(8, 6))

	for i, degree in enumerate(degree_list):
		color = COLORLIST[i]
		ax.loglog(2**cuts_list, error_list[i, :, k], color=color, marker='o', markerfacecolor='w',
					markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))

	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}$')
	if IS1DIM: ax.set_xlabel('No. of elements')
	else: ax.set_xlabel('No. of elements in main direction')
	ax.set_ylim(top=1e1, bottom=1e-9)
	ax.set_xlim(left=1e0, right=80)
	# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(subfolder + 'FigConvergIncrHeat' + str(j+1) +  '.pdf')
	plt.close(fig)

# ------------------------------------------------------------------
filename = folder + 'spacetimeheat' + lastsufix
if not IS1DIM:
	if TODOSIMU:
		error_list = np.ones((len(degree_list), len(cuts_list)))
		for j, cuts in enumerate(cuts_list):
			for i, degree in enumerate(degree_list):
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, 
													powerDensitySquare_spt)

				_, error_list[i, j] = problem_spt.normOfError(temp_spt, 
													normArgs={'type':'L2',
															'exactFunction':exactTemperatureSquare_spt})
		np.save(filename, error_list)

	error_list = np.load(filename+'.npy')
	fig, ax = plt.subplots(figsize=(8, 6))

	for i, degree in enumerate(degree_list):
		color = COLORLIST[i]
		ax.loglog(2**cuts_list, error_list[i, :], color=color, marker='o', markerfacecolor='w',
					markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))
		
	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$')
	if IS1DIM: ax.set_xlabel('No. of spatial elements')
	else: ax.set_xlabel('No. of elements in main spatial direction')
	ax.set_ylim(top=1e1, bottom=1e-9)
	ax.set_xlim(left=1e0, right=80)
	# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(folder + 'FigConvergSptHeat' +  '.pdf')
	plt.close(fig)