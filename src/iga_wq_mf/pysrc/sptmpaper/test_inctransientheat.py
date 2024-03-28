"""
.. Test of 1D transient heat transfer 
.. Author: Joaquin Cornejo
"""

from pysrc.lib.__init__ import *
from pysrc.sptmpaper.input_data import *

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/inctransient/'
subfolder = folder +  'steps/'
if not os.path.isdir(folder): os.mkdir(folder)
if not os.path.isdir(subfolder): os.mkdir(subfolder)

# Set global variables
TODOSIMU = True

if IS1DIM:
	# Trace material properties
	XX, TIME = np.meshgrid(np.linspace(0, 1, 201), np.linspace(0, 1, 129))
	TEMPERATURE = exactTemperature_inc(args={'position':XX, 'time':TIME})

	CONDUCTIVITY = conductivityProperty(args={'temperature':TEMPERATURE})
	for matfield, figname in zip([CONDUCTIVITY, TEMPERATURE], 
								['conductivity', 'temperature']):
		fig, ax = plt.subplots(figsize=(10, 4))
		im = ax.contourf(XX, TIME, matfield, 21, cmap='viridis')
		cbar = plt.colorbar(im)
		cbar.set_label(figname.capitalize())

		ax.grid(False)
		ax.set_ylabel('Time')
		ax.set_xlabel('Position')
		fig.tight_layout()
		fig.savefig(folder + 'transHeat_' + figname)

	problem_inc, time_inc, temperature_inc = simulate_incremental(2, 3, 
											powerdensity=powerDensity_inc, is1dim=IS1DIM)
	TEMPERATURE_INTERP = problem_inc.interpolateMeshgridField(temperature_inc, sampleSize=201)[0]
	CONDUCTIVITY_INTERP = conductivityProperty(args={'temperature':TEMPERATURE_INTERP.T})
	TEMPDIFF = np.abs(TEMPERATURE - TEMPERATURE_INTERP.T)
	CONDDIFF = np.abs(CONDUCTIVITY - CONDUCTIVITY_INTERP)

	for matfield, figname in zip([CONDUCTIVITY_INTERP, TEMPERATURE_INTERP.T, CONDDIFF, TEMPDIFF], 
								['conductivity', 'temperature', 'error conductivity', 'error temperature']):
		fig, ax = plt.subplots(figsize=(10, 4))
		im = ax.contourf(XX, TIME, matfield, 21, cmap='viridis')
		cbar = plt.colorbar(im)
		cbar.set_label(figname.capitalize())

		ax.grid(False)
		ax.set_ylabel('Time')
		ax.set_xlabel('Position')
		fig.tight_layout()
		fig.savefig(folder + 'transHeatinterp_' + figname)
	
	# for k, t in enumerate(time_inc[1:-1]):
	# 	error = problem_inc.normOfError(temperature_inc[:, k+1], 
	# 								normArgs={'type':'L2',
	# 										'exactFunction':exactTemperature_inc,
	# 										'exactExtraArgs':{'time':t}})
	# 	print('Step:%d, RelError:%.3e' %(k, error[-1]))

lastsufix = 'linear' if ISLINEAR else 'nonlin'
filename = folder + 'incrementalheat' + lastsufix

degree_list = np.array([1, 2, 3, 4, 5])
cuts_list   = np.arange(1, 6)

if TODOSIMU:
	error_list = np.ones((len(degree_list), len(cuts_list), 2**CUTS_TIME))
	for j, cuts in enumerate(cuts_list):
		for i, degree in enumerate(degree_list):
			problem_inc, time_inc, TEMPERATURE_INTERP = simulate_incremental(degree, cuts, 
											powerdensity=powerDensity_inc, is1dim=IS1DIM)
			for k, t in enumerate(time_inc[1:-1]):
				error_list[i, j, k], _ = problem_inc.normOfError(TEMPERATURE_INTERP[:, k+1], 
															normArgs={'type':'L2',
																	'exactFunction':exactTemperature_inc,
																	'exactExtraArgs':{'time':t}})
	np.save(filename, error_list)

error_list = np.load(filename+'.npy')
for j, k in enumerate(range(0, np.size(error_list, axis=2), 4)):
	fig, ax = plt.subplots(figsize=(9, 6))

	for i, degree in enumerate(degree_list):
		color = COLORLIST[i]
		ax.loglog(5*2**(cuts_list-1), error_list[i, :, k], color=color, marker='o', markerfacecolor='w',
					markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))

	ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L_2(\Omega)}$')
	ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
	ax.set_ylim(top=1e2, bottom=1e-5)
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	fig.tight_layout()
	fig.savefig(subfolder + 'FigConvergIncrHeat' + str(j+1) +  '.pdf')
	plt.close(fig)
