from thesis.Incremental.__init__ import *
from thesis.SpaceTime.input_data import *

def exactTemperature_inc(args):
	func = None
	if GEONAME == 'QA': func = exactTemperatureRing_inc(args)
	elif GEONAME == 'TP': func = exactTemperatureTrap_inc(args)
	else: raise Warning('Not possible')
	return func

def powerDensity_inc(args):
	func = None
	if GEONAME == 'QA': func = powerDensityRing_inc(args)
	elif GEONAME == 'TP': func = powerDensityTrap_inc(args)
	else: raise Warning('Not possible')
	return func

def exactTemperature_spt(args):
	func = None
	if GEONAME == 'QA': func = exactTemperatureRing_spt(args)
	elif GEONAME == 'TP': func = exactTemperatureTrap_spt(args)
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
EXTENSION = '.dat'

if RUNSIMU: assert (not IS1DIM), 'Try 2D methods'

degree, cuts = 8, 6
quadArgs = {'quadrule':'wq', 'type':2}
nbelincList = np.array([2**cuts for cuts in range(1, 7)])

if RUNSIMU:
	for k, alpha in enumerate([0.5, 1.0]):
		abserrorInc, relerrorInc = np.ones(len(nbelincList)), np.ones(len(nbelincList))

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
														geoArgs=geoArgs, nbel_time=nbelinc, quadArgs=quadArgs, alpha=alpha)
			
			abserrorInc[i], relerrorInc[i] = problem_spt_inc.normOfError(np.ravel(temp_inc, order='F'), 
														normArgs={'type':'L2',
																'exactFunction':exactTemperature_spt})
							
			np.savetxt(FOLDER2DATA+'2abserrorstag_inc'+str(k)+SUFIX+EXTENSION, abserrorInc)
			np.savetxt(FOLDER2DATA+'2relerrorstag_inc'+str(k)+SUFIX+EXTENSION, relerrorInc)

fig, ax = plt.subplots()

if PLOTRELATIVE: errorList1 = np.loadtxt(FOLDER2DATA+'2relerrorstag_inc0'+SUFIX+EXTENSION)
else: errorList1 = np.loadtxt(FOLDER2DATA+'2abserrorstag_inc0'+SUFIX+EXTENSION)

nbctrlpts = nbelincList+1
ax.loglog(nbctrlpts, errorList1, marker=CONFIGLINE4['marker'], markerfacecolor='w',
				markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'], 
				label='IGA-GL '+r'$\alpha=0.5$')
slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[3:]), 1)[0]
slope = round(slope, 1)
annotation.slope_marker((nbctrlpts[-2], errorList1[-2]), slope, 
				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

if PLOTRELATIVE: errorList1 = np.loadtxt(FOLDER2DATA+'2relerrorstag_inc1'+SUFIX+EXTENSION)
else: errorList1 = np.loadtxt(FOLDER2DATA+'2abserrorstag_inc1'+SUFIX+EXTENSION)

nbctrlpts = nbelincList+1
ax.loglog(nbctrlpts, errorList1, marker=CONFIGLINE5['marker'], markerfacecolor='w',
				markersize=CONFIGLINE5['markersize'], linestyle=CONFIGLINE5['linestyle'], 
				label='IGA-GL '+r'$\alpha=1.0$')
slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[3:]), 1)[0]
slope = round(slope, 1)
annotation.slope_marker((nbctrlpts[-2], errorList1[-2]), slope, 
				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

if PLOTRELATIVE: 
	ax.set_ylabel('Relative stagnation error in time')
	ax.set_ylim(top=1e0, bottom=1e-4)
else: 
	ax.set_ylabel('Stagnation error in time')
	ax.set_ylim(top=1e1, bottom=1e-2)

ax.set_xlabel('Number of increments in time')
ax.set_xlim(left=2, right=100)
ax.legend(loc='lower left')
fig.tight_layout()
fig.savefig(FOLDER2SAVE+'StagnationError'+SUFIX+'.pdf')
