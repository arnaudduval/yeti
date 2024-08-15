from thesis.Incremental.__init__ import *
from thesis.SpaceTime.input_data import *

# Set global variables
SUFIX = ('lin' if ISLINEAR else 'nonlin') + '1d'
RUNSIMU = False
EXTENSION = '.dat'

if RUNSIMU: assert IS1DIM, 'Try 1D methods'

degree, cuts = 8, 7
quadArgs = {'quadrule':'iga', 'type':'leg'}
nbelincList = np.array([2**cuts for cuts in range(2, 8)])
degsptList  = np.arange(1, 4)
abserrorInc, relerrorInc = np.ones(len(nbelincList)), np.ones(len(nbelincList))
abserrorInc2, relerrorInc2 = np.ones(len(nbelincList)), np.ones(len(nbelincList))
abserrorSpt, relerrorSpt = np.ones((len(degsptList), len(nbelincList))), np.ones((len(degsptList), len(nbelincList)))

if RUNSIMU:
	for i, nbelsinc in enumerate(nbelincList):

		# Incremental
		dirichlet_table = np.ones((2, 2))
		problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensitySquare_inc, 
													dirichlet_table=dirichlet_table,
													nbel_time=nbelsinc, 
													quadArgs={'quadrule':'iga'}, 
													is1dim=IS1DIM, alpha=0.5)
		
		abserrorInc[i], relerrorInc[i] = problem_inc.normOfError(temp_inc[:, -1], 
													normArgs={'type':'L2',
															'exactFunction':exactTemperatureSquare_inc,
															'exactExtraArgs':{'time':time_inc[-1]}})
		
		np.savetxt(FOLDER2DATA+'2abserrorstag_inc'+SUFIX+EXTENSION, abserrorInc)
		np.savetxt(FOLDER2DATA+'2relerrorstag_inc'+SUFIX+EXTENSION, relerrorInc)

		# # Space time
		# for j, degspt in enumerate(degsptList):
		# 	dirichlet_table = np.ones((2, 2)); dirichlet_table[-1, 1] = 0
		# 	problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensitySquare_spt, 
		# 											dirichlet_table=dirichlet_table,
		# 											degree_time=degspt, isfull=True,
		# 											nbel_time=nbelsinc+1-degspt, 
		# 											quadArgs=quadArgs, is1dim=IS1DIM)
				
		# 	newtemp_spt = np.reshape(temp_spt, newshape=(problem_spt.part.nbctrlpts_total, problem_spt.time.nbctrlpts_total), order='F')
		# 	abserrorSpt[j, i], relerrorSpt[j, i] = problem_inc.normOfError(newtemp_spt[:, -1], 
		# 											normArgs={'type':'L2',
		# 														'exactFunction':exactTemperatureSquare_inc,
		# 														'exactExtraArgs':{'time':time_inc[-1]}})

		# 	np.savetxt(FOLDER2DATA+'2abserrorstag_spt'+SUFIX+EXTENSION, abserrorSpt)
		# 	np.savetxt(FOLDER2DATA+'2relerrorstag_spt'+SUFIX+EXTENSION, relerrorSpt)
		
		# # enablePrint()

fig, ax = plt.subplots()
errorList1 = np.loadtxt(FOLDER2DATA+'2relerrorstag_spt'+SUFIX+EXTENSION)
for i, deg in enumerate(degsptList):
	nbctrlpts = nbelincList+deg
	ax.loglog(nbctrlpts, errorList1[i, :], color=COLORLIST[i], marker='s', markerfacecolor='w',
				label='SPT-IGA deg. '+str(int(deg)))
	slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[i, 3:]), 1)[0]
	slope = round(slope, 1)
	if deg == 3:
		annotation.slope_marker((nbctrlpts[2], errorList1[i, 2]), slope, 
					poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
	else:
		annotation.slope_marker((nbctrlpts[-2], errorList1[i, -2]), slope, 
						poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

nbctrlpts = nbelincList+1
errorList1 = np.loadtxt(FOLDER2DATA+'2relerrorstag_inc'+SUFIX+EXTENSION)
ax.loglog(nbctrlpts, errorList1, marker=CONFIGLINE4['marker'], markerfacecolor='w', color='k',
				markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'], label='INC-IGA '+r'$\theta=0.5$')
slope = np.polyfit(np.log10(nbctrlpts[3:]),np.log10(errorList1[3:]), 1)[0]
slope = round(slope, 1)
annotation.slope_marker((nbctrlpts[-2], errorList1[-2]), slope, 
				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)

ax.set_ylabel('Relative ' + r'$L^2$'+' error at last time-step')
ax.set_xlabel('Number of control points on time \n(or time-steps)')
ax.set_xlim(left=2, right=300)
ax.set_ylim(top=1e-1, bottom=1e-8)
ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.2), loc='upper center')
fig.tight_layout()
fig.savefig(FOLDER2SAVE+'StagnationError'+SUFIX+'.pdf')