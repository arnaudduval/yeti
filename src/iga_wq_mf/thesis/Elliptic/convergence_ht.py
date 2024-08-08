"""
.. Test of steady heat transfer 2D
.. The geometry is a quart of a annulus
.. All the boundary conditions are considered Dirichlet-like 
.. The heat source is computed by f = -grad(k grad T) where T is a given function
.. We compute the relative L2 norm using the exact solution T
.. The convergence curves are traced for IGA-Legendre and IGA-WQ
.. Joaquin Cornejo 
"""

from thesis.Elliptic.__init__ import *
from pysrc.lib.lib_base import vtk2png

# Set global variables
RUNSIMU = False
degList = np.arange(1, 5)
cutList = np.arange(1, 8)

if RUNSIMU:

	AbserrorList = np.ones((len(degList), len(cutList)))
	RelerrorList = np.copy(AbserrorList)

	for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
		quadArgs = {'quadrule': quadrule, 'type': quadtype}
		
		for i, degree in enumerate(degList):
			for j, cuts in enumerate(cutList):
				problem, temperature, _ = simulate_ht(degree, cuts, quadArgs)
				AbserrorList[i, j], RelerrorList[i, j] = problem.normOfError(temperature, 
														normArgs={'type':'L2',
														'exactFunction':exactTemperature_quartCircle,
														'exactFunctionDers':exactDiffTemperature_quartCircle})
			np.savetxt(FOLDER2DATA+'AbsError_ht_'+quadrule+'_'+str(quadtype)+'.dat', AbserrorList)
			np.savetxt(FOLDER2DATA+'RelError_ht_'+quadrule+'_'+str(quadtype)+'.dat', RelerrorList)

problem, temperature, _ = simulate_ht(4, 6, {'quadrule': 'wq', 'type': 2})
problem.part.postProcessingPrimal(fields={'temp':temperature}, folder=FOLDER2SAVE, name='ellipticht')
vtk2png(folder=FOLDER2SAVE, filename='ellipticht', fieldname='temp', cmap='coolwarm', title='Temperature')

fig, ax = plt.subplots(figsize=(6, 5))
figname = FOLDER2SAVE + 'ConvergenceL2_ht'
for quadrule, quadtype, plotpars in zip(['iga', 'wq', 'wq'], ['leg', 1, 2], [CONFIGLINE0, CONFIGLINE1, CONFIGLINE2]):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	errorList = np.loadtxt(FOLDER2DATA+'RelError_ht_'+quadrule+'_'+str(quadtype)+'.dat')

	for i, degree in enumerate(degList):
		color = COLORLIST[i]
		nbelList = 2**cutList

		if quadrule == 'iga': 
			ax.loglog(nbelList, errorList[i, :], label='IGA-GL deg. '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			slope = round(np.polyfit(np.log(nbelList), np.log(errorList[i, :]), 1)[0], 1)
			annotation.slope_marker((nbelList[-2],  errorList[i, -2]), slope, 
							poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
					
		else: 
			ax.loglog(2**cutList, errorList[i, :], color=color, marker=plotpars['marker'], markerfacecolor='w',
				markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])

ax.loglog([], [], color='k', marker=CONFIGLINE1['marker'], markerfacecolor='w',
				markersize=CONFIGLINE1['markersize'], linestyle=CONFIGLINE1['linestyle'], label="IGA-WQ 1")
ax.loglog([], [], color='k', marker=CONFIGLINE2['marker'], markerfacecolor='w',
		markersize=CONFIGLINE2['markersize'], linestyle=CONFIGLINE2['linestyle'], label="IGA-WQ 2")

ax.set_ylabel('Relative '+r'$L^2$'+' error')
ax.set_xlabel('Number of elements by dimension')
ax.set_ylim(top=1e0, bottom=1e-10)
ax.set_xlim(left=1, right=200)
ax.legend(loc='lower left')
fig.tight_layout()
fig.savefig(figname+'.pdf')