"""
.. Test of elasticity 2D
.. Infinite plate with a hole under uniaxial traction. 
.. The analytical solution of this problem is given by Timoshenko
.. The convergence curves are traced for IGA-Legendre and IGA-WQ
.. Joaquin Cornejo 
"""

from thesis.Elliptic.__init__ import *
import pickle
from pysrc.lib.lib_base import vtk2png

RUNSIMU = False
degList = np.arange(1, 5)
cutList = np.arange(1, 8)

if RUNSIMU:

	degree, cuts = 8, 8
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}
	problem, displacement = simulate_el(degree, cuts, quadArgs)[:2]
	np.save(FOLDER2DATA + 'dispel', displacement)
	with open(FOLDER2DATA + 'refpartel.pkl', 'wb') as outp:
		pickle.dump(problem.part, outp, pickle.HIGHEST_PROTOCOL)
	
	disp_ref = np.load(FOLDER2DATA + 'dispel.npy')
	with open(FOLDER2DATA + 'refpartel.pkl', 'rb') as inp:
		part_ref = pickle.load(inp)

	AbserrorList = np.ones((len(degList), len(cutList)))
	RelerrorList = np.copy(AbserrorList)

	for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
		quadArgs = {'quadrule': quadrule, 'type': quadtype}
		
		for i, degree in enumerate(degList):
			for j, cuts in enumerate(cutList):
				problem, displacement, _ = simulate_el(degree, cuts, quadArgs)
				AbserrorList[i, j], RelerrorList[i, j] = problem.normOfError(displacement, 
														normArgs={'type':'H1', 
																'part_ref':part_ref, 
																'u_ref':disp_ref})
			np.savetxt(FOLDER2DATA+'AbsError_el_'+quadrule+'_'+str(quadtype)+'.dat', AbserrorList)
			np.savetxt(FOLDER2DATA+'RelError_el_'+quadrule+'_'+str(quadtype)+'.dat', RelerrorList)	

problem, displacement, _ = simulate_el(4, 6, {'quadrule': 'wq', 'type': 2})
straintmp = problem.interpolate_strain(displacement)
strain = np.zeros((6, problem.part.nbqp_total))
strain[:2,:]=straintmp[:2,:]; strain[3,:]=straintmp[-1,:]
vonmises = problem.mechamaterial.evalElasticStress(strain)
problem.part.postProcessingDual(fields={'vms':vonmises}, folder=FOLDER2SAVE, name='ellipticel')
vtk2png(folder=FOLDER2SAVE, filename='ellipticel', fieldname='vms', cmap='coolwarm', title='Von Mises stress')

fig, ax = plt.subplots()
figname = FOLDER2SAVE + 'ConvergenceH1_el'
for quadrule, quadtype, plotpars in zip(['iga', 'wq'], ['leg', 2], [CONFIGLINE0, CONFIGLINE2]):
	quadArgs = {'quadrule': quadrule, 'type': quadtype}
	errorList = np.loadtxt(FOLDER2DATA+'RelError_el_'+quadrule+'_'+str(quadtype)+'.dat')

	for i, degree in enumerate(degList):
		color = COLORLIST[i]
		nbelList = 2**cutList

		if quadrule == 'iga': 
			ax.loglog(nbelList, errorList[i, :], label='Gauss quad. $p=$ '+str(degree), color=color, marker=plotpars['marker'], markerfacecolor='w',
					markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])
			slope = round(np.polyfit(np.log(nbelList), np.log(errorList[i, :]), 1)[0], 1)
			annotation.slope_marker((nbelList[-2],  errorList[i, -2]), slope, 
							poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)
					
		else: 
			ax.loglog(2**cutList, errorList[i, :], color=color, marker=plotpars['marker'], markerfacecolor='w',
				markersize=plotpars['markersize'], linestyle=plotpars['linestyle'])

# ax.loglog([], [], color='k', marker=CONFIGLINE1['marker'], markerfacecolor='w',
# 				markersize=CONFIGLINE1['markersize'], linestyle=CONFIGLINE1['linestyle'], label="IGA-WQ 1")
ax.loglog([], [], color='k', marker=CONFIGLINE2['marker'], markerfacecolor='w',
		markersize=CONFIGLINE2['markersize'], linestyle=CONFIGLINE2['linestyle'], label="Weighted quad.")

ax.set_ylabel('Relative '+r'$H^1$'+' error')
ax.set_xlabel('Number of elements by dimension')
ax.set_ylim(top=1e0, bottom=1e-10)
ax.set_xlim(left=1, right=200)
ax.legend(loc='lower left')
fig.tight_layout()
fig.savefig(figname+EXTENSION)