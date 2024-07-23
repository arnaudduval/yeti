from thesis.SpaceTime.__init__ import *
from thesis.SpaceTime.input_data import *

# Set global variables
SUFIX = ('lin' if ISLINEAR else 'nonlin') + '1d'
RUNSIMU = False
FIG_CASE = 2
EXTENSION = '.dat'

if RUNSIMU: assert IS1DIM, 'Try 1D methods'

if FIG_CASE == 0:
	if RUNSIMU:
		degList = np.array([1, 2, 3, 4, 5])
		cutList = np.arange(1, 6)	
		for quadrule, quadtype in zip(['iga', 'wq', 'wq'], ['leg', 1, 2]):
			sufix = '_' + quadrule + '_' + str(quadtype) + '_' + SUFIX
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			Abserror = np.zeros((len(degList)+1, len(cutList)+1))
			relerror = np.zeros((len(degList)+1, len(cutList)+1))
			Abserror[0, 1:] = cutList; relerror[0, 1:] = cutList
			Abserror[1:, 0] = degList; relerror[1:, 0] = degList
			filenameA1 = FOLDER2DATA+'L2abserror'+sufix+EXTENSION
			filenameR1 = FOLDER2DATA+'L2relerror'+sufix+EXTENSION
			for j, cuts in enumerate(cutList):
				for i, degree in enumerate(degList):
					blockPrint()
					dirichlet_table = np.ones((2, 2)); dirichlet_table[-1, 1] = 0
					problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensitySquare_spt, 
													dirichlet_table=dirichlet_table, quadArgs=quadArgs, 
													degree_time=degree, nbel_time=2**cuts, is1dim=IS1DIM)
					enablePrint()					
					Abserror[i+1, j+1], relerror[i+1, j+1] = problem_spt.normOfError(temp_spt, 
																	normArgs={'type':'L2', 
																			'exactFunction':exactTemperatureSquare_spt},)

					np.savetxt(filenameA1, Abserror)
					np.savetxt(filenameR1, relerror)

	plotoptions = [CONFIGLINE0, CONFIGLINE1, CONFIGLINE2]
	figname = FOLDER2SAVE + 'L2Convergence' + SUFIX + '.pdf'
	filenames = ['L2abserror_iga_leg_', 'L2abserror_wq_1_', 'L2abserror_wq_2_']

	fig, ax = plt.subplots(figsize=(8, 6))
	for filename, plotops in zip(filenames, plotoptions):
		quadrule = filename.split('_')[1]
		table = np.loadtxt(FOLDER2DATA + filename + SUFIX + EXTENSION)	
		nbels = 2**(table[0, 1:])
		degList = table[1:, 0]
		errList = table[1:, 1:]
		for i, degree in enumerate(degList):
			color = COLORLIST[i]
			if quadrule == 'iga': 
				ax.loglog(nbels, errList[i, :], label='IGA-GL deg. '+str(int(degree)), color=color, marker=plotops['marker'],
							markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])		
			else: 
				ax.loglog(nbels, errList[i, :], color=color, marker=plotops['marker'], markerfacecolor='w',
						markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
	ax.loglog([], [], color='k', marker=CONFIGLINE1['marker'], markerfacecolor='w',
					markersize=CONFIGLINE1['markersize'], linestyle=CONFIGLINE1['linestyle'], label='IGA-WQ 1')
	ax.loglog([], [], color='k', marker=CONFIGLINE2['marker'], markerfacecolor='w',
			markersize=CONFIGLINE2['markersize'], linestyle=CONFIGLINE2['linestyle'], label='IGA-WQ 2')

	ax.set_ylabel(r'$L^2$' + ' error')
	ax.set_xlabel('Number of elements by space-time direction')
	ax.set_xlim(left=1, right=50)
	ax.set_ylim(top=1e2, bottom=1e-8)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(figname)

elif FIG_CASE == 1:

	filenameA1 = FOLDER2DATA + '1incheatAbs'+SUFIX
	filenameR1 = FOLDER2DATA + '1incheatRel'+SUFIX
	filenameT1 = FOLDER2DATA + '1incheatTim'+SUFIX

	filenameA2 = FOLDER2DATA + '1sptheatAbs'+SUFIX
	filenameR2 = FOLDER2DATA + '1sptheatRel'+SUFIX
	filenameT2 = FOLDER2DATA + '1sptheatTim'+SUFIX

	degList = np.array([1, 2, 3, 4])
	cutList = np.arange(1, 8)

	if RUNSIMU:

		A1errorList = np.ones((len(degList), len(cutList)))
		R1errorList = np.ones((len(degList), len(cutList)))
		T1timeList = np.ones((len(degList), len(cutList)))

		A2errorList = np.ones((len(degList), len(cutList)))
		R2errorList = np.ones((len(degList), len(cutList)))
		T2timeList = np.ones((len(degList), len(cutList)))
		
		for j, cuts in enumerate(cutList):
			for i, degree in enumerate(degList):

				blockPrint()
				start = time.process_time()
				dirichlet_table = np.ones((2, 2))
				problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensitySquare_inc, 
													dirichlet_table=dirichlet_table, nbel_time=2**4, is1dim=IS1DIM)
				finish = time.process_time()
				T1timeList[i, j] = finish - start
				
				# Error of last step
				A1errorList[i, j], R1errorList[i, j] = problem_inc.normOfError(temp_inc[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperatureSquare_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})
				

				start = time.process_time()
				dirichlet_table = np.ones((2, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensitySquare_spt, 
													dirichlet_table=dirichlet_table, 
													quadArgs={'quadrule':'wq', 'type':2},
													degree_time=degree, nbel_time=2**4+1-degree, is1dim=IS1DIM)
					
				finish = time.process_time()
				T2timeList[i, j] = finish - start

				# Error of last "step"
				newtemp_spt = np.reshape(temp_spt, newshape=(problem_spt.part.nbctrlpts_total, problem_spt.time.nbctrlpts_total), order='F')
				A2errorList[i, j], R2errorList[i, j] = problem_inc.normOfError(newtemp_spt[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperatureSquare_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})

				enablePrint()

			np.savetxt(filenameA1+EXTENSION, A1errorList)
			np.savetxt(filenameR1+EXTENSION, R1errorList)
			np.savetxt(filenameT1+EXTENSION, T1timeList)
			np.savetxt(filenameA2+EXTENSION, A2errorList)
			np.savetxt(filenameR2+EXTENSION, R2errorList)
			np.savetxt(filenameT2+EXTENSION, T2timeList)

	errorList1 = np.loadtxt(filenameA1+EXTENSION)
	errorList2 = np.loadtxt(filenameA2+EXTENSION)
	fig, ax = plt.subplots(figsize=(8, 6))
	for i, degree in enumerate(degList):
		color = COLORLIST[i]
		ax.loglog(2**cutList, errorList2[i, :], color=color, marker=CONFIGLINE0['marker'], markerfacecolor='w',
					markersize=CONFIGLINE0['markersize'], linestyle=CONFIGLINE0['linestyle'], )
		
		ax.loglog(2**cutList, errorList1[i, :], color=color, marker=CONFIGLINE4['marker'], markerfacecolor='w',
					markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'], label='INC-IGA deg. '+str(degree))
		
	ax.loglog([], [], color='k', marker=CONFIGLINE0['marker'], markerfacecolor='w', 
					markersize=CONFIGLINE0['markersize'], linestyle=CONFIGLINE0['linestyle'], label='ST-IGA-GL')

	ax.set_ylabel(r'$L^2$' + ' error')
	ax.set_xlabel('Number of elements by spatial direction')
	ax.set_xlim(left=1, right=200)
	ax.set_ylim(top=1e2, bottom=1e-6)
	ax.legend(loc='upper right')
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE + 'INCL2Convergence' + SUFIX +  '.pdf')
	plt.close(fig)

elif FIG_CASE == 2:

	degree, cuts = 9, 6
	quadArgs = {'quadrule':'iga', 'type':'leg'}
	nbelincList = np.arange(5, 65, 5)
	degsptList  = np.arange(1, 5)
	abserrorInc, relerrorInc = np.ones(len(nbelincList)), np.ones(len(nbelincList))
	abserrorInc2, relerrorInc2 = np.ones(len(nbelincList)), np.ones(len(nbelincList))
	abserrorSpt, relerrorSpt = np.ones((len(degsptList), len(nbelincList))), np.ones((len(degsptList), len(nbelincList)))

	if RUNSIMU:
		for i, nbelsinc in enumerate(nbelincList):

			blockPrint()
			# Incremental
			dirichlet_table = np.ones((2, 2))
			problem_inc, time_inc, temp_inc = simulate_incremental(degree, cuts, powerDensitySquare_inc, 
														dirichlet_table=dirichlet_table,
														nbel_time=nbelsinc, 
														quadArgs={'quadrule':'iga'}, 
														is1dim=IS1DIM)
			
			abserrorInc[i], relerrorInc[i] = problem_inc.normOfError(temp_inc[:, -1], 
														normArgs={'type':'L2',
																'exactFunction':exactTemperatureSquare_inc,
																'exactExtraArgs':{'time':time_inc[-1]}})
			
			np.savetxt(FOLDER2DATA+'2abserrorstag_inc'+SUFIX+EXTENSION, abserrorInc)
			np.savetxt(FOLDER2DATA+'2relerrorstag_inc'+SUFIX+EXTENSION, relerrorInc)

			# Space time
			for j, degspt in enumerate(degsptList):
				dirichlet_table = np.ones((2, 2)); dirichlet_table[-1, 1] = 0
				problem_spt, time_spt, temp_spt = simulate_spacetime(degree, cuts, powerDensitySquare_spt, 
														dirichlet_table=dirichlet_table,
														degree_time=degspt, 
														nbel_time=nbelsinc+1-degspt, 
														quadArgs=quadArgs, is1dim=IS1DIM)
					
				newtemp_spt = np.reshape(temp_spt, newshape=(problem_spt.part.nbctrlpts_total, problem_spt.time.nbctrlpts_total), order='F')
				abserrorSpt[j, i], relerrorSpt[j, i] = problem_inc.normOfError(newtemp_spt[:, -1], 
														normArgs={'type':'L2',
																	'exactFunction':exactTemperatureSquare_inc,
																	'exactExtraArgs':{'time':time_inc[-1]}})

				np.savetxt(FOLDER2DATA+'2abserrorstag_spt'+SUFIX+EXTENSION, abserrorSpt)
				np.savetxt(FOLDER2DATA+'2relerrorstag_spt'+SUFIX+EXTENSION, relerrorSpt)
			
			enablePrint()

	fig, ax = plt.subplots(figsize=(8, 6))
	errorList1 = np.loadtxt(FOLDER2DATA+'2abserrorstag_spt'+SUFIX+EXTENSION)
	for i, deg in enumerate(degsptList):
		ax.loglog(nbelincList+1, errorList1[i, :], color=COLORLIST[i], marker=CONFIGLINE0['marker'], markerfacecolor='w',
					markersize=CONFIGLINE0['markersize'], linestyle=CONFIGLINE0['linestyle'], label='SPT-IGA deg. '+str(int(deg)))
		
	errorList1 = np.loadtxt(FOLDER2DATA+'2abserrorstag_inc'+SUFIX+EXTENSION)
	ax.loglog(nbelincList+1, errorList1, marker=CONFIGLINE4['marker'], markerfacecolor='w', color='k',
					markersize=CONFIGLINE4['markersize'], linestyle=CONFIGLINE4['linestyle'], label='INC-IGA')

	ax.set_ylabel('Stagnation error')
	ax.set_xlabel('Number of control points on time')
	ax.set_xlim(left=4, right=100)
	ax.set_ylim(top=1e1, bottom=1e-9)
	ax.legend(loc='lower left')
	fig.tight_layout()
	fig.savefig(FOLDER2SAVE+'StagnationError'+SUFIX+'.pdf')