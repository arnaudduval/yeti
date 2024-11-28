"""
.. Test of mecanical displacement 1D
.. Author: Fabio MADIE
.. Joaquin Cornejo added some corrections 28 nov. 2024
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformOpenCurve
from pysrc.lib.lib_part import part1D
from pysrc.lib.lib_job1d import mechaproblem1D
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_material import mechamat

FORCEVAL = 1e8
MASS, AREA = 380, 0.1
YOUNG, RHO, LENGTH = 210e9, 3800, 1
MECHAMATERIAL = mechamat({'elastic_modulus':YOUNG, 
						'elastic_limit':1e6, 
						'poisson_ratio':0.3,
						'isoHardLaw': {'name':'None'}})
MECHAMATERIAL.addDensity(RHO, isIsotropic=True)
WAVEVEL = np.sqrt(YOUNG/RHO)
COEF = 100 # (20*np.pi/WAVEVEL)**2

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

def dichotomy(upperbound, lowerbound, Nstep, gamma):

	def solveproblem(Nstep_guess, gamma):
		timeList = np.linspace(0, 0.001, Nstep_guess)
		FextList = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
		FextList[-1, :] = FORCEVAL
		displacement = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
		velocity = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
		acceleration = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
		problem.solveNewmarkresolution(displacement, velocity, acceleration, FextList, timeList, gamma=gamma)
		return displacement[-1,-1]

	if ((np.abs(upperbound - Nstep) <= 1) or ((np.abs(lowerbound - Nstep) <= 1))):
		nconvergence = upperbound
	else:
		displacement = solveproblem(Nstep,gamma)
		if np.abs(displacement) > 1 or np.isnan(displacement):
			new_Nstep = (Nstep + upperbound)//2
			nconvergence = dichotomy(upperbound, Nstep, new_Nstep, gamma)
		elif ((np.abs(displacement) < 1 ) and not (np.isnan(displacement))) :
			new_Nstep = (Nstep + lowerbound)//2 + 1
			nconvergence = dichotomy(Nstep, lowerbound, new_Nstep, gamma)	
	return nconvergence

def exactDisplacement(args:dict):
	t = args['time']
	x = args['position']
	uf = FORCEVAL*x/(YOUNG*AREA)
	for i in range(1, 201):
		factor = -8*FORCEVAL*LENGTH/(((np.pi)**2)*YOUNG*AREA)
		uf += factor*((((-1)**(i-1))/((2*i-1)**2))*np.sin((2*i-1)*np.pi*x/(2*LENGTH))*(np.cos((2*i-1)*np.pi*WAVEVEL*t/(2*LENGTH))))
	return uf

def simulate(curve, ntime, quadArgs):
	# Create geometry
	modelIGA = part1D(curve, kwargs={'quadArgs': quadArgs})
	timeList = np.linspace(0, 1e-3, ntime)

	# Create boundary condition
	boundary = boundaryCondition(modelIGA.nbctrlpts)
	boundary.add_DirichletConstTemperature(table=np.array([[1, 0]]))

	# Set mecanical problem
	problem = mechaproblem1D(MECHAMATERIAL, modelIGA, boundary)

	# Create external force	
	FextList = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
	FextList[-1, 1:] = FORCEVAL
	
    #Solve problem
	displacement = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
	velocity = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))
	acceleration = np.zeros((modelIGA.nbctrlpts_total, len(timeList)))	
	problem.solveNewmarkresolution(displacement, velocity, acceleration, FextList, timeList)

	return problem, displacement, timeList

# Convergence test
label_list = ['Linear', 'Quadratic C0', 'Quadratic C1']
pourcentage_list = [0.05, 0.1, 0.2, 0.5, 1.]
cutsList = np.arange(1, 5)

AbserrorList = np.zeros(len(cutsList))
RelerrorList = np.zeros(len(cutsList))

GLplot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
Lobplot = {'marker': 'D', 'linestyle': '-.', 'markersize': 4}
WQ1plot = {'marker': 'o', 'linestyle': '--', 'markersize': 4}
WQ2plot = {'marker': 'x', 'linestyle': ':', 'markersize': 4}
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

for k, pourcentage in enumerate(pourcentage_list):
	fig, ax = plt.subplots(figsize=(9, 6))
	name = 'Erreur_pascritique_' + str(pourcentage) + str('pourcent.png')
	figname = folder + name
	for quadrule, quadtype, plotpars in zip(['iga', 'wq'], ['leg', 2], [GLplot, WQ2plot]):
		quadArgs = {'quadrule': quadrule, 'type': quadtype}
		totaldof_list = np.zeros((len(cutsList), 4))
		for j, cut in enumerate(cutsList):
			nel = 2**cut
			curve_list = [createUniformOpenCurve(1, nel, 1), createUniformOpenCurve(2, nel, 1, multiplicity=2), createUniformOpenCurve(2, nel, 1)]
			for n, curve in enumerate(curve_list):
				print(str(curve.degree) +' '+ str(nel)+' '+ str(pourcentage))		
				modelIGA = part1D(curve, kwargs={'quadArgs': quadArgs}) 

                # Create boundary condition
				boundary = boundaryCondition(modelIGA.nbctrlpts)
				boundary.add_DirichletConstTemperature(table=np.array([[1 , 0]]))

                # Set mechanical problem
				problem = mechaproblem1D(MECHAMATERIAL, modelIGA, boundary)
				start_upperbound = 6000
				start_lowerbound = 3
				if start_lowerbound < 3: start_lowerbound = 3
				Nstep = (start_upperbound + start_lowerbound) // 2 #(start_hboundery + start_lboundery)//2
				nconvergence = dichotomy(start_upperbound, start_lowerbound, Nstep, 0.5)
				print('Nstep min = ' + str(nconvergence))
				n_conv = round(nconvergence * (1 + pourcentage))
				print('Ntime = ' + str(n_conv))
				problem, displacement, timeList = simulate(curve, n_conv, quadArgs)
				totaldof_list[j, n] = problem.part.nbctrlpts_total
				print(problem.part.nbctrlpts_total)
				AbserrorList[j], RelerrorList[j] = problem.normOfError(displacement[:, -1], 
															normArgs={
																	'exactFunction':exactDisplacement,
																	'exactExtraArgs':{'time': 1e-3}
																	}
																	)

		for i in range(len(curve_list)):
			ax.loglog(totaldof_list[:, i], AbserrorList, color=cycle[i], marker=plotpars['marker'])
		fig.savefig(figname)
		
	for i in range(len(curve_list)):
		ax.loglog([], [], label=label_list[i],color=cycle[i])
	ax.loglog([], [], color='k', marker=GLplot['marker'], label="IGA-GL")      
	ax.loglog([], [], color='k', marker=WQ2plot['marker'], label="IGA-WQ 2")
	ax.set_ylabel(r'$||u-u^h||_{L^2(\Omega)}$')
	ax.set_xlabel('Number of elements')
	ax.legend()
	fig.tight_layout()
	fig.savefig(figname)