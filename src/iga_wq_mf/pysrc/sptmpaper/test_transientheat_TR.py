from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import createUniformCurve
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part, part1D
from pysrc.lib.lib_material import heatmat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import heatproblem
from pysrc.lib.lib_stjob import stheatproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/transient/'
if not os.path.isdir(folder): os.mkdir(folder)

extension = '.dat'
FIG_CASE  = 3
DATAEXIST = False
ISLINEAR  = False
c = 0.001

def conductivityProperty(args:dict):
	temperature = args['temperature']
	Kref  = np.array([[1., 0.5],[0.5, 2.0]])
	Kprop = np.zeros((2, 2, len(temperature)))
	for i in range(2): 
		for j in range(2):
			if ISLINEAR: Kprop[i, j, :] = Kref[i, j]
			else: Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0*np.exp(-np.abs(temperature)))
	return Kprop 

def capacityProperty(args:dict):
	temperature = args['temperature']
	if ISLINEAR: Cprop = np.ones(len(temperature))
	else: Cprop = (1.0 + np.exp(-np.abs(temperature)))
	return Cprop

def exactTemperature_inc(args:dict):
	qpPhy = args.get('position')
	t = args.get('time')
	x = qpPhy[0, :]; y = qpPhy[1, :]
	u = c*(-5*x + 6*y + 45)*(5*x + 6*y - 45)*np.sin(np.pi*x)*np.sin(np.pi*t)
	return u

def exactTemperature_st(qpPhy):
	x = qpPhy[0, :]; y = qpPhy[1, :]; t = qpPhy[2, :]
	u = c*(-5*x + 6*y + 45)*(5*x + 6*y - 45)*np.sin(np.pi*x)*np.sin(np.pi*t)
	return u

def powerDensity(args:dict):
	position = args['position']; t = args['time']
	x = position[0, :]; y = position[1, :]
	nc_sp = np.size(position, axis=1); f = np.zeros(nc_sp)
	if ISLINEAR:
		f = (4*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(5*x + 6*y - 45) 
					- 16*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45) 
					- 94*c*np.sin(np.pi*t)*np.sin(np.pi*x) 
					+ c*np.pi*np.cos(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45)*(5*x + 6*y - 45) 
					+ c*np.pi**2*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
		)
	else: 
		u = c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
		f = (
			(2*np.exp(-np.abs(u)) + 1)*(
				50*c*np.sin(np.pi*t)*np.sin(np.pi*x) 
				- 10*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45) 
				+ 10*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(5*x + 6*y - 45) 
				+ c*np.pi**2*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45)*(5*x + 6*y - 45)) 
			- 2*(np.exp(-np.abs(u)) + 1/2)*(
				6*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45) 
				+ 6*c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(5*x + 6*y - 45)) 
			+ 2*np.exp(-np.abs(u))*np.sign(u)*(
				5*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45) 
				- 5*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(5*x + 6*y - 45) 
				+ c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45)*(5*x + 6*y - 45))**2 
			+ 4*np.exp(-np.abs(u))*np.sign(u)*(
				6*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45) 
				+ 6*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(5*x + 6*y - 45))**2 
			+ 2*np.exp(-np.abs(u))*np.sign(u)*(
				6*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45) 
				+ 6*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(5*x + 6*y - 45))*(
					5*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(6*y - 5*x + 45) 
					- 5*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(5*x + 6*y - 45) 
					+ c*np.pi*np.cos(np.pi*x)*np.sin(np.pi*t)*(6*y - 5*x + 45)*(5*x + 6*y - 45)) 
			- 72*c*np.sin(np.pi*t)*np.sin(np.pi*x)*(4*np.exp(-np.abs(u)) + 2) 
			+ c*np.pi*np.cos(np.pi*t)*np.sin(np.pi*x)*(np.exp(-np.abs(u)) + 1)*(6*y - 5*x + 45)*(5*x + 6*y - 45)
		)
	return f

def simulate(degree, cuts, quadArgs, cuts_time=None):
	geoArgs = {'name': 'TP', 'degree': degree*np.ones(3, dtype=int), 
						'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	if cuts_time is None: cuts_time = np.copy(cuts)
	timespan  = 1
	nbsteps   = 2**cuts_time
	time_inc  = np.linspace(0, timespan, nbsteps+1) 
	time_crv  = part1D(createUniformCurve(1, nbsteps, timespan), {'quadArgs': quadArgs})

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	boundary_inc = boundaryCondition(modelPhy.nbctrlpts)
	boundary_inc.add_DirichletConstTemperature(table=np.ones((2, 2)))

	dirichlet_table = np.ones((3, 2)); dirichlet_table[-1, 1] = 0
	stnbctrlpts = np.array([*modelPhy.nbctrlpts[:modelPhy.dim], time_crv.nbctrlpts])
	boundary_st = boundaryCondition(stnbctrlpts)
	boundary_st.add_DirichletConstTemperature(table=dirichlet_table)

	# Transient model
	problem_inc = heatproblem(material, modelPhy, boundary_inc)
	problem_st  = stheatproblem(material, modelPhy, time_crv, boundary_st)

	# Add external force 
	Fext_list = np.zeros((problem_inc.part.nbctrlpts_total, len(time_inc)))
	for i, t in enumerate(time_inc):
		Fext_list[:, i] = problem_inc.compute_volForce(powerDensity, args={'position':problem_inc.part.qpPhy, 'time':t})

	# Solve
	Tinout = np.zeros((modelPhy.nbctrlpts_total, len(time_inc)))
	problem_inc.solveFourierTransientProblem(Tinout=Tinout, Fext_list=Fext_list, time_list=time_inc, alpha=0.5)

	return problem_inc, problem_st, Tinout

if not DATAEXIST:
	if FIG_CASE == 1:
		lastsufix = 'linear' if ISLINEAR else 'nonlin'
		degree_list = np.array([1, 2, 3, 4])
		cuts_list   = np.arange(1, 6)
		for quadrule, quadtype in zip(['iga'], ['leg']):
			sufix = '_' + quadrule + '_' + quadtype + '_' + lastsufix
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			L2errorTable = np.zeros((len(degree_list)+1, len(cuts_list)+1))
			L2relerrorTable = np.zeros((len(degree_list)+1, len(cuts_list)+1))
			L2errorTable[0, 1:] = cuts_list; L2relerrorTable[0, 1:] = cuts_list
			L2errorTable[1:, 0] = degree_list; L2relerrorTable[1:, 0] = degree_list
			filename1 = folder+'L2error_meshpar'+sufix+extension
			filename2 = folder+'L2relerror_meshpar'+sufix+extension
			for j, cuts in enumerate(cuts_list):
				for i, degree in enumerate(degree_list):
					nbels = 2**cuts_list
					problem_inc, problem_st, output = simulate(degree, cuts, quadArgs)
					L2errorTable[i+1, j+1], L2relerrorTable[i+1, j+1] = problem_st.normOfError(np.ravel(output, order='F'), 
																								normArgs={'type':'L2', 
																								'exactFunction':exactTemperature_st},)
					np.savetxt(filename1, L2errorTable)
					np.savetxt(filename2, L2relerrorTable)

	elif FIG_CASE == 3:
		lastsufix = 'linear' if ISLINEAR else 'nonlin'
		degree_list = np.array([1, 2, 3, 4])
		cuts_list   = np.arange(1, 6)
		for quadrule, quadtype in zip(['iga'], ['leg']):
			quadArgs = {'quadrule': quadrule, 'type': quadtype}
			error_list = np.ones((len(degree_list), len(cuts_list), 2**np.max(cuts_list)))
			for j, cuts in enumerate(cuts_list):
				for i, degree in enumerate(degree_list):
					nbels = 2**cuts_list
					problem_inc, problem_st, output = simulate(degree, cuts, quadArgs, cuts_time=np.max(cuts_list))
					for k, step in enumerate(problem_st.time.ctrlpts[1:-1]):
						_, error_list[i, j, k] = problem_inc.normOfError(output[:, k+1], 
																		normArgs={'type':'L2',
																				'exactFunction':exactTemperature_inc,
																				'exactExtraArgs':{'time':step}})
			np.save(folder + 'incrementalheat', error_list)
			error_list = np.load(folder + 'incrementalheat.npy')

			for k in range(np.size(error_list, axis=2)):
				fig, ax = plt.subplots(figsize=(9, 6))
				for i, degree in enumerate(degree_list):
					color = COLORLIST[i]
					ax.loglog(2**cuts_list, error_list[i, :, k], color=color, marker='o', markerfacecolor='w',
								markersize=10, linestyle='-', label='degree ' + r'$p=\,$' + str(degree))
				ax.set_ylabel(r'$\displaystyle\frac{||u - u^h||_{L_2(\Omega)}}{||u||_{L_2(\Omega)}}$')
				ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
				ax.set_ylim(top=1e2, bottom=1e-4)
				ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
				fig.tight_layout()
				fig.savefig(folder + 'steps1/FigConvergenceIncrHeat' + str(k+1) +  '.pdf')
				plt.close(fig)

else:

	if FIG_CASE == 1:
		normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
		onlyMarker1 = {'marker': 'o', 'linestyle': '--', 'markersize': 6}
		onlyMarker2 = {'marker': 'x', 'linestyle': ':', 'markersize': 6}
		plotoptions = [normalPlot, onlyMarker1, onlyMarker2]

		lastsufix = 'linear' if ISLINEAR else 'nonlin'
		figname = folder + 'SPTNonLinearConvergenceL2'+lastsufix+'.pdf'
		filenames = ['L2relerror_meshpar_iga_leg_']

		normalPlot  = {'marker': 's', 'linestyle': '-', 'markersize': 10}
		fig, ax = plt.subplots(figsize=(8, 6))

		for filename, plotops in zip(filenames, plotoptions):
			quadrule = filename.split('_')[2]
			table = np.loadtxt(folder+filename+lastsufix+extension)	
			nbels   = 2**(table[0, 1:])
			degrees = table[1:, 0]
			errors  = table[1:, 1:]
			for i, degree in enumerate(degrees):
				color = COLORLIST[i]
				if quadrule == 'iga': 
					ax.loglog(nbels, errors[i, :], label='IGA-GL deg. '+str(int(degree)), color=color, marker=plotops['marker'],
								markerfacecolor='w', markersize=plotops['markersize'], linestyle=plotops['linestyle'])
					
					# slope = np.polyfit(np.log10(nbels[2:]),np.log10(errors[i, 2:]), 1)[0]
					# slope = round(slope, 1)
					# annotation.slope_marker((nbels[-1], errors[i, -1]), slope, 
					# 				poly_kwargs={'facecolor': (0.73, 0.8, 1)}, ax=ax)			
				else: 
					ax.loglog(nbels, errors[i, :], color=color, marker=plotops['marker'], markerfacecolor='w',
							markersize=plotops['markersize'], linestyle=plotops['linestyle'])
						
				fig.savefig(figname)

		# ax.loglog([], [], color='k', marker=onlyMarker1['marker'], markerfacecolor='w',
		# 				markersize=onlyMarker1['markersize'], linestyle=onlyMarker1['linestyle'], label="IGA-WQ 2")
		# ax.loglog([], [], color='k', marker=onlyMarker2['marker'], markerfacecolor='w',
		# 		markersize=onlyMarker2['markersize'], linestyle=onlyMarker2['linestyle'], label="IGA-WQ 4")

		ax.set_ylabel(r'$\displaystyle ||u - u^h||_{L^2(\Pi)}/||u||_{L^2(\Pi)}$')
		ax.set_xlabel('Mesh discretization ' + r'$h^{-1}$')
		ax.set_xlim(left=1, right=100)
		ax.set_ylim(top=1e1, bottom=1e-7)
		ax.legend(loc='lower left')
		fig.tight_layout()
		fig.savefig(figname)