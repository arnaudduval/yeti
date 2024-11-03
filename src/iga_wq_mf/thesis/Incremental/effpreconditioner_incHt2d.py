from thesis.Incremental.__init__ import *
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
from pysrc.lib.lib_base import sigmoid, vtk2png

folder = FOLDER2SAVE + '/animation_heat/'
if not os.path.isdir(folder): os.mkdir(folder)

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((3, 3, len(temperature)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j] 
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j]*(3.0 + 2.0*np.tanh(temperature/50))
	return Kprop 

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	Cprop = np.ones(len(temperature))
	return Cprop

# Set global variables
NBSTEPS = 65
TIMELIST = np.linspace(0., 1., NBSTEPS) 
geonameList = ['cb', 'vb']
ITERMETHODS = ['C', 'TDC', 'JMC']
RUNSIMU = False

def frompoints2hull(a, b, color, factor=1.0):
	points = np.vstack((a, b)).T
	hull = ConvexHull(points)
	cent = np.mean(points, axis=0)
	pts = points[hull.vertices]
	poly = Polygon(factor*(pts - cent) + cent, closed=True,
			capstyle='round', facecolor=color, 
			edgecolor=color, linewidth=2, alpha=0.5)
	return poly

def simulation_heattransfer(geoname='cb', preconditioner='TDC'):
	# Create model 
	geoArgs = {'name': geoname, 'degree': DEGREE*np.ones(3, dtype=int), 
				'nb_refinementByDirection': CUTS*np.ones(3, dtype=int)}
	quadArgs  = {'quadrule': 'wq', 'type': 1}

	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)

	# Add material 
	material = heatmat()
	material.addConductivity(conductivityProperty, isIsotropic=False) 
	material.addCapacity(capacityProperty, isIsotropic=False) 

	# Block boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	boundary.add_DirichletConstTemperature(table=np.array([[1, 0], [0, 0], [0, 0]]))
	boundary.add_DirichletConstTemperature(table=np.array([[0, 1], [0, 0], [0, 0]]), temperature=10.0)

	problem = heatproblem(material, modelPhy, boundary)
	problem.addSolverConstraints(solverArgs={'preconditioner': preconditioner})
	problem._thresLin = 1e-12

	# Create a Dirichlet condition
	Tinout = np.zeros((modelPhy.nbctrlpts_total, len(TIMELIST)))
	for i in range(1, len(TIMELIST)): Tinout[boundary.thdod, i] = TIMELIST[i]/TIMELIST[-1]*boundary.thDirichletBound[boundary.thdod]

	# Add external force 
	Fend = np.zeros((problem.part.nbctrlpts_total, 1))
	Fext = np.kron(Fend, sigmoid(TIMELIST))

	# Solve
	AllresLin = problem.solveFourierTransientProblem(Tinout=Tinout, 
													Fext_list=Fext, 
													time_list=TIMELIST, 
													alpha=1.0)
	return problem, Tinout, AllresLin

if RUNSIMU:

	DEGREE, CUTS = 3, 4
	quadArgs = {'quadrule': 'wq', 'type': 1}

	problem, Tinout, _ = simulation_heattransfer()
	for k, i in enumerate(range(0, np.size(Tinout, axis=1), 4)):
		problem.part.postProcessingPrimal(fields={'temp':np.atleast_2d(Tinout[:, i])}, 
										name='out_'+str(k), folder=folder)
		
	run(folder=folder, filename='out_', nbFiles=k)

	for i in range(17): vtk2png(folder, filename='out_'+str(i), title='Temperature', clim=[0, 10], n_colors=21, camera_position='xz')

	for precond in ITERMETHODS:
		for geoname in geonameList:

			filename = FOLDER2DATA + 'ResidualHt_' + geoname + '_' + precond + '.dat'  
			
			# blockPrint()
			AllresLin = simulation_heattransfer(geoname, precond)[-1]
			# enablePrint()
			np.savetxt(filename, AllresLin)

# for geoname in geonameList:
# 	fig1, ax1 = plt.subplots()
# 	fig2, ax2 = plt.subplots()
# 	axs = [ax1, ax2]
# 	for i, precond in enumerate(ITERMETHODS):
# 		filename = FOLDER2DATA + 'ResidualHt_' + geoname + '_' + precond + '.dat'  
# 		AllresLin = np.loadtxt(filename)
# 		color = COLORLIST[i+1]

# 		if precond == "C": labelmethod = 'Classic FD'
# 		elif precond == "JMC": labelmethod = 'This work'
# 		elif precond == "TDC": labelmethod = 'Literature'

# 		stepsMax = int(np.max(AllresLin[:, 0]))
# 		enum_ax1, enum_ax2 = [], []
# 		points_ax1, points_ax2 = [], []
# 		for j in range(1, stepsMax):
# 			indices = np.where(AllresLin[:, 0]==j)
# 			newresidue = AllresLin[np.min(indices), 2:]; newresidue = newresidue[newresidue>0]
# 			enum_ax1.extend(np.arange(len(newresidue))); points_ax1.extend(newresidue)
# 			newresidue = AllresLin[np.max(indices)-1, 2:]; newresidue = newresidue[newresidue>0]
# 			enum_ax2.extend(np.arange(len(newresidue))); points_ax2.extend(newresidue)
		
# 		poly = frompoints2hull(enum_ax1, np.log10(points_ax1), color)
# 		axs[0].add_patch(poly)
# 		poly = frompoints2hull(enum_ax2, np.log10(points_ax2), color)
# 		axs[1].add_patch(poly)

# 		axs[0].plot([], [], marker='s', color=color, label=labelmethod, linewidth=0.5)
# 		axs[1].plot([], [], marker='s', color=color, label=labelmethod, linewidth=0.5)

# 	for ax in axs:
# 		ax.set_xlim(left=0, right=50)
# 		ax.set_ylim(top=0, bottom=-12)
# 		ax.set_xlabel('Number of iterations (GMRES)')
# 		ax.set_ylabel('Log. of relative residue')

# 	axs[0].legend()

# 	for fig, sufix in zip([fig1, fig2], ['first', 'last']):
# 		filename = FOLDER2SAVE + 'HtPrecond_' + sufix + '_' + geoname  + '.pdf'
# 		fig.tight_layout()
# 		fig.savefig(filename)
# 		plt.close(fig=fig)

for geoname in geonameList:
	fig, ax = plt.subplots()
	for i, precond in enumerate(ITERMETHODS):
		filename = FOLDER2DATA + 'ResidualHt_' + geoname + '_' + precond + '.dat'  
		AllresLin = np.loadtxt(filename)
		color = COLORLIST[i+1]

		if precond == "C": labelmethod = 'Standard Fast Diag.'
		elif precond == "JMC": labelmethod = 'My contribution'
		elif precond == "TDC": labelmethod = 'Proposition in literature'

		stepsMax = int(np.max(AllresLin[:, 0]))
		enum_ax = []
		points_ax = []
		for j in range(1, stepsMax):
			indices = np.where(AllresLin[:, 0]==j)
			newresidue = AllresLin[np.min(indices), 2:]; newresidue = newresidue[newresidue>0]
			enum_ax.extend(np.arange(len(newresidue))); points_ax.extend(newresidue)
			newresidue = AllresLin[np.max(indices)-1, 2:]; newresidue = newresidue[newresidue>0]
			enum_ax.extend(np.arange(len(newresidue))); points_ax.extend(newresidue)
		
		poly = frompoints2hull(enum_ax, np.log10(points_ax), color)
		ax.add_patch(poly)

		ax.plot([], [], marker='s', color=color, label=labelmethod, linewidth=0.5)


	ax.set_xlim(left=0, right=50)
	ax.set_ylim(top=0, bottom=-12)
	ax.set_xlabel('Number of iterations (GMRES)')
	ax.set_ylabel('Log. of relative residue')
	# ax.legend()
	ax.legend(ncol=2, bbox_to_anchor=(0.5, 1.2), loc='upper center')

	filename = FOLDER2SAVE + 'HtPrecond_' + 'All' + '_' + geoname  + '.png'
	fig.tight_layout()
	fig.savefig(filename)
	plt.close(fig=fig)