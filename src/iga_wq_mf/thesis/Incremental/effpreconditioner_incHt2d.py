from thesis.Incremental.__init__ import *
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
from pysrc.lib.lib_base import sigmoid

def conductivityProperty(args:dict):
	temperature = args.get('temperature')
	Kref  = np.array([[1, 0.5, 0.1],[0.5, 2, 0.25], [0.1, 0.25, 3]])
	Kprop = np.zeros((3, 3, len(temperature)))
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j] 
	for i in range(3): 
		for j in range(3):
			Kprop[i, j, :] = Kref[i, j]*(1.0 + 2.0/(1.0 + np.exp(-5.0*(temperature-1.0))))
	return Kprop 

def capacityProperty(args:dict):
	temperature = args.get('temperature')
	Cprop = np.ones(len(temperature))
	return Cprop

# Set global variables
nbsteps = 26
time_list = np.linspace(0, 0.25, nbsteps) 
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

if RUNSIMU:

	DEGREE, CUTS = 6, 4
	quadArgs = {'quadrule': 'wq', 'type': 2}

	for precond in ITERMETHODS:
		for geoname in geonameList:

			filename = FOLDER2SAVE + 'ResidualHt_' + geoname + '_' + precond + '.dat'  
			
			# blockPrint()
			# Create model 
			geoArgs = {'name': geoname, 'degree': DEGREE*np.ones(3, dtype=int), 
						'nb_refinementByDirection': CUTS*np.ones(3, dtype=int)}
			quadArgs  = {'quadrule': 'wq', 'type': 2}

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
			boundary.add_DirichletConstTemperature(table=np.array([[0, 1], [0, 0], [0, 0]]), temperature=1.0)
		
			problem = heatproblem(material, modelPhy, boundary)
			problem.addSolverConstraints(solverArgs={'preconditioner': precond})

			# Create a Dirichlet condition
			Tinout = np.zeros((modelPhy.nbctrlpts_total, len(time_list)))
			for i in range(1, len(time_list)): Tinout[boundary.thdod, i] = boundary.thDirichletBound[boundary.thdod]

			# Add external force 
			Fend = np.zeros((problem.part.nbctrlpts_total, 1))
			Fext = np.kron(Fend, sigmoid(time_list))

			# Solve
			AllresLin = problem.solveFourierTransientProblem(Tinout=Tinout, 
															Fext_list=Fext, 
															time_list=time_list, 
															alpha=1.0)
			# enablePrint()
		
			np.savetxt(filename, AllresLin)

else:

	for geoname in geonameList:
		fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
		for i, precond in enumerate(ITERMETHODS):
			filename = FOLDER2SAVE + 'ResidualHt_' + geoname + '_' + precond + '.dat'  
			AllresLin = np.loadtxt(filename)
			color = COLORLIST[i+1]

			if precond == "C": labelmethod = 'Classic FD'
			elif precond == "JMC": labelmethod = 'This work'
			elif precond == "TDC": labelmethod = 'Literature'

			stepsMax = int(np.max(AllresLin[:, 0]))
			enum_ax1, enum_ax2 = [], []
			points_ax1, points_ax2 = [], []
			for j in range(1, stepsMax):
				indices = np.where(AllresLin[:, 0]==j)
				newresidue = AllresLin[np.min(indices), 2:]; newresidue = newresidue[newresidue>0]
				enum_ax1.extend(np.arange(len(newresidue))); points_ax1.extend(newresidue)
				newresidue = AllresLin[np.max(indices), 2:]; newresidue = newresidue[newresidue>0]
				enum_ax2.extend(np.arange(len(newresidue))); points_ax2.extend(newresidue)
			
			poly = frompoints2hull(enum_ax1, np.log10(points_ax1), color)
			axs[0].add_patch(poly)
			axs[0].scatter(enum_ax1, np.log10(points_ax1), s=1.5, c=color, alpha=0.2)
			poly = frompoints2hull(enum_ax2, np.log10(points_ax2), color)
			axs[1].add_patch(poly)
			axs[1].scatter(enum_ax2, np.log10(points_ax2), s=1., c=color, alpha=0.2)

			axs[0].plot([], [], marker='s', color=color, label=labelmethod, linewidth=0.5)
			axs[1].plot([], [], marker='s', color=color, label=labelmethod, linewidth=0.5)

		axs[0].set_title('First NR iterations')
		axs[1].set_title('Last NR iterations')
		for ax in axs:
			ax.set_xlim(left=0, right=50)
			ax.set_xlabel('Number of iterations (GMRES)')
			ax.set_ylabel('Log. of relative residue')

		axs[0].legend()
		filename = FOLDER2SAVE + 'HtPreconditioner_' + geoname  + '.pdf'
		fig.tight_layout()
		fig.savefig(filename)
		plt.close(fig=fig)