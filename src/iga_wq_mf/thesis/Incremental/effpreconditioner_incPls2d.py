from thesis.Incremental.__init__ import *
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon

# Set global variables
ITERMETHODS = ['C', 'JMC', 'TDC']
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

	DEGREE, CUTS = 2, 5
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}

	for precond in ITERMETHODS:
		filename = FOLDER2DATA + 'Residualpls_' + precond + '.dat'        
		AllresLin = simulate_2d(DEGREE, CUTS, quadArgs, precond=precond)[2]
		np.savetxt(filename, AllresLin)

else:

	fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
	for i, precond in enumerate(ITERMETHODS):
		filename = FOLDER2DATA + 'Residualpls_' + precond + '.dat' 
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
	filename = FOLDER2SAVE + 'PlsPreconditioner' + '.pdf'
	fig.tight_layout()
	fig.savefig(filename)