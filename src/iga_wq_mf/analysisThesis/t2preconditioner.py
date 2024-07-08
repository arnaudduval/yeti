from pysrc.lib.__init__ import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
from pysrc.lib.lib_material import mechamat
from pysrc.lib.lib_boundary import boundaryCondition
from pysrc.lib.lib_job import mechaproblem

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/'
if not os.path.isdir(folder): os.mkdir(folder)

# Set global variables
GEONAME = 'SQ'
TRACTION = 400.0
YOUNG, POISSON = 2500, 0.25
NBSTEPS = 101
TIME_LIST = np.linspace(0, np.pi/2, NBSTEPS)
MATARGS = {'elastic_modulus':YOUNG, 'elastic_limit':5, 'poisson_ratio': POISSON, 
			'isoHardLaw': {'name':'linear', 'Eiso':0.0}, 
			'kineHardLaw':{'parameters':np.array([[500, 0]])}
			}
ITERMETHODS = ['C', 'JMC', 'TDC']
runSimu = False

def forceSurf_infPlate(P:list):
	x = P[0, :]; nnz = np.size(P, axis=1)
	tmp = np.zeros((2, nnz)); tmp[1, :] = x**2-1/4
	F = np.zeros((2, nnz))
	F[1, :] = -TRACTION*(np.min(tmp, axis=0))**2
	return F

def simulate(degree, cuts, quadArgs, precond='JMC'):
	geoArgs = {'name': 'SQ', 'degree': degree*np.ones(3, dtype=int), 
				'nb_refinementByDirection': cuts*np.ones(3, dtype=int), 
				'extra':{'XY':np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])}
			}
	blockPrint()
	material = mechamat(MATARGS)
	modelGeo = Geomdl(geoArgs)
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs=quadArgs)
	meshparam = modelPhy.compute_global_mesh_parameter()

	# Set Dirichlet boundaries
	boundary = boundaryCondition(modelPhy.nbctrlpts)
	table = np.zeros((2, 2, 2), dtype=int); table[1, 0, :] = 1
	boundary.add_DirichletDisplacement(table=table)
	enablePrint()

	# Solve elastic problem
	problem = mechaproblem(material, modelPhy, boundary)
	problem._thresNL = 1e-6; problem._itersNL = 100; problem._linPreCond = precond
	Fref = problem.compute_surfForce(forceSurf_infPlate, nbFacePosition=3)[0]
	Fext_list = np.zeros((2, modelPhy.nbctrlpts_total, NBSTEPS))
	for k in range(len(TIME_LIST)): Fext_list[:, :, k] = np.sin(TIME_LIST[k])*Fref
	displacement = np.zeros(np.shape(Fext_list))
	AllresLin, internalVars = problem.solveElastoPlasticityProblem(displacement, Fext_list)
	return problem, AllresLin, displacement, meshparam, internalVars

if runSimu:

	DEGREE, CUTS = 2, 5
	quadArgs = {'quadrule': 'iga', 'type': 'leg'}

	for PCGmethod in ITERMETHODS:
		filename = folder + 'ResPCGpls_' + GEONAME + '_' + PCGmethod + '.dat'        
		problem, AllresLin, displacement, meshparam, _= simulate(DEGREE, CUTS, quadArgs, precond=PCGmethod)
		np.savetxt(filename, AllresLin)

else:

	from scipy.spatial import ConvexHull
	from matplotlib.patches import Polygon

	def frompoints2hull(a, b, color):
		points = np.vstack((a, b)).T
		hull = ConvexHull(points)
		cent = np.mean(points, axis=0)
		pts = points[hull.vertices]
		poly = Polygon(factor*(pts - cent) + cent, closed=True,
				capstyle='round', facecolor=color, alpha=0.5)
		return poly

	factor = 1.00
	fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
	for i, PCGmethod in enumerate(ITERMETHODS):
		filename = folder + 'ResPCGpls_' + GEONAME + '_' + PCGmethod + '.dat'
		AllresLin   = np.loadtxt(filename)

		if PCGmethod == "C": labelmethod = 'Classic FD\nmethod'
		elif PCGmethod == "JMC": labelmethod = 'This work'
		elif PCGmethod == "TDC": labelmethod = 'Literature'

		nbstepssimu = int(np.max(AllresLin[:, 0]))
		enum_ax1, enum_ax2 = [], []
		points_ax1, points_ax2 = [], []
		for j in range(2, nbstepssimu):
			indices = np.where(AllresLin[:, 0]==j)
			newresidue = AllresLin[np.min(indices), 2:]; newresidue = newresidue[newresidue>0]
			enum_ax1.extend(np.arange(len(newresidue))); points_ax1.extend(newresidue)
			newresidue = AllresLin[np.max(indices), 2:]; newresidue = newresidue[newresidue>0]
			enum_ax2.extend(np.arange(len(newresidue))); points_ax2.extend(newresidue)
		
		poly = frompoints2hull(enum_ax1, np.log10(points_ax1), COLORLIST[i+1])
		axs[0].add_patch(poly)
		poly = frompoints2hull(enum_ax2, np.log10(points_ax2), COLORLIST[i+1])
		axs[1].add_patch(poly)

		axs[0].plot([], [], marker='s', color=COLORLIST[i+1], label=labelmethod, linewidth=0.5)
		axs[1].plot([], [], marker='s', color=COLORLIST[i+1], label=labelmethod, linewidth=0.5)

	# fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 5))
	# for i, PCGmethod in enumerate(ITERMETHODS):
	# 	filename = folder + 'ResPCGpls_' + GEONAME + '_' + PCGmethod + '.dat'
	# 	AllresLin   = np.loadtxt(filename)

	# 	if PCGmethod == "C": labelmethod = 'Classic FD\nmethod'
	# 	elif PCGmethod == "JMC": labelmethod = 'This work'
	# 	elif PCGmethod == "TDC": labelmethod = 'Literature'
		
	# 	maxStep = int(np.max(AllresLin[:, 0]))
	# 	for j in range(2, maxStep):
	# 		ind = np.where(AllresLin[:, 0]==j)

	# 		newresidue = AllresLin[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 		axs[0].semilogy(np.arange(len(newresidue)), newresidue, 
	# 					color=COLORLIST[i+1], linewidth=0.5)
			
	# 		newresidue = AllresLin[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 		axs[1].semilogy(np.arange(len(newresidue)), newresidue, 
	# 					color=COLORLIST[i+1], linewidth=0.5)
			
	# 	ind = np.where(AllresLin[:, 0]==1)
	# 	newresidue = AllresLin[np.min(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 	axs[0].semilogy(np.arange(len(newresidue)), newresidue, marker='s', color=COLORLIST[i+1], label=labelmethod, linewidth=0.5)
		
	# 	newresidue = AllresLin[np.max(ind), 2:]; newresidue = newresidue[newresidue>0]
	# 	axs[1].semilogy(np.arange(len(newresidue)), newresidue, marker='s', color=COLORLIST[i+1], label=labelmethod, linewidth=0.5)

	axs[0].set_title('First NR iterations')
	axs[1].set_title('Last NR iterations')
	for ax in axs:
		ax.set_xlim(left=0, right=50)
		ax.set_xlabel('Number of iterations of iterative solver')
		ax.set_ylabel('Log. relative residue')
		# ax.set_ybound(lower=1e-8, upper=10)

	axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))
	filename = folder + 'PlasticityPreconditioner' + '.png'
	fig.tight_layout()
	fig.savefig(filename)