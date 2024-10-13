from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import (createUniformKnotvector_Rmultiplicity, evalDersBasisCSRPy, cropImage)
from pysrc.lib.lib_quadrules import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part
import imageio

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/videos/'

def plot2DGeo(model:part, sampleSize=101):
	"Plots a 2D geometry "

	def plot_mesh(pts, shape, ax):
		"Plots mesh of control points"

		if pts.shape[0] == 3: pts = pts[:2, :]

		pts2D = []
		for j in range(shape[1]):
			pts2D_temp = []
			for i in range(shape[0]):
				pos = i + j*shape[0]
				pts2D_temp.append(pts[:, pos].tolist())
			pts2D.append(pts2D_temp)
		pts2D = np.asarray(pts2D)

		# In the first direction
		for _ in range(shape[1]): 
			x = pts2D[_, :, 0]; y = pts2D[_, :, 1]
			ax.plot(x, y, color=COLORLIST[1], linestyle='--')

		# In the second direction
		for _ in range(shape[0]):
			x = pts2D[:, _, 0]; y = pts2D[:, _, 1]
			ax.plot(x, y, color=COLORLIST[1], linestyle='--')

		return

	ctrlpts = model.ctrlpts
	evalpts = model.interpolateMeshgridField(sampleSize=sampleSize)[0]

	# Get basis using fortran
	basis, indices = [], []
	for i in range(model.dim):
		dersb, indi, indj = evalDersBasisCSRFortran(model.degree[i], model.knotvector[i], np.unique(model.knotvector[i]))
		basis.append(dersb); indices.append(indi); indices.append(indj)

	# Get position and determinant 
	nbknots = [len(np.unique(model.knotvector[i])) for i in range(model.dim)]
	inpts = [*nbknots, *indices, *basis, model.ctrlpts]
	knotsPhy = geophy.interpolate_meshgrid_2d(*inpts)

	X = np.asarray(evalpts[0, :].reshape((sampleSize, sampleSize)).tolist())
	Y = np.asarray(evalpts[1, :].reshape((sampleSize, sampleSize)).tolist())
	Z = np.zeros(X.shape)

	fig, ax = plt.subplots(figsize=(5, 5))
	ax.grid(None)
	plot_mesh(ctrlpts, model.nbctrlpts, ax)
	ax.pcolormesh(X, Y, Z, cmap=plt.cm.Paired, shading='gouraud')
	ax.plot([], [], label='B-Spline surface')
	ax.plot(ctrlpts[0, :], ctrlpts[1, :], 'o', label='Control points net')
	ax.plot(knotsPhy[0, :], knotsPhy[1, :], color='k', marker='s', linestyle='', label='Knots')
	
	ax.set_xticks(np.arange(0, np.ceil(max(evalpts[:, 0]))+2, .5))
	ax.set_yticks(np.arange(0, np.ceil(max(evalpts[:, 1]))+1, .5))
	ax.set_xlim(left=-0.1, right=1.25)
	ax.set_ylim(bottom=-0.1, top=1.25)

	# ax.legend(loc='lower left')
	ax.legend(loc='upper right')
	ax.set_aspect('equal', adjustable='box')
	ax.set_xlabel(r'$x_1$')
	ax.set_ylabel(r'$x_2$')
	fig.tight_layout()
	
	return fig


# Set global variables
CASE      = 0
extension = '.png'

if CASE == 0: # B-spline curve

	def case0(folder, extension, yoffset=0):
		# Set filename
		filename = folder + 'BSplinecurve' + extension

		# Create the curve 
		crv            = BSpline.Curve()
		crv.degree     = 3
		crv.ctrlpts    = [[-1, 1, 0], [-0.5, 0.25, 0], [0, 2+yoffset, 0], [0.5, 1., 0.],
							[0.75, -0.5, 0], [1.5, 1, 0], [2, 0, 0]]
		crv.knotvector = np.array([0., 0., 0., 0., 0.25, 0.75, 0.75, 1., 1., 1., 1.])

		# Get data
		evalpts = np.asarray(crv.evalpts); ctrlpts = np.asarray(crv.ctrlpts)
		basisKnots = evalDersBasisDensePy(crv.degree, crv.knotvector, np.unique(crv.knotvector))
		knotsPhy = basisKnots[0].T @ np.array(crv.ctrlpts)
		evalpts = np.asarray(crv.evalpts); ctrlpts = np.asarray(crv.ctrlpts)
		
		fig, ax = plt.subplots()
		ax.plot(evalpts[:, 0], evalpts[:, 1], label='B-Spline curve')
		ax.plot(ctrlpts[:, 0], ctrlpts[:, 1], 'o--', markersize=10, label='Control points net')
		ax.plot(knotsPhy[:, 0], knotsPhy[:, 1], color='k', marker='s', linestyle='', label='Knots')
		ax.set_xlim([-1.5, 2.5]); ax.set_ylim([-1., 2.5])
		ax.set_xlabel(r'$x_1$'); ax.set_ylabel(r'$x_2$') 
		ax.legend(loc='upper right')
		fig.tight_layout()
		fig.savefig(filename, dpi=300)
		return filename

	images = []
	for i, yoffset in enumerate(np.linspace(-2.5, 0., 16)):
		filename = case0(folder, str(i)+extension, yoffset=yoffset)
		images.append(imageio.imread(filename))
	images += images[-2::-1] 
	imageio.mimsave(folder + 'bspline1D.gif', images)
	
	def case1(folder, extension): 
		# Set filename
		filename = folder + 'BSplinebasis1D' + extension

		# B-spline properties 
		degree = 3
		knotvector = np.array([0., 0., 0., 0., 0.25, 0.75, 0.75, 1., 1., 1., 1.])

		quadRule   = QuadratureRules(degree, knotvector)
		basis, knots = quadRule.getSampleBasis(sampleSize=201)
		B0 = basis[0].toarray(); B1 = basis[1].toarray()

		fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(6, 4))
		for i in range(np.shape(B0)[0]): 
			ax1.plot(knots, B0[i, :], linewidth=1, color='k')
		ax1.plot([], [], linewidth=1, color='k', label='B-Spline basis')
		ax1.plot(quadRule._uniqueKV, np.zeros(len(quadRule._uniqueKV)), color='k', marker='s', linestyle='', label='Knots')

		ax1.set_xlabel(r'$\xi$')
		ax1.set_xticks(np.linspace(0, 1, 5))
		ax1.set_yticks([0, 0.5, 1])
		ax1.set_ylabel(r'$\hat{b}_{A,\,p}(\xi)$')
		ax1.legend()
		fig.tight_layout()
		fig.savefig(filename, dpi=300)
		return
	
	case1(folder, extension)

elif CASE == 1: # Bivariate functions
	
	def case2(folder, extension):
		# Set filename
		filename = folder + 'BSplinebasis2D' + extension

		# B-Spline properties
		degree1, degree2, nbel = 1, 3, 2
		knotvector1   = createUniformKnotvector_Rmultiplicity(degree1, nbel)
		knotvector2   = createUniformKnotvector_Rmultiplicity(degree2, nbel)
		quadRule1     = QuadratureRules(degree1, knotvector1)
		quadRule2     = QuadratureRules(degree2, knotvector2)
		basis1, knots1 = quadRule1.getSampleBasis()
		basis2, knots2 = quadRule2.getSampleBasis()
		
		# B-Spline 2D
		X, Y = np.meshgrid(np.unique(quadRule1.knotvector), np.unique(quadRule2.knotvector))
		Xb, Yb = np.meshgrid(knots1, knots2)
		Zb = np.kron(basis1[0].todense()[1, :], basis2[0].todense()[2, :]).reshape((len(knots1), len(knots2)))

		from mpl_toolkits.axes_grid1 import make_axes_locatable
		fig, axs = plt.subplots(2, 2, sharex="col", sharey="row", 
								gridspec_kw=dict(height_ratios=[1, 3.2],
												width_ratios=[3.2, 1]), figsize=(5, 5))

		axs[0,1].set_visible(False)
		axs[0,0].set_box_aspect(1/3)
		axs[1,0].set_box_aspect(1)
		axs[1,1].set_box_aspect(3/1)
		axs[1,0].grid(None)
		axs[1,0].pcolormesh(Xb, Yb, Zb.T, cmap='GnBu', shading='gouraud', rasterized=True)
		axs[1,0].plot(X, Y, color='k', marker='s', linestyle='--')
		axs[1,0].plot(X.T, Y.T, color='k', marker='s', linestyle='--')

		axs[1,0].set_yticks([0, 0.5, 1])
		axs[1,0].set_xticks([0, 0.5, 1])

		for i in range(degree1+nbel): 
			axs[0, 0].plot(knots1, np.ravel(basis1[0].todense()[i, :]), linewidth=1, color='k', alpha=0.8)
		uvk = np.unique(quadRule1.knotvector)
		axs[0, 0].plot(uvk, np.zeros(len(uvk)), color='k', marker='s', linestyle='')

		for i in range(degree2+nbel): 
			axs[1, 1].plot(np.ravel(basis2[0].todense()[i, :]), knots2, linewidth=1, color='k', alpha=0.8)
		uvk = np.unique(quadRule2.knotvector)
		axs[1, 1].plot(np.zeros(len(uvk)), uvk, color='k', marker='s', linestyle='')

		axs[0, 0].axis(ymin=0, ymax=1)
		axs[0,0].set_xlabel(r'$\xi_1$')
		axs[0,0].set_ylabel(r'$\hat{b}_{A_1,\,p_1}$')
		axs[1, 1].axis(xmin=0, xmax=1)
		axs[1,1].set_ylabel(r'$\xi_2$')
		axs[1,1].set_xlabel(r'$\hat{b}_{A_2,\,p_2}$')
		axs[1,1].set_xticks([0, 1])
		axs[0,0].set_yticks([0, 1])
		fig.tight_layout()
		fig.savefig(filename, dpi=300) 
		return
	
	case2(folder, extension)

	# Set filename
	filename = folder + 'BSplinesurface'

	# Surface properties
	modelGeo = Geomdl(geoArgs={'name':'quarter_annulus', 'degree': np.array([1, 3, 1]), 
				'nb_refinementByDirection': np.array([1, 1, 1])})
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs={'quadrule':'iga'})
	XoriginalPos = modelPhy.ctrlpts[1, 8]
	YoriginalPos = modelPhy.ctrlpts[0, 8]

	images = []
	for i, offset in enumerate(np.linspace(-0.25, 0.1, 16)):
		modelPhy.ctrlpts[0, 8] = XoriginalPos + offset
		modelPhy.ctrlpts[1, 8] = YoriginalPos + offset
		fig = plot2DGeo(modelPhy)
		newfilename = filename+str(i)+extension
		fig.savefig(newfilename, dpi=300) 
		images.append(imageio.imread(newfilename))

	images += images[-2::-1] 
	imageio.mimsave(folder + 'bspline2D.gif', images)
	