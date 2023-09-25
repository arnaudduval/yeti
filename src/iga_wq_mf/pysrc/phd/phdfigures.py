"""
.. Isogeometric analysis figures
.. Select case: 
.. CASE 0: B-spline curve
.. CASE 1: Univariate B-spline functions in parametric space
.. CASE 2: Bivariate B-spline functions in parametric space
.. CASE 3: Quadrature points in IGA-Galerkin approach
.. CASE 4: Quadrature points in IGA-WQ approach
.. CASE 5: B-spline surface
.. CASE 6: Quadrature rules W00 or W11
.. CASE 7: Plot 3D geometries
.. CASE 8: Plot example of results
"""

from pysrc.lib.__init__ import *
from pysrc.lib.lib_base import (createUniformKnotvector_Rmultiplicity, evalDersBasisPy, cropImage)
from pysrc.lib.lib_quadrules import *
from pysrc.lib.lib_geomdl import Geomdl
from pysrc.lib.lib_part import part

# Select folder
full_path = os.path.realpath(__file__)
folder = os.path.dirname(full_path) + '/results/phdfig/'

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
			ax.plot(x, y, 'k--')

		# In the second direction
		for _ in range(shape[0]):
			x = pts2D[:, _, 0]; y = pts2D[:, _, 1]
			ax.plot(x, y, 'k--')

		return

	ctrlpts = model.ctrlpts
	evalpts = model.interpolateMeshgridField()[0]

	X = np.asarray(evalpts[0, :].reshape((sampleSize, sampleSize)).tolist())
	Y = np.asarray(evalpts[1, :].reshape((sampleSize, sampleSize)).tolist())
	Z = np.zeros(X.shape)

	fig, ax = plt.subplots(figsize=(4.5, 4.5))
	ax.grid(None)
	plot_mesh(ctrlpts, model.nbctrlpts, ax)
	ax.pcolormesh(X, Y, Z, cmap=plt.cm.Pastel1, shading='gouraud')
	ax.plot(ctrlpts[0, :], ctrlpts[1, :], 'o', label='Control\n Points')
	
	ax.set_xticks(np.arange(0, np.ceil(max(evalpts[:, 0]))+2, 1.0))
	ax.set_yticks(np.arange(0, np.ceil(max(evalpts[:, 1]))+1, 1.0))
	ax.set_xlim(left=-0.25, right=2.25)
	ax.set_ylim(bottom=-0.25, top=2.25)

	ax.legend(loc='lower left')
	ax.set_aspect('equal', adjustable='box')
	ax.set_xlabel(r'$x_1$')
	ax.set_ylabel(r'$x_2$')
	fig.tight_layout()
	
	return fig

def plotVerticalLine(x, y, ax=None, color='k'):
	for xi, yi in zip(x, y):
		ax.plot([xi, xi], [0, yi], color=color)
	return

# Set global variables
CASE      = 8
extension = '.png'

if CASE == 0: # B-spline curve

	def case0(folder, extension):
		# Set filename
		filename = folder + 'BsplineCurve' + extension

		# Create the curve 
		degree, nbel = 2, 4
		knotvector   = createUniformKnotvector_Rmultiplicity(degree, nbel)
		crv            = BSpline.Curve()
		crv.degree     = degree
		crv.ctrlpts    = [[-1, 1, 0], [-0.5, 0.25, 0], [0, 2, 0], 
							[0.75, -0.5, 0], [1.5, 1, 0], [2, 0, 0]]
		crv.knotvector = knotvector
		crv.delta      = 0.01

		# Get data
		evalpts = np.asarray(crv.evalpts)
		ctrlpts = np.asarray(crv.ctrlpts)
		x = evalpts[:, 0]; y = evalpts[:, 1]
		
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.plot(x, y, label='B-Spline curve')
		ax.plot(ctrlpts[:, 0], ctrlpts[:, 1], 'o--', markersize=10, label='Control Points')

		ax.legend()
		ax.set_xlabel(r'$x_1$')
		ax.set_ylabel(r'$x_2$')
		ax.axis('equal')
		fig.tight_layout()
		fig.savefig(filename, dpi=300)
		return

	case0(folder, extension)
	
elif CASE == 1: # Univariate functions

	def case1(folder, extension): 
		# Set filename
		filename = folder + 'UnivariateFunctions' + extension

		# B-spline properties 
		degree, nbel = 3, 4
		multiplicity = 1
		knotvector   = createUniformKnotvector_Rmultiplicity(degree, nbel, multiplicity=multiplicity)
		quadRule     = QuadratureRules(degree, knotvector)
		basis, knots = quadRule.getSampleBasis()
		B0 = basis[0].toarray(); B1 = basis[0].toarray()

		fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
		for i in range(np.shape(B0)[0]): 
			ax1.plot(knots, B0[i, :], linewidth=2)
			ax2.plot(knots, B1[i, :], linewidth=2)

		for ax in [ax1, ax2]:
			ax.set_xlabel(r'$\xi$')
			ax.set_xticks(np.linspace(0, 1, nbel+1))
		ax1.set_ylabel(r'$\hat{b}_{A,\,p}(\xi)$')
		ax2.set_ylabel(r"${\hat{b}'}_{A,\,p}(\xi)$")
		fig.tight_layout()
		fig.savefig(filename, dpi=300)
		return
	
	case1(folder, extension)

elif CASE == 2: # Bivariate functions

	def case2(folder, extension, is2D=True):
		# Set filename
		filename = folder + 'BivariateFunctions' + extension

		# B-Spline properties
		degree, nbel = 3, 4
		knotvector   = createUniformKnotvector_Rmultiplicity(degree, nbel)
		quadRule     = QuadratureRules(degree, knotvector)
		basis, knots = quadRule.getSampleBasis()
		B0 = basis[0].toarray()
		B02plot = B0[1, :]

		# B-Spline 2D
		X, Y = np.meshgrid(knots, knots)
		Z = np.kron(B02plot, B02plot).reshape((len(knots), len(knots)))

		if is2D:
			fig, axs = plt.subplots(2, 2, sharex="col", sharey="row", 
									gridspec_kw=dict(height_ratios=[1,3],
													width_ratios=[3,1]))
			axs[0,1].set_visible(False)
			axs[0,0].set_box_aspect(1/3)
			axs[1,0].set_box_aspect(1)
			axs[1,1].set_box_aspect(3/1)
			axs[1,0].grid(None)
			axs[1,0].pcolormesh(X, Y, Z, cmap='GnBu', shading='gouraud', rasterized=True)
			axs[1,0].set_yticks([0, 0.5, 1])

			for i in range(degree+nbel): 
				axs[0, 0].plot(knots, B0[i, :], color="0.8")
				axs[1, 1].plot(B0[i, :], knots, color="0.8")

			axs[0,0].plot(knots, B02plot); axs[0, 0].axis(ymin=0, ymax=1)
			axs[0,0].set_xlabel(r'$\xi_1$')
			axs[0,0].set_ylabel(r'$\hat{b}_{A_1,\,p_1}$')
			axs[1,1].plot(B02plot, knots); axs[1, 1].axis(xmin=0, xmax=1)
			axs[1,1].set_ylabel(r'$\xi_2$')
			axs[1,1].set_xlabel(r'$\hat{b}_{A_2,\,p_2}$')
			axs[1,1].set_xticks([0, 1])
			fig.tight_layout()
			fig.savefig(filename, dpi=300) 

		else:
			fig = plt.figure()
			ax = plt.axes(projection='3d')
			ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='GnBu', edgecolor='none')
			ax.grid(False)
			ax.set_xticks([0, 0.5, 1])
			ax.set_yticks([0, 0.5, 1])
			ax.set_zticks([0, 0.1, 0.2])
			ax.set_xlabel(r'$\xi_1$')
			ax.set_ylabel(r'$\xi_2$')
			ax.set_zlabel(r'$B_{i,p}(\xi)$')
			fig.tight_layout()
			fig.savefig(filename, dpi=300) 
		return
	
	case2(folder, extension)

elif CASE == 3: # Quadrature points in IGA

	def case3(folder, extension):
		# Set filename
		filename = folder + 'QuadPtsIGA' + extension

		fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
		degree = 4
		for ax, nbel in zip([ax1, ax2], [8, 32]):
			knotvector  = createUniformKnotvector_Rmultiplicity(degree, nbel)
			quadRule    = GaussQuadrature(degree, knotvector, quadArgs={})
			quadRule.getQuadratureRulesInfo()
			XX, YY      = np.meshgrid(quadRule.quadPtsPos, quadRule.quadPtsPos)
			ax.plot(XX, YY, 'ko', markersize=1.2)

			grid = np.linspace(0, 1, nbel+1)
			for i in grid:
				ax.plot([i, i], [0, 1], 'grey', linewidth=0.5, alpha=0.8)
				ax.plot([0, 1], [i, i], 'grey', linewidth=0.5, alpha=0.8)

			ax.set_xticks([0, 0.5, 1])
			ax.set_yticks([0, 0.5, 1])
			ax.axis('equal')
			ax.set_ylabel(r'$\xi_2$')
			ax.set_xlabel(r'$\xi_1$')   

		fig.tight_layout()
		fig.savefig(filename, dpi=300) 
		return

	case3(folder, extension)

elif CASE == 4: # Quadrature points in WQ

	def case4(folder, extension):
		# Set filename
		filename = folder + 'QuadPtsWQ' + extension

		fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))
		degree = 4
		for ax, nbel in zip([ax1, ax2], [8, 32]):
			knotvector  = createUniformKnotvector_Rmultiplicity(degree, nbel)
			quadRule    = WeightedQuadrature(degree, knotvector, quadArgs={})
			quadRule.getQuadratureRulesInfo()
			XX, YY      = np.meshgrid(quadRule.quadPtsPos, quadRule.quadPtsPos)
			ax.plot(XX, YY, 'ko', markersize=1.2)

			grid = np.linspace(0.,1,nbel+1)
			for i in grid:
				ax.plot([i, i], [0, 1], 'grey', linewidth=0.5, alpha=0.8)
				ax.plot([0, 1], [i, i], 'grey', linewidth=0.5, alpha=0.8)

			ax.set_xticks([0, 0.5, 1])
			ax.set_yticks([0, 0.5, 1])
			ax.axis('equal')
			ax.set_ylabel(r'$\xi_2$')
			ax.set_xlabel(r'$\xi_1$')

		fig.tight_layout()
		fig.savefig(filename, dpi=300) 
		return
	
	case4(folder, extension)

elif CASE == 5: # B-spline surface

	# Set filename
	filename = folder + 'BsplineSurface' + extension

	# Surface properties
	modelGeo = Geomdl(geoArgs={'name':'quarter_annulus'})
	modelIGA = modelGeo.getIGAParametrization()
	modelPhy = part(modelIGA, quadArgs={'quadrule':'iga'})
	fig = plot2DGeo(modelPhy)
	fig.savefig(filename, dpi=300) 

elif CASE == 6: # Weights W00 and W11

	def case6(folder, extension):
		WeightName = 0 # or 0
		if WeightName not in [0, 1]: raise Warning('Not possible')
		if WeightName == 0: ylim1 = [0, 1];  ylim2 = [0, 0.25]
		else:               ylim1 = [-3, 3]; ylim2 = [-0.6, 0.6]

		# Set filename
		filename = folder + 'WeightsW' + str(WeightName) + extension

		# B-spline properties 
		WeightPos    = 2
		degree, nbel = 2, 3
		knotvector   = createUniformKnotvector_Rmultiplicity(degree, nbel)
		quadRule     = WeightedQuadrature(degree, knotvector, {'type': 1, 'extra':{'r': 3, 's': 2}})
		quadRule.getQuadratureRulesInfo()
		basis, knots = quadRule.getSampleBasis()
		B = basis[WeightName].toarray()

		# Get weights
		quadRule.getDenseQuadRules()
		W00 = quadRule._denseWeights[-WeightName]
		wgt = W00.toarray()[WeightPos, :]

		if WeightName == 0:
			Bref = B
		else:
			kvref = createUniformKnotvector_Rmultiplicity(degree-1, nbel)
			Bref = evalDersBasisPy(degree-1, kvref, knots)[0]
			Bref = Bref.toarray()

		fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8,4))
		for i in range(np.shape(Bref)[0]): 
			ax1.plot(knots, Bref[i, :], linewidth=2)

		# Fill basis chosen
		ax1.fill_between(x=knots, y1=B[WeightPos, :], color="g", alpha=0.2)
		ax1.set_ylim(ylim1)

		ax2 = ax1.twinx()
		ax2.plot(quadRule.quadPtsPos, wgt, 'ko')
		plotVerticalLine(quadRule.quadPtsPos, wgt, ax2)
		ax2.set_ylim(ylim2)
		ax2.grid(None)

		ax1.set_xlabel(r'$\xi$')
		ax1.set_ylabel('Basis')
		ax2.set_ylabel('Weights')
		ax1.set_xticks(np.linspace(0, 1, nbel+1), ['0', '1/3', '2/3', '1'])
		fig.tight_layout()
		fig.savefig(filename, dpi=300)
		return
	
	case6(folder, extension)

elif CASE == 7: # 2D Geometries

	def case7(folder):
		# Create model
		name    = 'QA'
		filename = folder + 'VTK_' + name + '.png'
		degree, cuts = 5, 5
		geoArgs   = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
		quadArgs  = {'quadrule': 'wq', 'type': 1}

		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)
		modelPhy.exportResultsCP(folder=folder)
		
		# Read data
		fileVTK   = folder + 'IGAparametrization'
		grid      = pv.read(fileVTK + '.vts')
		
		sargs = dict(
				title = 'normalized det J',
				title_font_size=50,
				label_font_size=40,
				shadow=True,
				n_labels=2,
				fmt="%.1f",
				position_x=0.2, 
				position_y=0.1,
		)
		pv.start_xvfb()
		plotter = pv.Plotter(off_screen=True)
		if name == 'SQ':
			boring_cmap = plt.cm.get_cmap("GnBu", 1)
			plotter.add_mesh(grid, cmap=boring_cmap, scalar_bar_args=sargs)
		else:
			plotter.add_mesh(grid, cmap='viridis', scalar_bar_args=sargs)
		
		plotter.camera_position  = 'xy'
		plotter.camera.zoom(0.8)

		plotter.background_color = 'white'
		plotter.window_size = [1600, 1600]
		plotter.screenshot(filename)
		os.remove(fileVTK + '.vts')
		cropImage(filename)

	case7(folder)

elif CASE == 8: # 3D Geometries

	def case8(folder):
		# Create model
		name    = 'TR'
		filename = folder + 'VTK_' + name + '.png'
		degree, cuts = 5, 5
		geoArgs   = {'name': name, 'degree': degree*np.ones(3, dtype=int), 
					'nb_refinementByDirection': cuts*np.ones(3, dtype=int)}
		quadArgs  = {'quadrule': 'wq', 'type': 1}

		modelGeo = Geomdl(geoArgs)
		modelIGA = modelGeo.getIGAParametrization()
		modelPhy = part(modelIGA, quadArgs=quadArgs)
		modelPhy.exportResultsCP(folder=folder)

		# Read data
		fileVTK   = folder + 'IGAparametrization'
		grid      = pv.read(fileVTK + '.vts')
	
		sargs = dict(
				title = 'normalized det J',
				title_font_size=50,
				label_font_size=40,
				shadow=True,
				n_labels=2,
				fmt="%.1f",
				position_x=0.2, 
				position_y=0.1,
		)
		pv.start_xvfb()
		plotter = pv.Plotter(off_screen=True)
		if name == 'CB': 
			boring_cmap = plt.cm.get_cmap("GnBu", 1)
			plotter.add_mesh(grid, cmap=boring_cmap, scalar_bar_args=sargs)
			plotter.camera.zoom(0.6)
		else: 
			plotter.add_mesh(grid, cmap='viridis', scalar_bar_args=sargs)

		if name == 'VB': 
			plotter.camera_position  = 'yz'
			plotter.camera.elevation = 45
		elif name == 'RQA': 
			plotter.camera_position  = 'xz'
			plotter.camera.azimuth   = -45
			plotter.camera.zoom(0.8)

		plotter.background_color = 'white'
		plotter.window_size = [1600, 1600]
		plotter.screenshot(filename)
		os.remove(fileVTK + '.vts')
		cropImage(filename)
		return
	
	case8(folder)

else: raise Warning('Case unkwnon')